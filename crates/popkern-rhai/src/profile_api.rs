//! Planes, analytic profiles, extrude/revolve ends, and axes (spec
//! §§1.4, 3.1, 5). A profile is a value the script builds up and
//! records with `.sketch(label)`.

use popkern_script::{AxisSpec, EdgeSelector, End, FaceSelector, Plane, Scalar};
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args;
use crate::body_api::ScriptBody;
use crate::error::{IntoScriptResult, argument_error};
use crate::model_api::ScriptModel;

/// An extrude or revolve end condition value.
#[derive(Debug, Clone)]
pub enum ScriptEnd {
    /// `blind(d)`, `through_all()`, `up_to_surface(sel)`,
    /// `offset_from_surface(sel, d)`, `mid_plane(d)`.
    Extrude(End),
    /// `full_circle()` (revolve only).
    FullCircle,
}

/// A revolve/rotation axis value.
#[derive(Debug, Clone)]
pub enum ScriptAxis {
    /// `axis_at(origin, direction)`.
    Literal {
        /// A point on the axis.
        origin: [f64; 3],
        /// The axis direction.
        direction: [f64; 3],
    },
    /// `axis(edge_selection)`: resolved against the active body at
    /// record time.
    Edge(EdgeSelector),
}

impl ScriptAxis {
    /// The recorder's axis argument.
    pub fn spec(&self) -> AxisSpec {
        match self {
            ScriptAxis::Literal { origin, direction } => AxisSpec::Literal {
                origin: *origin,
                direction: *direction,
            },
            ScriptAxis::Edge(selection) => AxisSpec::Edge(selection.clone()),
        }
    }
}

#[derive(Debug, Clone)]
enum Shape {
    Rect(Scalar, Scalar),
    Circle(Scalar),
    Polygon(Vec<(f64, f64)>),
    Path {
        start: (Scalar, Scalar),
        segments: Vec<Segment>,
    },
}

#[derive(Debug, Clone)]
enum Segment {
    Line(Scalar, Scalar),
    Arc((f64, f64), (f64, f64)),
}

/// A profile under construction (`model.profile(plane)`).
#[derive(Debug, Clone)]
pub struct ScriptProfile {
    model: ScriptModel,
    plane: Plane,
    shape: Option<Shape>,
    center: (Scalar, Scalar),
}

impl ScriptProfile {
    pub(crate) fn new(model: ScriptModel, plane: Plane) -> ScriptProfile {
        ScriptProfile {
            model,
            plane,
            shape: None,
            center: (Scalar::Literal(0.0), Scalar::Literal(0.0)),
        }
    }

    fn with_shape(&self, shape: Shape) -> ScriptProfile {
        let mut next = self.clone();
        next.shape = Some(shape);
        next
    }

    fn push_segment(&self, segment: Segment) -> Result<ScriptProfile, Box<EvalAltResult>> {
        let mut next = self.clone();
        match &mut next.shape {
            Some(Shape::Path { segments, .. }) => {
                segments.push(segment);
                Ok(next)
            }
            _ => Err(argument_error(
                "path",
                "line_to/arc_to need a path started with .path(x, y)",
            )),
        }
    }

    fn sketch(&self, label: &str) -> Result<ScriptBody, Box<EvalAltResult>> {
        let shape = self.shape.clone().ok_or_else(|| {
            argument_error("profile", "no shape given (rect/circle/polygon/path)")
        })?;
        let model = self.model.clone();
        let body = {
            let mut guard = model.borrow_mut();
            let builder = guard.profile(self.plane.clone());
            let builder = match shape {
                Shape::Rect(w, h) => builder.rect(w, h),
                Shape::Circle(r) => builder.circle(r),
                Shape::Polygon(points) => builder.polygon(&points),
                Shape::Path { start, segments } => {
                    let mut b = builder.path(start.0, start.1);
                    for segment in segments {
                        b = match segment {
                            Segment::Line(x, y) => b.line_to(x, y),
                            Segment::Arc(via, to) => b.arc_to(via, to),
                        };
                    }
                    b
                }
            };
            builder
                .at(self.center.0.clone(), self.center.1.clone())
                .sketch(label)
                .script()?
        };
        Ok(ScriptBody::new(model, body))
    }
}

fn polygon_points(value: &Dynamic) -> Result<Vec<(f64, f64)>, Box<EvalAltResult>> {
    let items = value
        .clone()
        .into_array()
        .map_err(|_| argument_error("points", "expected an array of [x, y] pairs"))?;
    items
        .iter()
        .map(|p| args::point("points", p, 2).map(|v| (v[0], v[1])))
        .collect()
}

fn pair(argument: &str, value: &Dynamic) -> Result<(f64, f64), Box<EvalAltResult>> {
    let v = args::point(argument, value, 2)?;
    Ok((v[0], v[1]))
}

/// Register planes, profiles, ends, and axes.
pub fn register(engine: &mut Engine) {
    engine
        .register_type_with_name::<Plane>("Plane")
        .register_fn("xy", || Plane::XY)
        .register_fn("xz", || Plane::XZ)
        .register_fn("yz", || Plane::YZ)
        .register_fn("xy_at", |z: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(Plane::xy_at(args::scalar("z", &z)?))
        })
        .register_fn("xz_at", |y: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(Plane::xz_at(args::scalar("y", &y)?))
        })
        .register_fn("yz_at", |x: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(Plane::yz_at(args::scalar("x", &x)?))
        })
        .register_fn("on_face", Plane::on_face)
        .register_fn("on_face_offset", |sel: FaceSelector, d: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(Plane::on_face_offset(sel, args::scalar("distance", &d)?))
        })
        .register_fn("plane_at", |origin: Dynamic, normal: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(Plane::custom(
                args::vec3("origin", &origin)?,
                args::vec3("normal", &normal)?,
            ))
        });
    engine
        .register_type_with_name::<ScriptProfile>("Profile")
        .register_fn("rect", |p: &mut ScriptProfile, w: Dynamic, h: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(p.with_shape(Shape::Rect(
                args::scalar("width", &w)?,
                args::scalar("height", &h)?,
            )))
        })
        .register_fn("circle", |p: &mut ScriptProfile, r: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(p.with_shape(Shape::Circle(args::scalar("radius", &r)?)))
        })
        .register_fn("polygon", |p: &mut ScriptProfile, points: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(p.with_shape(Shape::Polygon(polygon_points(&points)?)))
        })
        .register_fn("path", |p: &mut ScriptProfile, x: Dynamic, y: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(p.with_shape(Shape::Path {
                start: (args::scalar("x", &x)?, args::scalar("y", &y)?),
                segments: Vec::new(),
            }))
        })
        .register_fn(
            "line_to",
            |p: &mut ScriptProfile, x: Dynamic, y: Dynamic| {
                p.push_segment(Segment::Line(
                    args::scalar("x", &x)?,
                    args::scalar("y", &y)?,
                ))
            },
        )
        .register_fn(
            "arc_to",
            |p: &mut ScriptProfile, via: Dynamic, to: Dynamic| {
                p.push_segment(Segment::Arc(pair("via", &via)?, pair("to", &to)?))
            },
        )
        .register_fn("at", |p: &mut ScriptProfile, x: Dynamic, y: Dynamic| {
            let mut next = p.clone();
            next.center = (args::scalar("x", &x)?, args::scalar("y", &y)?);
            Ok::<_, Box<EvalAltResult>>(next)
        })
        .register_fn("sketch", |p: &mut ScriptProfile| p.sketch(""))
        .register_fn("sketch", |p: &mut ScriptProfile, label: &str| {
            p.sketch(label)
        });
    engine
        .register_type_with_name::<ScriptEnd>("End")
        .register_fn("blind", |d: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(ScriptEnd::Extrude(End::Blind(args::scalar(
                "distance", &d,
            )?)))
        })
        .register_fn("through_all", || ScriptEnd::Extrude(End::ThroughAll))
        .register_fn("up_to_surface", |sel: FaceSelector| {
            ScriptEnd::Extrude(End::UpToSurface(sel))
        })
        .register_fn("offset_from_surface", |sel: FaceSelector, d: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(ScriptEnd::Extrude(End::OffsetFromSurface(
                sel,
                args::scalar("distance", &d)?,
            )))
        })
        .register_fn("mid_plane", |d: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(ScriptEnd::Extrude(End::MidPlane(args::scalar(
                "distance", &d,
            )?)))
        })
        .register_fn("full_circle", || ScriptEnd::FullCircle);
    engine
        .register_type_with_name::<ScriptAxis>("Axis")
        .register_fn("axis_at", |origin: Dynamic, direction: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(ScriptAxis::Literal {
                origin: args::vec3("origin", &origin)?,
                direction: args::vec3("direction", &direction)?,
            })
        })
        .register_fn("axis", ScriptAxis::Edge);
}
