//! Analytic profile builder (spec §3.1): records `Feature::Sketch`
//! with the analytic arm; every dimension is a scalar binding.

use popkern_history::feature::{Feature, SketchPlaneDesc};

use crate::plane::Plane;
use crate::selection::{EntityRef, exactly_one, resolve_faces};
use popkern_history::feature_sketch::{ScalarBinding, SketchPayload};
use popkern_sketch::profile::{AnalyticProfile, PathSegment, ProfileShape};

use crate::body::Body;
use crate::error::ScriptError;
use crate::model::Model;
use crate::record::record_feature;
use crate::scalar::{Scalar, bind_scalar};

/// A profile under construction.
pub struct ProfileBuilder<'m> {
    model: &'m mut Model,
    plane: Plane,
    shape: Option<Shape>,
    center: [Scalar; 2],
}

enum Shape {
    Rect(Scalar, Scalar),
    Circle(Scalar),
    Ellipse(Scalar, Scalar),
    Polygon(Vec<[Scalar; 2]>),
    Spline(Vec<[Scalar; 2]>, SplineOptions),
    Path {
        start: [Scalar; 2],
        segments: Vec<Segment>,
    },
}

enum Segment {
    Line([Scalar; 2]),
    Arc([Scalar; 2], [Scalar; 2]),
    Conic([Scalar; 2], [Scalar; 2], f64),
    Spline(Vec<[Scalar; 2]>, [Scalar; 2], SplineOptions),
}

/// A spline's shape beyond its control points (spec §3.1): degree,
/// rational weights and, for a path segment, a clamped knot vector.
/// Weights and knots are dimensionless.
#[derive(Debug, Clone, PartialEq)]
pub struct SplineOptions {
    /// Degree (default 3).
    pub degree: usize,
    /// One weight per control point (empty: all 1).
    pub weights: Vec<f64>,
    /// Clamped knot vector (empty: clamped uniform). Path segments
    /// only: a closed spline's knots are uniform.
    pub knots: Vec<f64>,
}

impl Default for SplineOptions {
    fn default() -> Self {
        SplineOptions {
            degree: 3,
            weights: Vec::new(),
            knots: Vec::new(),
        }
    }
}

fn literal(p: &(f64, f64)) -> [Scalar; 2] {
    [Scalar::Literal(p.0), Scalar::Literal(p.1)]
}

impl<'m> ProfileBuilder<'m> {
    pub(crate) fn new(model: &'m mut Model, plane: Plane) -> Self {
        ProfileBuilder {
            model,
            plane,
            shape: None,
            center: [Scalar::Literal(0.0), Scalar::Literal(0.0)],
        }
    }
    /// A rectangle centred on the placement.
    pub fn rect(mut self, width: impl Into<Scalar>, height: impl Into<Scalar>) -> Self {
        self.shape = Some(Shape::Rect(width.into(), height.into()));
        self
    }
    /// A circle centred on the placement.
    pub fn circle(mut self, radius: impl Into<Scalar>) -> Self {
        self.shape = Some(Shape::Circle(radius.into()));
        self
    }
    /// An axis-aligned ellipse centred on the placement: semi-axis
    /// `rx` along the plane's u, `ry` along its v.
    pub fn ellipse(mut self, rx: impl Into<Scalar>, ry: impl Into<Scalar>) -> Self {
        self.shape = Some(Shape::Ellipse(rx.into(), ry.into()));
        self
    }
    /// A closed (periodic) B-spline by its control polygon (literal
    /// coordinates, offset by the placement).
    pub fn spline(mut self, points: &[(f64, f64)], options: SplineOptions) -> Self {
        self.shape = Some(Shape::Spline(points.iter().map(literal).collect(), options));
        self
    }
    /// A closed polygon (literal coordinates).
    pub fn polygon(mut self, points: &[(f64, f64)]) -> Self {
        self.shape = Some(Shape::Polygon(
            points
                .iter()
                .map(|(x, y)| [Scalar::Literal(*x), Scalar::Literal(*y)])
                .collect(),
        ));
        self
    }
    /// Start a line/arc path at `(x, y)`.
    pub fn path(mut self, x: impl Into<Scalar>, y: impl Into<Scalar>) -> Self {
        self.shape = Some(Shape::Path {
            start: [x.into(), y.into()],
            segments: Vec::new(),
        });
        self
    }
    /// Straight segment to `(x, y)`.
    pub fn line_to(mut self, x: impl Into<Scalar>, y: impl Into<Scalar>) -> Self {
        if let Some(Shape::Path { segments, .. }) = &mut self.shape {
            segments.push(Segment::Line([x.into(), y.into()]));
        }
        self
    }
    /// Arc through `via` to `to`.
    pub fn arc_to(mut self, via: (f64, f64), to: (f64, f64)) -> Self {
        if let Some(Shape::Path { segments, .. }) = &mut self.shape {
            segments.push(Segment::Arc(
                [Scalar::Literal(via.0), Scalar::Literal(via.1)],
                [Scalar::Literal(to.0), Scalar::Literal(to.1)],
            ));
        }
        self
    }
    /// Conic arc to `to` whose end tangents meet at `control`, with
    /// middle weight `weight` (below 1 an ellipse arc, 1 a parabola,
    /// above 1 a hyperbola).
    pub fn conic_to(mut self, control: (f64, f64), to: (f64, f64), weight: f64) -> Self {
        if let Some(Shape::Path { segments, .. }) = &mut self.shape {
            segments.push(Segment::Conic(literal(&control), literal(&to), weight));
        }
        self
    }
    /// Spline segment through the control points `points` (the last
    /// is the end point), starting at the current point. An empty
    /// `points` adds nothing.
    pub fn spline_to(mut self, points: &[(f64, f64)], options: SplineOptions) -> Self {
        if let Some(Shape::Path { segments, .. }) = &mut self.shape
            && let Some((to, control)) = points.split_last()
        {
            segments.push(Segment::Spline(
                control.iter().map(literal).collect(),
                literal(to),
                options,
            ));
        }
        self
    }
    /// Close the path (spec §3.1 `.close()`). Paths are always closed
    /// back to their start when compiled, so this records nothing; it
    /// exists so the canonical spelling reads as written.
    pub fn close(self) -> Self {
        self
    }
    /// Place the shape at `(x, y)` on the plane.
    pub fn at(mut self, x: impl Into<Scalar>, y: impl Into<Scalar>) -> Self {
        self.center = [x.into(), y.into()];
        self
    }

    /// Record the sketch feature; returns the body it starts (or
    /// belongs to, for on-face planes).
    pub fn sketch(self, label: &str) -> Result<Body, ScriptError> {
        let ProfileBuilder {
            model,
            plane,
            shape,
            center,
        } = self;
        let shape = shape.ok_or_else(|| ScriptError::Argument {
            argument: "profile".into(),
            reason: "no shape given (rect/circle/ellipse/polygon/spline/path)".into(),
        })?;
        let unit = model.unit;
        let mut bind = |role: &str, s: &Scalar| -> Result<ScalarBinding, ScriptError> {
            Ok(ScalarBinding::Parameter(bind_scalar(
                &mut model.tree,
                &mut model.pending_counter,
                role,
                s,
                unit,
            )?))
        };
        let shape = match shape {
            Shape::Rect(w, h) => ProfileShape::Rect {
                width: bind("width", &w)?,
                height: bind("height", &h)?,
            },
            Shape::Circle(r) => ProfileShape::Circle {
                radius: bind("radius", &r)?,
            },
            Shape::Ellipse(rx, ry) => ProfileShape::Ellipse {
                rx: bind("rx", &rx)?,
                ry: bind("ry", &ry)?,
            },
            Shape::Polygon(pts) => ProfileShape::Polygon {
                points: bind_points(&mut bind, "p", &pts)?,
            },
            Shape::Spline(pts, options) => ProfileShape::Spline {
                points: bind_points(&mut bind, "p", &pts)?,
                degree: options.degree,
                weights: options.weights,
            },
            Shape::Path { start, segments } => ProfileShape::Path {
                start: [bind("start_x", &start[0])?, bind("start_y", &start[1])?],
                segments: segments
                    .iter()
                    .enumerate()
                    .map(|(i, seg)| match seg {
                        Segment::Line(to) => Ok(PathSegment::Line {
                            to: [
                                bind(&format!("s{i}_x"), &to[0])?,
                                bind(&format!("s{i}_y"), &to[1])?,
                            ],
                        }),
                        Segment::Arc(via, to) => Ok(PathSegment::Arc {
                            via: [
                                bind(&format!("s{i}_via_x"), &via[0])?,
                                bind(&format!("s{i}_via_y"), &via[1])?,
                            ],
                            to: [
                                bind(&format!("s{i}_x"), &to[0])?,
                                bind(&format!("s{i}_y"), &to[1])?,
                            ],
                        }),
                        Segment::Conic(control, to, weight) => Ok(PathSegment::Conic {
                            control: [
                                bind(&format!("s{i}_ctl_x"), &control[0])?,
                                bind(&format!("s{i}_ctl_y"), &control[1])?,
                            ],
                            to: [
                                bind(&format!("s{i}_x"), &to[0])?,
                                bind(&format!("s{i}_y"), &to[1])?,
                            ],
                            weight: *weight,
                        }),
                        Segment::Spline(control, to, options) => Ok(PathSegment::Spline {
                            control: bind_points(&mut bind, &format!("s{i}_c"), control)?,
                            to: [
                                bind(&format!("s{i}_x"), &to[0])?,
                                bind(&format!("s{i}_y"), &to[1])?,
                            ],
                            degree: options.degree,
                            weights: options.weights.clone(),
                            knots: options.knots.clone(),
                        }),
                    })
                    .collect::<Result<Vec<_>, ScriptError>>()?,
            },
        };
        let center = [bind("center_x", &center[0])?, bind("center_y", &center[1])?];
        let profile = AnalyticProfile { shape, center };
        let mut records = Vec::new();
        let plane = match plane.simple_desc(unit) {
            Some(desc) => desc,
            None => resolve_plane(model, &plane, &mut records)?,
        };
        let feature = Feature::Sketch {
            plane,
            sketch_data: SketchPayload::Analytic(profile),
        };
        let id = record_feature(
            &mut model.tree,
            model.tolerance,
            feature,
            label,
            None,
            records,
            &mut model.geometry_current,
        )?;
        model.note_source(id);
        let body = model
            .tree
            .feature(id)
            .map(|n| n.body)
            .ok_or(popkern_history::HistoryError::FeatureNotFound)?;
        model.active = Some(body);
        Ok(Body::new(body))
    }
}

/// Resolve the planes that need a parameter or a face of the active body.
fn resolve_plane(
    model: &mut Model,
    plane: &Plane,
    records: &mut Vec<popkern_history::SelectionRecord>,
) -> Result<SketchPlaneDesc, ScriptError> {
    let unit = model.unit;
    match plane {
        Plane::Offset(base, d) => Ok(SketchPlaneDesc::BaseOffset {
            base: *base,
            distance: bind_scalar(
                &mut model.tree,
                &mut model.pending_counter,
                "plane_offset",
                d,
                unit,
            )?,
        }),
        Plane::OnFace(sel) | Plane::OnFaceOffset(sel, _) => {
            let body = model.active()?;
            let brep = model.body_brep(body.id)?;
            let labels = |r: EntityRef| model.resolve_entity_ref(r);
            let snapshot = resolve_faces(&brep, sel, &labels, unit)?;
            records.push(Body::selection_record("plane", &snapshot));
            let face = exactly_one(snapshot, crate::error::SelectionKind::Faces)?;
            match plane {
                Plane::OnFaceOffset(_, d) => Ok(SketchPlaneDesc::OnFaceOffset(
                    face,
                    bind_scalar(
                        &mut model.tree,
                        &mut model.pending_counter,
                        "offset",
                        d,
                        unit,
                    )?,
                )),
                _ => Ok(SketchPlaneDesc::OnFace(face)),
            }
        }
        _ => plane
            .simple_desc(unit)
            .ok_or_else(|| ScriptError::Argument {
                argument: "plane".into(),
                reason: "unsupported plane".into(),
            }),
    }
}

/// Bind a run of points as `{prefix}{i}_x` / `{prefix}{i}_y`.
fn bind_points(
    bind: &mut impl FnMut(&str, &Scalar) -> Result<ScalarBinding, ScriptError>,
    prefix: &str,
    points: &[[Scalar; 2]],
) -> Result<Vec<[ScalarBinding; 2]>, ScriptError> {
    points
        .iter()
        .enumerate()
        .map(|(i, [x, y])| {
            Ok([
                bind(&format!("{prefix}{i}_x"), x)?,
                bind(&format!("{prefix}{i}_y"), y)?,
            ])
        })
        .collect()
}
