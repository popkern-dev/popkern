//! Analytic profile builder (spec §3.1): records `Feature::Sketch`
//! with the analytic arm; every dimension is a scalar binding.

use popkern_history::feature::{Feature, SketchPlaneDesc};

use crate::plane::Plane;
use crate::selection::{exactly_one, resolve_faces};
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
    Polygon(Vec<[Scalar; 2]>),
    Path {
        start: [Scalar; 2],
        segments: Vec<Segment>,
    },
}

enum Segment {
    Line([Scalar; 2]),
    Arc([Scalar; 2], [Scalar; 2]),
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
            reason: "no shape given (rect/circle/polygon/path)".into(),
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
            Shape::Polygon(pts) => ProfileShape::Polygon {
                points: pts
                    .iter()
                    .enumerate()
                    .map(|(i, [x, y])| {
                        Ok([bind(&format!("p{i}_x"), x)?, bind(&format!("p{i}_y"), y)?])
                    })
                    .collect::<Result<Vec<_>, ScriptError>>()?,
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
                    })
                    .collect::<Result<Vec<_>, ScriptError>>()?,
            },
        };
        let center = [bind("center_x", &center[0])?, bind("center_y", &center[1])?];
        let profile = AnalyticProfile { shape, center };
        let plane = match plane.simple_desc() {
            Some(desc) => desc,
            None => resolve_plane(model, &plane)?,
        };
        let feature = Feature::Sketch {
            plane,
            sketch_data: SketchPayload::Analytic(profile),
        };
        let id = record_feature(&mut model.tree, model.tolerance, feature, label)?;
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
fn resolve_plane(model: &mut Model, plane: &Plane) -> Result<SketchPlaneDesc, ScriptError> {
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
            let labels = |l: &str| model.entity_by_label(l);
            let face = exactly_one(
                resolve_faces(&brep, sel, &labels)?,
                crate::error::SelectionKind::Faces,
            )?;
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
        _ => plane.simple_desc().ok_or_else(|| ScriptError::Argument {
            argument: "plane".into(),
            reason: "unsupported plane".into(),
        }),
    }
}
