//! Analytic sketch profiles: closed 2-D outlines described by shape
//! and placement, compiled to EXACT geometry without a constraint
//! solve.
//!
//! This is the second sketch dialect (scripted-modeling design D2):
//! the constraint dialect is a [`Sketch`](crate::Sketch) re-solved on
//! replay; the analytic dialect is an [`AnalyticProfile`] compiled
//! directly. Exactness matters — the sweep kernel's canonical
//! recognition turns an extruded rational-quadratic circle into a true
//! `Cylinder` (and lines into `Plane`s), and the extrude operator
//! promotes an exact `Arc` edge the same way; an approximated profile
//! would silently degrade every downstream operation (booleans,
//! fillets, STEP) to NURBS. Everything here compiles to the exact
//! forms those recognizers accept (see [`profile_compile`]).
//!
//! The description is generic over its scalar type `S` so a document
//! can bind dimensions to parameters (`S` = a parameter binding) and
//! resolve them to `f64` at replay through
//! [`AnalyticProfile::map_scalars`]; the compiler works on
//! `AnalyticProfile<f64>`.
//!
//! [`profile_compile`]: crate::profile_compile

use serde::{Deserialize, Serialize};

/// A closed profile shape in sketch coordinates. Every shape is
/// centred on the profile's placement (see [`AnalyticProfile::center`]).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ProfileShape<S> {
    /// Axis-aligned rectangle: `width` along X, `height` along Y,
    /// centred on the placement.
    Rect {
        /// Extent along sketch X.
        width: S,
        /// Extent along sketch Y.
        height: S,
    },
    /// Circle of the given radius, centred on the placement.
    Circle {
        /// Radius.
        radius: S,
    },
    /// Closed simple polygon (first point NOT repeated); winding may
    /// be either direction.
    Polygon {
        /// Vertices in order.
        points: Vec<[S; 2]>,
    },
    /// Mixed line/arc outline. An outline whose last segment does not
    /// end at `start` is closed by a straight segment back to it.
    Path {
        /// Outline start point.
        start: [S; 2],
        /// Segments in order.
        segments: Vec<PathSegment<S>>,
    },
}

/// One segment of a [`ProfileShape::Path`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum PathSegment<S> {
    /// Straight segment to `to`.
    Line {
        /// End point.
        to: [S; 2],
    },
    /// Circular arc from the current point THROUGH `via` to `to`
    /// (three-point form: unambiguous, no sweep-direction flags).
    Arc {
        /// A point on the arc strictly between its ends.
        via: [S; 2],
        /// End point.
        to: [S; 2],
    },
}

/// A closed analytic profile: one shape and its in-plane placement.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AnalyticProfile<S> {
    /// The outline.
    pub shape: ProfileShape<S>,
    /// In-plane placement of the shape's centre (sketch coordinates).
    pub center: [S; 2],
}

/// Why an analytic profile cannot be compiled.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ProfileError {
    /// A polygon needs at least three vertices.
    #[error("polygon profile needs at least 3 points, got {got}")]
    PolygonTooShort {
        /// Vertices supplied.
        got: usize,
    },
    /// A path needs at least two segments to enclose area.
    #[error("path profile needs at least 2 segments")]
    PathTooShort,
    /// The three points of an arc segment are collinear (no circle
    /// passes through them).
    #[error("arc_to points are collinear")]
    CollinearArc,
    /// A dimension that must be strictly positive is not.
    #[error("{what} must be positive and finite, got {value}")]
    NonPositiveDimension {
        /// Which dimension.
        what: &'static str,
        /// The offending value.
        value: f64,
    },
}

impl<S> ProfileShape<S> {
    /// Every scalar in the shape, in a stable order.
    pub fn scalars(&self) -> Vec<&S> {
        match self {
            ProfileShape::Rect { width, height } => vec![width, height],
            ProfileShape::Circle { radius } => vec![radius],
            ProfileShape::Polygon { points } => points.iter().flatten().collect(),
            ProfileShape::Path { start, segments } => {
                let mut out: Vec<&S> = start.iter().collect();
                for seg in segments {
                    match seg {
                        PathSegment::Line { to } => out.extend(to.iter()),
                        PathSegment::Arc { via, to } => {
                            out.extend(via.iter());
                            out.extend(to.iter());
                        }
                    }
                }
                out
            }
        }
    }

    /// Rebuild the shape with every scalar mapped through `f`,
    /// stopping at the first error.
    pub fn map_scalars<T, E>(
        &self,
        mut f: impl FnMut(&S) -> Result<T, E>,
    ) -> Result<ProfileShape<T>, E> {
        let mut pt = |p: &[S; 2]| -> Result<[T; 2], E> { Ok([f(&p[0])?, f(&p[1])?]) };
        Ok(match self {
            ProfileShape::Rect { width, height } => ProfileShape::Rect {
                width: f(width)?,
                height: f(height)?,
            },
            ProfileShape::Circle { radius } => ProfileShape::Circle { radius: f(radius)? },
            ProfileShape::Polygon { points } => ProfileShape::Polygon {
                points: points.iter().map(&mut pt).collect::<Result<_, E>>()?,
            },
            ProfileShape::Path { start, segments } => {
                let start = pt(start)?;
                let mut out = Vec::with_capacity(segments.len());
                for seg in segments {
                    out.push(match seg {
                        PathSegment::Line { to } => PathSegment::Line { to: pt(to)? },
                        PathSegment::Arc { via, to } => PathSegment::Arc {
                            via: pt(via)?,
                            to: pt(to)?,
                        },
                    });
                }
                ProfileShape::Path {
                    start,
                    segments: out,
                }
            }
        })
    }
}

impl<S> AnalyticProfile<S> {
    /// A profile centred at the sketch origin.
    pub fn new(shape: ProfileShape<S>) -> Self
    where
        S: Default,
    {
        Self {
            shape,
            center: [S::default(), S::default()],
        }
    }

    /// Every scalar (placement first, then the shape's), in a stable
    /// order — for dependency reporting on bound scalars.
    pub fn scalars(&self) -> Vec<&S> {
        let mut out: Vec<&S> = self.center.iter().collect();
        out.extend(self.shape.scalars());
        out
    }

    /// The role name of every scalar, in [`Self::scalars`] order: the
    /// parameter names a feature tree registers for an analytic
    /// sketch (`center_x`, `center_y`, then the shape's own roles —
    /// `width`/`height`, `radius`, `p{i}_x`/`p{i}_y`, or
    /// `start_x`/`start_y` followed by `s{i}_via_x`/`s{i}_via_y`/
    /// `s{i}_x`/`s{i}_y` per path segment).
    pub fn scalar_roles(&self) -> Vec<String> {
        let mut out = vec!["center_x".to_string(), "center_y".to_string()];
        match &self.shape {
            ProfileShape::Rect { .. } => out.extend(["width".to_string(), "height".to_string()]),
            ProfileShape::Circle { .. } => out.push("radius".to_string()),
            ProfileShape::Polygon { points } => {
                for i in 0..points.len() {
                    out.push(format!("p{i}_x"));
                    out.push(format!("p{i}_y"));
                }
            }
            ProfileShape::Path { segments, .. } => {
                out.extend(["start_x".to_string(), "start_y".to_string()]);
                for (i, seg) in segments.iter().enumerate() {
                    if matches!(seg, PathSegment::Arc { .. }) {
                        out.push(format!("s{i}_via_x"));
                        out.push(format!("s{i}_via_y"));
                    }
                    out.push(format!("s{i}_x"));
                    out.push(format!("s{i}_y"));
                }
            }
        }
        out
    }

    /// Rebuild the profile with every scalar mapped through `f`.
    pub fn map_scalars<T, E>(
        &self,
        mut f: impl FnMut(&S) -> Result<T, E>,
    ) -> Result<AnalyticProfile<T>, E> {
        let center = [f(&self.center[0])?, f(&self.center[1])?];
        Ok(AnalyticProfile {
            shape: self.shape.map_scalars(f)?,
            center,
        })
    }
}

impl AnalyticProfile<f64> {
    /// Human-readable description (feature logs and errors), e.g.
    /// `circle r=0.35 at (1, 0)`.
    pub fn describe(&self) -> String {
        let at = if self.center == [0.0, 0.0] {
            String::new()
        } else {
            format!(" at ({}, {})", self.center[0], self.center[1])
        };
        match &self.shape {
            ProfileShape::Rect { width, height } => format!("rect {width}×{height}{at}"),
            ProfileShape::Circle { radius } => format!("circle r={radius}{at}"),
            ProfileShape::Polygon { points } => format!("polygon[{}]{at}", points.len()),
            ProfileShape::Path { segments, .. } => {
                format!("path[{} segs]{at}", segments.len())
            }
        }
    }
}

#[cfg(test)]
#[path = "profile_tests.rs"]
mod tests;
