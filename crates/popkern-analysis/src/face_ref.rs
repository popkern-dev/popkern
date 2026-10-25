//! References from a setup to the faces of a model.
//!
//! A face is named by a stable popkern name (from integration onward),
//! by a geometric selector (works for any source, including STEP), or
//! by raw id (fragile; warned in the report). Selectors match BRep
//! surface geometry, not triangles. A reference that matches the wrong
//! number of faces is a typed decline, never a silent pick.

use crate::model::FaceId;

/// How a setup refers to one or more faces.
#[derive(Debug, Clone, PartialEq)]
pub enum FaceRef {
    /// A stable popkern face name.
    Name(String),
    /// A geometric selector.
    Select(Selector),
    /// A raw face id. Fragile: ids change when the model changes.
    Id(FaceId),
}

/// A geometric face selector. Optional tolerances override the setup's
/// default [`Tolerance`]; angles are in degrees, distances in model
/// units.
#[derive(Debug, Clone, PartialEq)]
pub enum Selector {
    /// Planar faces with the given outward normal through a point.
    Plane {
        /// Outward normal (need not be unit length; must be non-zero).
        normal: [f64; 3],
        /// A point on the plane.
        through: [f64; 3],
        /// Angle tolerance in degrees.
        angle_tol: Option<f64>,
        /// Distance tolerance.
        dist_tol: Option<f64>,
    },
    /// Cylindrical faces with the given axis direction and radius,
    /// optionally with the axis passing through a point.
    Cylinder {
        /// Axis direction (either sense; must be non-zero).
        axis: [f64; 3],
        /// Radius.
        radius: f64,
        /// A point on the axis.
        through: Option<[f64; 3]>,
        /// Angle tolerance in degrees.
        angle_tol: Option<f64>,
        /// Distance tolerance.
        dist_tol: Option<f64>,
    },
    /// Spherical faces with the given centre and radius.
    Sphere {
        /// Centre.
        centre: [f64; 3],
        /// Radius.
        radius: f64,
        /// Distance tolerance.
        dist_tol: Option<f64>,
    },
    /// The face nearest to a point.
    NearestTo([f64; 3]),
    /// Faces matched by every selector.
    AllOf(Vec<Selector>),
    /// Faces matched by any selector.
    AnyOf(Vec<Selector>),
}

/// How many faces a reference must match.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Expect {
    /// Exactly this many.
    Exactly(u32),
    /// At least this many.
    AtLeast(u32),
}

impl Default for Expect {
    fn default() -> Self {
        Expect::Exactly(1)
    }
}

/// Default selector tolerances, stated explicitly in the setup.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize)]
pub struct Tolerance {
    /// Angle tolerance in degrees.
    pub angle_deg: f64,
    /// Distance tolerance in model units.
    pub distance: f64,
}

impl Default for Tolerance {
    /// 0.5 degrees and 0.01 mm, as in the boundary doc's example.
    fn default() -> Self {
        Tolerance {
            angle_deg: 0.5,
            distance: 0.01,
        }
    }
}

/// A non-finite or non-positive quantity found in a selector.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SelectorDefect {
    pub(crate) what: &'static str,
    pub(crate) value: f64,
}

fn finite(v: &[f64; 3]) -> bool {
    v.iter().all(|c| c.is_finite())
}

fn nonzero(v: &[f64; 3]) -> bool {
    finite(v) && v.iter().any(|&c| c != 0.0)
}

fn positive(what: &'static str, value: f64) -> Result<(), SelectorDefect> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(SelectorDefect { what, value })
    }
}

fn optional_positive(what: &'static str, value: Option<f64>) -> Result<(), SelectorDefect> {
    value.map_or(Ok(()), |v| positive(what, v))
}

fn vector(what: &'static str, v: &[f64; 3], must_be_nonzero: bool) -> Result<(), SelectorDefect> {
    let ok = if must_be_nonzero {
        nonzero(v)
    } else {
        finite(v)
    };
    if ok {
        Ok(())
    } else {
        Err(SelectorDefect { what, value: 0.0 })
    }
}

impl Selector {
    /// Checks that every vector is finite (directions non-zero), every
    /// radius and tolerance positive, and every combinator non-empty.
    pub(crate) fn check(&self) -> Result<(), SelectorDefect> {
        match self {
            Selector::Plane {
                normal,
                through,
                angle_tol,
                dist_tol,
            } => {
                vector("plane normal", normal, true)?;
                vector("plane point", through, false)?;
                optional_positive("angle tolerance", *angle_tol)?;
                optional_positive("distance tolerance", *dist_tol)
            }
            Selector::Cylinder {
                axis,
                radius,
                through,
                angle_tol,
                dist_tol,
            } => {
                vector("cylinder axis", axis, true)?;
                positive("cylinder radius", *radius)?;
                if let Some(p) = through {
                    vector("cylinder point", p, false)?;
                }
                optional_positive("angle tolerance", *angle_tol)?;
                optional_positive("distance tolerance", *dist_tol)
            }
            Selector::Sphere {
                centre,
                radius,
                dist_tol,
            } => {
                vector("sphere centre", centre, false)?;
                positive("sphere radius", *radius)?;
                optional_positive("distance tolerance", *dist_tol)
            }
            Selector::NearestTo(point) => vector("nearest-to point", point, false),
            Selector::AllOf(parts) | Selector::AnyOf(parts) => {
                if parts.is_empty() {
                    return Err(SelectorDefect {
                        what: "empty selector list",
                        value: 0.0,
                    });
                }
                parts.iter().try_for_each(Selector::check)
            }
        }
    }
}
