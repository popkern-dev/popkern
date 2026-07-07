//! 2D profiles: exact analytic profiles on a plane, compiled to the
//! kernel's NURBS wires.
//!
//! Exactness matters: the sweep kernel's canonical recognition turns
//! an extruded rational-quadratic circle into a true `Cylinder`
//! surface (and extruded lines into `Plane`s) — approximated
//! profiles would silently degrade every downstream operation
//! (booleans, fillets, STEP) to NURBS. Profiles here always compile
//! to the exact wire forms the kernel recognizes.

use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_sweep::{SweepOptions, SweepProfile, sweep};

use crate::error::PartError;
use crate::part::Part;

/// A sketch plane. v1 supports planes parallel to global XY (which
/// the flagship-class parts need); general placements arrive with
/// sketch-on-face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    /// Height of the plane along +Z.
    z: f64,
}

impl Plane {
    /// The global XY plane (z = 0).
    #[allow(non_upper_case_globals)]
    pub const XY: Plane = Plane { z: 0.0 };

    /// A plane parallel to XY at the given height.
    pub fn xy_at(z: f64) -> Plane {
        Plane { z }
    }

    /// The plane's height along +Z.
    pub fn z(&self) -> f64 {
        self.z
    }
}

/// The profile shapes a profile can hold (one closed shape per
/// profile in v1).
#[derive(Debug, Clone)]
enum Shape {
    /// Axis-aligned rectangle: width along X, height along Y.
    Rect { width: f64, height: f64 },
    /// Circle of the given radius.
    Circle { radius: f64 },
    /// Closed simple polygon (planar, first point NOT repeated).
    Polygon { points: Vec<(f64, f64)> },
}

/// A 2D sketch: one closed profile on a plane, positioned in-plane.
///
/// ```no_run
/// # use popkern_part::prelude::*;
/// let solid = Profile::on(Plane::XY).circle(0.35).at(1.0, 0.0).extrude(0.8)?;
/// # Ok::<(), PartError>(())
/// ```
#[derive(Debug, Clone)]
pub struct Profile {
    plane: Plane,
    shape: Option<Shape>,
    center: (f64, f64),
    tolerance: Tolerance,
}

impl Profile {
    /// Start a sketch on a plane with the default tolerance.
    pub fn on(plane: Plane) -> Profile {
        Profile {
            plane,
            shape: None,
            center: (0.0, 0.0),
            tolerance: Tolerance::default(),
        }
    }

    /// Override the tolerance carried into the part.
    pub fn with_tolerance(mut self, tolerance: Tolerance) -> Profile {
        self.tolerance = tolerance;
        self
    }

    /// A rectangle centred on the sketch position (see [`Profile::at`]).
    pub fn rect(mut self, width: f64, height: f64) -> Profile {
        self.shape = Some(Shape::Rect { width, height });
        self
    }

    /// A circle centred on the sketch position.
    pub fn circle(mut self, radius: f64) -> Profile {
        self.shape = Some(Shape::Circle { radius });
        self
    }

    /// A closed simple polygon in sketch coordinates (do not repeat
    /// the first point). Winding may be either direction.
    pub fn polygon(mut self, points: &[(f64, f64)]) -> Profile {
        self.shape = Some(Shape::Polygon {
            points: points.to_vec(),
        });
        self
    }

    /// Position the profile centre in the sketch plane.
    pub fn at(mut self, x: f64, y: f64) -> Profile {
        self.center = (x, y);
        self
    }

    /// Extrude the profile by `depth` along +Z (negative depths
    /// extrude downward), producing a validated solid [`Part`].
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] when the sketch has no profile;
    /// [`PartError::Sweep`] when the sweep kernel declines;
    /// [`PartError::Validation`] if the swept body fails validation.
    pub fn extrude(self, depth: f64) -> Result<Part, PartError> {
        let feature = self.describe();
        // The sweep kernel positions the profile relative to the
        // PATH frame (it recentres the wire onto the path start), so
        // in-plane placement must be expressed through the path —
        // a wire built at (cx, cy) would be silently recentred to
        // the origin (measured: a circle profile at x=1 produced a
        // solid at x=0).
        let wire = self.wire_local()?;
        let (cx, cy) = self.center;
        let (z0, z1) = if depth >= 0.0 {
            (self.plane.z, self.plane.z + depth)
        } else {
            (self.plane.z + depth, self.plane.z)
        };
        let path = line_path(cx, cy, z0, z1);
        let result = sweep(
            &SweepProfile::Wire(wire),
            &path,
            &SweepOptions::default(),
            self.tolerance,
            None,
        )
        .map_err(|source| PartError::Sweep {
            feature: feature.clone(),
            source,
        })?;
        Part::from_validated(
            result.brep,
            self.tolerance,
            "extrude",
            format!("{feature} depth {depth}"),
        )
    }

    /// Human-readable profile description (feature log / errors).
    fn describe(&self) -> String {
        let at = if self.center == (0.0, 0.0) {
            String::new()
        } else {
            format!(" at ({}, {})", self.center.0, self.center.1)
        };
        match &self.shape {
            Some(Shape::Rect { width, height }) => format!("rect {width}×{height}{at}"),
            Some(Shape::Circle { radius }) => format!("circle r={radius}{at}"),
            Some(Shape::Polygon { points }) => format!("polygon[{}]{at}", points.len()),
            None => "empty sketch".into(),
        }
    }

    /// Compile the profile to an exact closed NURBS wire in the
    /// profile's LOCAL frame (centred at the origin; placement is
    /// carried by the sweep path).
    fn wire_local(&self) -> Result<NurbsCurve, PartError> {
        let (cx, cy) = (0.0, 0.0);
        let z = self.plane.z;
        match &self.shape {
            Some(Shape::Rect { width, height }) => {
                let (hw, hh) = (width / 2.0, height / 2.0);
                Ok(polygon_wire(
                    &[
                        (cx - hw, cy - hh),
                        (cx + hw, cy - hh),
                        (cx + hw, cy + hh),
                        (cx - hw, cy + hh),
                    ],
                    z,
                ))
            }
            Some(Shape::Circle { radius }) => Ok(circle_wire(cx, cy, *radius, z)),
            Some(Shape::Polygon { points }) => {
                if points.len() < 3 {
                    return Err(PartError::NotApplicable(format!(
                        "polygon profile needs at least 3 points, got {}",
                        points.len()
                    )));
                }
                let placed: Vec<(f64, f64)> =
                    points.iter().map(|&(x, y)| (x + cx, y + cy)).collect();
                Ok(polygon_wire(&placed, z))
            }
            None => Err(PartError::NotApplicable(
                "sketch has no profile: call .rect(), .circle(), or .polygon() first".into(),
            )),
        }
    }
}

/// Closed degree-1 wire through the given 2D points at height `z`.
fn polygon_wire(points: &[(f64, f64)], z: f64) -> NurbsCurve {
    let mut pts: Vec<Point3> = points.iter().map(|&(x, y)| Point3::new(x, y, z)).collect();
    pts.push(pts[0]);
    let n = pts.len();
    let mut knots = vec![0.0];
    for i in 0..n {
        knots.push(i as f64 / (n - 1) as f64);
    }
    knots.push(1.0);
    NurbsCurve {
        degree: 1,
        knots,
        control_points: pts,
        weights: vec![1.0; n],
    }
}

/// Exact rational-quadratic circle wire (standard 9-CP form) — the
/// form the sweep kernel's canonical recognition maps to a true
/// `Cylinder`.
fn circle_wire(cx: f64, cy: f64, r: f64, z: f64) -> NurbsCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let pts = vec![
        Point3::new(cx + r, cy, z),
        Point3::new(cx + r, cy + r, z),
        Point3::new(cx, cy + r, z),
        Point3::new(cx - r, cy + r, z),
        Point3::new(cx - r, cy, z),
        Point3::new(cx - r, cy - r, z),
        Point3::new(cx, cy - r, z),
        Point3::new(cx + r, cy - r, z),
        Point3::new(cx + r, cy, z),
    ];
    let weights = vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
    let knots = vec![
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ];
    NurbsCurve {
        degree: 2,
        knots,
        control_points: pts,
        weights,
    }
}

/// Straight extrusion path from `z0` to `z1` along +Z at the
/// in-plane placement `(cx, cy)`.
fn line_path(cx: f64, cy: f64, z0: f64, z1: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(cx, cy, z0), Point3::new(cx, cy, z1)],
        weights: vec![1.0, 1.0],
    }
}
