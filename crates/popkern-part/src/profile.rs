//! 2D profiles: exact analytic profiles on a plane, compiled to the
//! kernel's NURBS wires.
//!
//! Exactness matters: the sweep kernel's canonical recognition turns
//! an extruded rational-quadratic circle into a true `Cylinder`
//! surface (and extruded lines into `Plane`s) — approximated
//! profiles would silently degrade every downstream operation
//! (booleans, fillets, STEP) to NURBS. The compilation itself lives in
//! `popkern-sketch` ([`AnalyticProfile`]) so the feature tree's
//! analytic sketch arm and this builder share one exact compiler
//! (PK-364); this module is the fluent builder over it plus the
//! extrude entry point.

use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_sketch::{AnalyticProfile, PathSegment, ProfileShape, SketchPlane};
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

    /// The plane as a sketch frame: global X/Y axes at height `z`.
    fn sketch_plane(&self) -> SketchPlane {
        SketchPlane {
            origin: Point3::new(0.0, 0.0, self.z),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        }
    }
}

/// Builder for a mixed line/arc sketch outline. Created by
/// [`Profile::path`]; finished by [`PathBuilder::close`]. An outline
/// whose last segment does not end at the start point is closed by a
/// straight segment back to it at compile time.
#[derive(Debug, Clone)]
pub struct PathBuilder {
    profile: Profile,
    start: (f64, f64),
    segments: Vec<PathSegment<f64>>,
}

impl PathBuilder {
    /// Straight segment to `(x, y)`.
    pub fn line_to(mut self, x: f64, y: f64) -> PathBuilder {
        self.segments.push(PathSegment::Line { to: [x, y] });
        self
    }

    /// Circular arc THROUGH `via` ending at `to`.
    pub fn arc_to(mut self, via: (f64, f64), to: (f64, f64)) -> PathBuilder {
        self.segments.push(PathSegment::Arc {
            via: [via.0, via.1],
            to: [to.0, to.1],
        });
        self
    }

    /// Close the outline and return the profile.
    pub fn close(self) -> Profile {
        let mut profile = self.profile;
        profile.shape = Some(ProfileShape::Path {
            start: [self.start.0, self.start.1],
            segments: self.segments,
        });
        profile
    }
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
    shape: Option<ProfileShape<f64>>,
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
        self.shape = Some(ProfileShape::Rect { width, height });
        self
    }

    /// A circle centred on the sketch position.
    pub fn circle(mut self, radius: f64) -> Profile {
        self.shape = Some(ProfileShape::Circle { radius });
        self
    }

    /// A closed simple polygon in sketch coordinates (do not repeat
    /// the first point). Winding may be either direction.
    pub fn polygon(mut self, points: &[(f64, f64)]) -> Profile {
        self.shape = Some(ProfileShape::Polygon {
            points: points.iter().map(|&(x, y)| [x, y]).collect(),
        });
        self
    }

    /// Start a mixed line/arc outline at `(x, y)` (see
    /// [`PathBuilder`]). Arc segments compile to exact rational
    /// quadratics; the sweep kernel recognizes their extruded walls
    /// as true `Cylinder` surfaces.
    pub fn path(self, x: f64, y: f64) -> PathBuilder {
        PathBuilder {
            profile: self,
            start: (x, y),
            segments: Vec::new(),
        }
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
        // The sweep recentres the wire onto ITS OWN sampled centroid
        // before placing it on the path. Rect and circle wires are
        // symmetric about the origin so this is a no-op, but a
        // POLYGON wire's centroid is wherever its outline puts it —
        // every polygon extrude was silently translated by
        // −centroid (measured: a trapezoid deck landed 23 mm from
        // its authored coordinates and every downstream tool
        // missed). Fold the same sampled centroid into the path so
        // the recentring cancels and authored coordinates are
        // preserved.
        let c = {
            use popkern_core::geometry::curve::CurveEval;
            let (t0, t1) = wire.parameter_range();
            let n = 64;
            let mut sx = 0.0;
            let mut sy = 0.0;
            for i in 0..n {
                let t = t0 + (t1 - t0) * (i as f64) / (n as f64);
                let p = wire.point_at(t);
                sx += p.x;
                sy += p.y;
            }
            (sx / n as f64, sy / n as f64)
        };
        let (cx, cy) = (self.center.0 + c.0, self.center.1 + c.1);
        if popkern_core::trace::trace_enabled("part") {
            eprintln!(
                "  [extrude] c=({:.4},{:.4}) center=({:.4},{:.4}) path=({cx:.4},{cy:.4})",
                c.0, c.1, self.center.0, self.center.1
            );
        }
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

    /// The analytic description this builder holds, if a shape was
    /// given — the same record the feature tree's analytic sketch arm
    /// stores.
    fn analytic(&self) -> Result<AnalyticProfile<f64>, PartError> {
        let shape = self.shape.clone().ok_or_else(|| {
            PartError::NotApplicable(
                "sketch has no profile: call .rect(), .circle(), or .polygon() first".into(),
            )
        })?;
        Ok(AnalyticProfile {
            shape,
            center: [self.center.0, self.center.1],
        })
    }

    /// Human-readable profile description (feature log / errors).
    fn describe(&self) -> String {
        match self.analytic() {
            Ok(p) => p.describe(),
            Err(_) => "empty sketch".into(),
        }
    }

    /// Compile the profile to an exact closed NURBS wire in the
    /// profile's LOCAL frame (centred at the origin; placement is
    /// carried by the sweep path).
    fn wire_local(&self) -> Result<NurbsCurve, PartError> {
        self.analytic()?
            .local_wire(&self.plane.sketch_plane())
            .map_err(|e| PartError::NotApplicable(e.to_string()))
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
