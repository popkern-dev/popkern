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
    /// Mixed line/arc path (closed; authored via [`Profile::path`]).
    Path {
        start: (f64, f64),
        segs: Vec<PathSeg>,
    },
}

/// One segment of a mixed sketch path.
#[derive(Debug, Clone)]
enum PathSeg {
    /// Straight segment to `to`.
    Line { to: (f64, f64) },
    /// Circular arc from the current point THROUGH `via` to `to`
    /// (three-point form: unambiguous, no sweep-direction flags).
    Arc { via: (f64, f64), to: (f64, f64) },
}

/// Builder for a mixed line/arc sketch outline. Created by
/// [`Profile::path`]; finished by [`PathBuilder::close`], which
/// closes the outline back to the start point with a straight
/// segment unless the last segment already ends there.
#[derive(Debug, Clone)]
pub struct PathBuilder {
    profile: Profile,
    start: (f64, f64),
    segs: Vec<PathSeg>,
}

impl PathBuilder {
    /// Straight segment to `(x, y)`.
    pub fn line_to(mut self, x: f64, y: f64) -> PathBuilder {
        self.segs.push(PathSeg::Line { to: (x, y) });
        self
    }

    /// Circular arc THROUGH `via` ending at `to`.
    pub fn arc_to(mut self, via: (f64, f64), to: (f64, f64)) -> PathBuilder {
        self.segs.push(PathSeg::Arc { via, to });
        self
    }

    /// Close the outline and return the profile.
    pub fn close(mut self) -> Profile {
        let mut profile = self.profile;
        // Ensure closure back to the start.
        let end = match self.segs.last() {
            Some(PathSeg::Line { to }) | Some(PathSeg::Arc { to, .. }) => *to,
            None => self.start,
        };
        let dx = end.0 - self.start.0;
        let dy = end.1 - self.start.1;
        if (dx * dx + dy * dy).sqrt() > 1.0e-12 {
            self.segs.push(PathSeg::Line { to: self.start });
        }
        profile.shape = Some(Shape::Path {
            start: self.start,
            segs: self.segs,
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

    /// Start a mixed line/arc outline at `(x, y)` (see
    /// [`PathBuilder`]). Arc segments compile to exact rational
    /// quadratics; the sweep kernel recognizes their extruded walls
    /// as true `Cylinder` surfaces.
    pub fn path(self, x: f64, y: f64) -> PathBuilder {
        PathBuilder {
            profile: self,
            start: (x, y),
            segs: Vec::new(),
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
            Some(Shape::Path { segs, .. }) => format!("path[{} segs]{at}", segs.len()),
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
            Some(Shape::Path { start, segs }) => {
                if segs.len() < 2 {
                    return Err(PartError::NotApplicable(
                        "path profile needs at least 2 segments".into(),
                    ));
                }
                mixed_path_wire(*start, segs, z)
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

/// Compile a mixed line/arc path into a single closed degree-2
/// rational NURBS wire. Every span joint carries FULL-multiplicity
/// knots; arc spans wider than 90° are subdivided so each rational
/// quadratic stays well-conditioned. Line spans are exact (collinear
/// control points, unit weights); arc spans are exact rational
/// quadratics (mid-weight cos(Δ/2)).
fn mixed_path_wire(start: (f64, f64), segs: &[PathSeg], z: f64) -> Result<NurbsCurve, PartError> {
    // Each entry: (P0, Pm, P2, w_mid) — one quadratic span.
    let mut spans: Vec<(Point3, Point3, Point3, f64)> = Vec::new();
    let mut cur = start;
    let p3 = |p: (f64, f64)| Point3::new(p.0, p.1, z);
    for seg in segs {
        match seg {
            PathSeg::Line { to } => {
                let a = p3(cur);
                let b = p3(*to);
                let m = Point3::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0, z);
                spans.push((a, m, b, 1.0));
                cur = *to;
            }
            PathSeg::Arc { via, to } => {
                let (ax, ay) = cur;
                let (vx, vy) = *via;
                let (bx, by) = *to;
                // Circumcenter of the three points.
                let d = 2.0 * (ax * (vy - by) + vx * (by - ay) + bx * (ay - vy));
                if d.abs() < 1.0e-12 {
                    return Err(PartError::NotApplicable(
                        "arc_to points are collinear".into(),
                    ));
                }
                let aa = ax * ax + ay * ay;
                let vv = vx * vx + vy * vy;
                let bb = bx * bx + by * by;
                let ox = (aa * (vy - by) + vv * (by - ay) + bb * (ay - vy)) / d;
                let oy = (aa * (bx - vx) + vv * (ax - bx) + bb * (vx - ax)) / d;
                let r = ((ax - ox).powi(2) + (ay - oy).powi(2)).sqrt();
                let a0 = (ay - oy).atan2(ax - ox);
                let a_v = (vy - oy).atan2(vx - ox);
                let a1 = (by - oy).atan2(bx - ox);
                // Choose the sweep direction that passes through via.
                let norm = |t: f64| t.rem_euclid(std::f64::consts::TAU);
                let ccw_v = norm(a_v - a0);
                let ccw_1 = norm(a1 - a0);
                let (sweep, dir) = if ccw_v <= ccw_1 {
                    (ccw_1, 1.0)
                } else {
                    (std::f64::consts::TAU - ccw_1, -1.0)
                };
                let n_sub = (sweep / std::f64::consts::FRAC_PI_2).ceil().max(1.0) as usize;
                let step = dir * sweep / n_sub as f64;
                for i in 0..n_sub {
                    let sa = a0 + step * i as f64;
                    let sb = sa + step;
                    let half = (step / 2.0).abs();
                    let w = half.cos();
                    let pa = Point3::new(ox + r * sa.cos(), oy + r * sa.sin(), z);
                    let pb = Point3::new(ox + r * sb.cos(), oy + r * sb.sin(), z);
                    let mid_ang = (sa + sb) / 2.0;
                    let pm = Point3::new(
                        ox + (r / w) * mid_ang.cos(),
                        oy + (r / w) * mid_ang.sin(),
                        z,
                    );
                    spans.push((pa, pm, pb, w));
                }
                cur = *to;
            }
        }
    }
    // Assemble: clamped degree-2, double interior knots at every
    // span joint.
    let n_spans = spans.len();
    let mut control_points: Vec<Point3> = Vec::with_capacity(2 * n_spans + 1);
    let mut weights: Vec<f64> = Vec::with_capacity(2 * n_spans + 1);
    for (i, (a, m, b, w)) in spans.iter().enumerate() {
        if i == 0 {
            control_points.push(*a);
            weights.push(1.0);
        }
        control_points.push(*m);
        weights.push(*w);
        control_points.push(*b);
        weights.push(1.0);
    }
    let mut knots: Vec<f64> = vec![0.0, 0.0, 0.0];
    for i in 1..n_spans {
        let k = i as f64 / n_spans as f64;
        knots.push(k);
        knots.push(k);
    }
    knots.extend_from_slice(&[1.0, 1.0, 1.0]);
    Ok(NurbsCurve {
        degree: 2,
        knots,
        control_points,
        weights,
    })
}
