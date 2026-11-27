//! A cone face whose boundary sweeps round the apex close beside it,
//! meshed in the cone's development (split from `boundary_first.rs`
//! under the file-size rule).
//!
//! In the face chart a cone's azimuth is stretched by `1 / r` toward
//! the apex: a ring passing a few tolerances from the apex and turning
//! radians there (a plane missing the apex by ten tolerances cuts a
//! conic whose vertex turns about π at radius under 1e-5) runs along
//! `v ≈ 0` for radians of `u`, nearly collinear with the apex row. The
//! constrained recovery declines on those segments, the encroachment
//! splits cannot help (each split lands on the same row), and the
//! fallback filled the complement (measured: the common of a box and a
//! cone, a ruling sliver from the apex, meshed 0.08 for 1.4e-8).
//!
//! The cone is developable: unrolled about its apex, a point at slant
//! distance `s` and azimuth `u` lands at `s (cos φ, sin φ)` with
//! `φ = (u_mid − u) sin α` (against `u`, so the development keeps the
//! face chart's handedness: `(u, v)` runs as `(φ, s)` would reverse it),
//! an isometry of the face onto the plane (a face spanning at most a
//! period of `u` spans at most `2π sin α` of `φ`, so the development is
//! one-to-one; a seam's two uses are two rays from the origin). There the sweep is a
//! curve of its true size round the origin and the shared constrained
//! triangulation and deviation refinement apply unchanged.
//!
//! The development serves only the rings it was made for: a cone ring
//! whose run of samples within a hundredth of the face's reach from the
//! apex, none on it, turns more than a quarter radian of chart azimuth.
//! Every other cone face keeps the face chart. On any decline the face
//! falls back to the face-chart meshing (`None`).

use popkern_math::PortableMath;

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};

use super::BoundaryRing;

#[path = "boundary_first_cone_development_mesh.rs"]
mod mesh;
pub(super) use mesh::mesh_developed;

/// The fraction of the face's reach from the apex within which a sweep
/// is near it: the face chart's azimuth there is stretched a hundredfold
/// against the face's far side.
const NEAR_APEX: f64 = 1.0e-2;

/// The chart azimuth a run of near-apex samples turns that marks a
/// sweep (a quarter radian: a ruling's samples keep one azimuth, and a
/// rim lies far from the apex).
const SWEEP_STEP: f64 = 0.25;

/// A cone's development about its apex: `(x, y) = s (cos φ, sin φ)`,
/// `s` the slant distance, `φ = (u_mid − u) sin α`.
pub(super) struct Development<'a> {
    surface: &'a SurfaceGeometry,
    apex: Point3,
    /// The azimuth frame (`u = 0` along `e1`) and the axis.
    frame: (Vec3, Vec3, Vec3),
    /// `sin α`, `cos α` of the half-angle.
    sin_a: f64,
    cos_a: f64,
    /// The face's mid azimuth (φ = 0 there).
    u_mid: f64,
    /// The nappe: `+1` where `v ≥ 0`, `−1` on the other.
    nappe: f64,
    domain: ((f64, f64), (f64, f64)),
}

impl Development<'_> {
    /// The development coordinates of a face point at `(u, v)`.
    fn develop(&self, (u, v): (f64, f64)) -> (f64, f64) {
        let s = (v * self.nappe).max(0.0) / self.cos_a;
        let phi = (self.u_mid - u) * self.sin_a;
        (s * phi.pcos(), s * phi.psin())
    }

    /// The development coordinates of a canonical position, its azimuth
    /// read from the position (the representative nearest `u_near`, the
    /// sample's chart u): a position on the apex lands on the origin.
    fn develop_point(&self, p: &Point3, u_near: f64) -> (f64, f64) {
        let d = *p - self.apex;
        let (e1, e2, axis) = self.frame;
        let s = (d.dot(&axis) * self.nappe).max(0.0) / self.cos_a;
        let radial = d - axis * d.dot(&axis);
        if radial.length() <= 0.0 {
            return (0.0, 0.0);
        }
        let az = radial.dot(&e2).patan2(radial.dot(&e1));
        let tau = std::f64::consts::TAU;
        let u = az + tau * ((u_near - az) / tau).round();
        self.develop((u, s * self.cos_a * self.nappe))
    }

    /// The face parameters at development coordinates `(x, y)`.
    pub(super) fn face_uv(&self, x: f64, y: f64) -> (f64, f64) {
        let s = x.phypot(y);
        let phi = if s > 0.0 { y.patan2(x) } else { 0.0 };
        (self.u_mid - phi / self.sin_a, self.nappe * s * self.cos_a)
    }
}

impl SurfaceEval for Development<'_> {
    fn point_at(&self, x: f64, y: f64) -> Point3 {
        let (u, v) = self.face_uv(x, y);
        self.surface.point_at(u, v)
    }

    fn partial_u(&self, x: f64, y: f64) -> Vec3 {
        self.first_partials(x, y).0
    }

    fn partial_v(&self, x: f64, y: f64) -> Vec3 {
        self.first_partials(x, y).1
    }

    /// The chain rule through the inverse: `∂(u, v)/∂(x, y)` from
    /// `u = u_mid − atan2(y, x) / sin α`, `v = ± s cos α`. At the apex the
    /// development is the tangent cone's isometry and has no partials;
    /// the frame there is the plane's own.
    fn first_partials(&self, x: f64, y: f64) -> (Vec3, Vec3) {
        let s2 = x * x + y * y;
        if s2 <= 0.0 {
            return (Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
        }
        let s = s2.sqrt();
        let (u, v) = self.face_uv(x, y);
        let (su, sv) = self.surface.first_partials(u, v);
        let (du_dx, du_dy) = (y / (s2 * self.sin_a), -x / (s2 * self.sin_a));
        let (dv_dx, dv_dy) = (
            self.nappe * self.cos_a * x / s,
            self.nappe * self.cos_a * y / s,
        );
        (su * du_dx + sv * dv_dx, su * du_dy + sv * dv_dy)
    }

    fn normal_at(&self, x: f64, y: f64) -> Vec3 {
        let (u, v) = self.face_uv(x, y);
        self.surface.normal_at(u, v)
    }

    fn parameter_range(&self) -> ((f64, f64), (f64, f64)) {
        self.domain
    }

    fn closest_point(&self, point: &Point3) -> Option<(f64, f64)> {
        let uv = self.surface.closest_point(point)?;
        // The representative within half a period of the mid azimuth.
        let p = std::f64::consts::TAU;
        let u = uv.0 - p * ((uv.0 - self.u_mid) / p).round();
        Some(self.develop((u, uv.1)))
    }

    fn bounding_box(&self, u_range: (f64, f64), v_range: (f64, f64)) -> Aabb {
        popkern_core::geometry::surface_helpers::bounding_box_by_sampling(self, u_range, v_range)
    }
}

/// The cone's development for a face's built rings and each ring
/// sample's azimuth (read from its position, unwrapped along its ring,
/// the holes placed on the outer ring's sheet), or `None` when the face
/// is not a cone face sweeping beside its apex, winds round the axis, or
/// reaches both nappes.
///
/// The face chart's own `u` is not used for the development: a ruling's
/// pcurve ending on the apex carries an arbitrary `u` there, and its
/// samples beside the apex sweep the chart along the apex row where the
/// positions stand on one ruling (the face-chart unwrap then reads a full
/// wrap). It only marks the sweep.
pub(super) fn development_for<'a>(
    surface: &'a SurfaceGeometry,
    rings: &[BoundaryRing],
) -> Option<(Development<'a>, Vec<Vec<f64>>)> {
    let SurfaceGeometry::Cone {
        apex,
        axis,
        half_angle,
    } = surface
    else {
        return None;
    };
    let (sin_a, cos_a) = (half_angle.psin(), half_angle.pcos());
    if !(sin_a > 0.0 && cos_a > 0.0) || rings.is_empty() {
        return None;
    }
    let tau = std::f64::consts::TAU;
    let positions = || rings.iter().flat_map(|r| r.pos.iter());
    let reach = positions().map(|p| p.distance_to(apex)).fold(0.0, f64::max);
    if reach.is_nan() || reach <= 0.0 {
        return None;
    }
    let near = NEAR_APEX * reach;
    // A sample on the apex has no azimuth: the apex row in space, a
    // millionth of the reach.
    let on_apex = 1.0e-6 * reach;
    // The azimuth a run of near-apex segments turns in the face chart,
    // summed (the chart may spread one position's samples along the row
    // in small steps), each step its nearest representative.
    let sweeps = rings.iter().any(|r| {
        let n = r.uv.len();
        let mut run = 0.0;
        (0..2 * n).any(|i| {
            let (k, j) = (i % n, (i + 1) % n);
            let (da, db) = (r.pos[k].distance_to(apex), r.pos[j].distance_to(apex));
            if da > on_apex && db > on_apex && da < near && db < near {
                let du = r.uv[j].0 - r.uv[k].0;
                run += (du - tau * (du / tau).round()).abs();
            } else {
                run = 0.0;
            }
            run > SWEEP_STEP
        })
    });
    if !sweeps {
        return None;
    }
    let declined = |reason: &str| {
        if popkern_core::trace::trace_enabled("boundary-first") {
            eprintln!("    [bf] cone face sweeping beside its apex: no development ({reason})");
        }
    };
    // One nappe: the development is then one-to-one.
    let heights: Vec<f64> = positions().map(|p| (*p - *apex).dot(axis)).collect();
    let nappe = if heights.iter().all(|&h| h >= -on_apex) {
        1.0
    } else if heights.iter().all(|&h| h <= on_apex) {
        -1.0
    } else {
        declined("both nappes");
        return None;
    };
    let (e1, e2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(*axis);
    let mut azimuths: Vec<Vec<f64>> = Vec::with_capacity(rings.len());
    for r in rings {
        let Some(mut az) = ring_azimuths(r, (*apex, e1, e2, *axis), on_apex) else {
            declined("a ring winds round the axis");
            return None;
        };
        // A hole on the outer ring's sheet.
        if let Some(outer) = azimuths.first() {
            let mean = |a: &[f64]| a.iter().sum::<f64>() / a.len().max(1) as f64;
            let shift = tau * ((mean(outer) - mean(&az)) / tau).round();
            az.iter_mut().for_each(|a| *a += shift);
        }
        azimuths.push(az);
    }
    let all = || azimuths.iter().flatten().copied();
    let (u_lo, u_hi) = (
        all().fold(f64::INFINITY, f64::min),
        all().fold(f64::NEG_INFINITY, f64::max),
    );
    // One sheet of the development: the azimuth span, developed, under
    // a full turn (a face spanning a whole period between its seam's two
    // uses develops to `2π sin α`).
    let developed = (u_hi - u_lo) * sin_a;
    if developed.is_nan() || developed >= tau {
        declined("the rings span a period");
        return None;
    }
    let (h_lo, h_hi) = (
        heights.iter().copied().fold(f64::INFINITY, f64::min),
        heights.iter().copied().fold(f64::NEG_INFINITY, f64::max),
    );
    Some((
        Development {
            surface,
            apex: *apex,
            frame: (e1, e2, *axis),
            sin_a,
            cos_a,
            u_mid: 0.5 * (u_lo + u_hi),
            nappe,
            domain: ((u_lo, u_hi), (h_lo, h_hi)),
        },
        azimuths,
    ))
}

/// A ring's sample azimuths about the axis from their positions,
/// unwrapped along the ring (a sample on the apex takes its
/// predecessor's), or `None` when every sample lies on the apex or the
/// ring winds round the axis other than once through the apex.
///
/// A step between two samples turning less than a quarter period is read
/// from the positions alone (where a ruling's pcurve carries an arbitrary
/// `u` beside the apex, its positions keep one azimuth). A larger step is
/// ambiguous from positions (a conic turning round the apex within a few
/// samples steps past half a period), and takes the whole-period
/// representative nearest the pcurve's own step, the curve's continuity
/// (along one edge; across edges, its nearest representative).
///
/// The walk starts at the first sample off the apex after one on it, so
/// a ring closing round the axis through the apex (a face spanning a
/// whole period between its seam's two uses, which meet there) takes its
/// period where it passes the apex: there the development is continuous,
/// the seam's two uses two rays from the origin.
fn ring_azimuths(
    ring: &BoundaryRing,
    (apex, e1, e2, axis): (Point3, Vec3, Vec3, Vec3),
    on_apex: f64,
) -> Option<Vec<f64>> {
    let tau = std::f64::consts::TAU;
    let n = ring.pos.len();
    let raw = |p: &Point3| -> Option<f64> {
        let d = *p - apex;
        let radial = d - axis * d.dot(&axis);
        (d.length() > on_apex && radial.length() > 0.0)
            .then(|| radial.dot(&e2).patan2(radial.dot(&e1)))
    };
    let nearest = |x: f64, to: f64| x + tau * ((to - x) / tau).round();
    let regular: Vec<bool> = ring.pos.iter().map(|p| raw(p).is_some()).collect();
    let after_apex = (0..n).find(|&i| regular[i] && !regular[(i + n - 1) % n]);
    let first = after_apex.or_else(|| (0..n).find(|&i| regular[i]))?;
    let mut az = vec![0.0; n];
    let start = raw(&ring.pos[first])?;
    let (mut prev, mut prev_i) = (start, first);
    for step in 1..=n {
        let i = (first + step) % n;
        if let Some(a) = raw(&ring.pos[i]) {
            let d = nearest(a - prev, 0.0);
            let d = if d.abs() < 0.25 * tau {
                d
            } else {
                // Along one edge the pcurve's u is continuous; across
                // edges their branches may differ by whole periods.
                let du = ring.uv[i].0 - ring.uv[prev_i].0;
                let one_edge = ring.sources[prev_i].0 == ring.sources[(i + n - 1) % n].0;
                nearest(d, if one_edge { du } else { nearest(du, 0.0) })
            };
            prev += d;
            prev_i = i;
        }
        if step < n {
            az[i] = prev;
        }
    }
    az[first] = start;
    // Back round to the first sample: the whole periods the ring turned.
    let turns = ((prev - start) / tau).round();
    (turns == 0.0 || (turns.abs() == 1.0 && after_apex.is_some())).then_some(az)
}

#[cfg(test)]
#[path = "boundary_first_cone_development_tests.rs"]
mod tests;
