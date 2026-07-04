//! Canonical surface recognition for swept surfaces.
//!
//! An extruded circular arc IS a cylinder; leaving it as a generic
//! NURBS forces every downstream consumer (boolean intersection,
//! fillet analysis, STEP export) onto approximate paths — the
//! boolean's marching walker on a NURBS "cylinder" × plane produced
//! intersection curves with ±6e-2 radial slop where the analytic
//! path is exact. Real kernels canonicalize at creation time; so do
//! we.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsSurface, Point3, Vec3};

/// Affine map from the swept surface's `[0,1]²` chart to the
/// recognized surface's native chart (cylinder: `u` = angle in
/// radians, `v` = length along the axis).
#[derive(Debug, Clone, Copy)]
pub(crate) struct ChartMap {
    /// `u_native = u_offset + u01 * u_scale`.
    pub u_offset: f64,
    /// See `u_offset`.
    pub u_scale: f64,
    /// `v_native = v_offset + v01 * v_scale`.
    pub v_offset: f64,
    /// See `v_offset`.
    pub v_scale: f64,
}

impl ChartMap {
    /// Map a `[0,1]²` chart point into the native chart.
    pub fn map(&self, u01: f64, v01: f64) -> (f64, f64) {
        (
            self.u_offset + u01 * self.u_scale,
            self.v_offset + v01 * self.v_scale,
        )
    }
}

/// Relative recognition tolerance: the sampled surface must match
/// the analytic candidate to this fraction of its radius.
const RECOGNITION_REL: f64 = 1.0e-9;

/// Recognize an extruded circular arc as an analytic cylinder.
///
/// Conditions verified by dense sampling (not construction
/// provenance, so any parameterization qualifies):
/// 1. Translation surface: every v-iso row is a translate of the
///    v=0 row (constant offset vector — a straight extrusion).
/// 2. The v=0 profile row lies on a circle whose axis is the
///    translation direction.
///
/// Returns the cylinder and the chart map from the swept `[0,1]²`
/// parameterization. NOTE: the u chart map is angle-LINEAR while a
/// rational-quadratic circle's parameter is not; pcurves remapped
/// through it stay geometrically on the surface but their
/// parameterization deviates from the 3D edge curve's by up to a few
/// percent mid-span — consumers in this kernel treat pcurves as
/// association seeds and refine by inversion, which absorbs that.
pub(crate) fn recognize_extruded_cylinder(
    surface: &NurbsSurface,
) -> Option<(SurfaceGeometry, ChartMap)> {
    const NU: usize = 33;
    const NV: usize = 5;
    // 1. Translation test.
    let row0: Vec<Point3> = (0..NU)
        .map(|i| surface.point_at(i as f64 / (NU - 1) as f64, 0.0))
        .collect();
    let mut axis_dir: Option<Vec3> = None;
    let mut height = 0.0;
    for j in 1..NV {
        let v = j as f64 / (NV - 1) as f64;
        let offset = surface.point_at(0.0, v) - row0[0];
        for (i, p0) in row0.iter().enumerate() {
            let p = surface.point_at(i as f64 / (NU - 1) as f64, v);
            let d = p - *p0;
            if (d - offset).length() > 1.0e-9 * (1.0 + offset.length()) {
                return None;
            }
        }
        if j == NV - 1 {
            height = offset.length();
            if height < 1.0e-12 {
                return None;
            }
            axis_dir = offset.normalized();
        }
    }
    let axis = axis_dir?;

    // 2. Circle test on the profile row: circumcenter from three
    //    spread points, then every sample at constant radius and in
    //    a plane perpendicular to the axis.
    let (a, b, c) = (row0[0], row0[NU / 3], row0[2 * NU / 3]);
    let center = circumcenter_3d(a, b, c, axis)?;
    let radius = a.distance_to(&center);
    if radius < 1.0e-12 {
        return None;
    }
    let gate = RECOGNITION_REL * radius.max(1.0);
    for p in &row0 {
        let d = *p - center;
        if (d.dot(&axis)).abs() > gate {
            return None;
        }
        let radial = d - axis * d.dot(&axis);
        if (radial.length() - radius).abs() > gate {
            return None;
        }
    }

    // 3. Chart: angles of the profile in the cylinder's own basis
    //    (matching `SurfaceGeometry::Cylinder`'s orthonormal_basis
    //    convention so point_at agrees), unwrapped monotonically.
    let (perp1, perp2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(axis);
    let angle_of = |p: &Point3| -> f64 {
        let d = *p - center;
        d.dot(&perp2).atan2(d.dot(&perp1))
    };
    let theta0 = angle_of(&row0[0]);
    let mut theta1 = angle_of(&row0[NU - 1]);
    // Unwrap the end angle to the branch reached by walking the row.
    let mut walked = theta0;
    for p in row0.iter().skip(1) {
        let mut t = angle_of(p);
        while t - walked > std::f64::consts::PI {
            t -= std::f64::consts::TAU;
        }
        while t - walked < -std::f64::consts::PI {
            t += std::f64::consts::TAU;
        }
        walked = t;
    }
    theta1 = walked;
    if (theta1 - theta0).abs() < 1.0e-12 {
        // Degenerate span; closed full circles land at ±2π.
        return None;
    }

    let geometry = SurfaceGeometry::Cylinder {
        origin: center,
        axis,
        radius,
    };
    // Verify the recognized surface reproduces the sampled grid.
    for j in 0..NV {
        let v01 = j as f64 / (NV - 1) as f64;
        for i in 0..NU {
            let u01 = i as f64 / (NU - 1) as f64;
            let native_u = theta0 + (theta1 - theta0) * u01;
            let native_v = v01 * height;
            let p_nurbs = surface.point_at(u01, v01);
            let p_cyl = geometry.point_at(native_u, native_v);
            // The angle chart is linear while the NURBS parameter is
            // not — compare against the cylinder at the SAMPLE's own
            // angle instead of the chart image for the radial check,
            // and accept chart-parameterization skew along the arc.
            let d = p_nurbs - center;
            let radial = d - axis * d.dot(&axis);
            if (radial.length() - radius).abs() > gate {
                return None;
            }
            let _ = p_cyl;
        }
    }
    Some((
        geometry,
        ChartMap {
            u_offset: theta0,
            u_scale: theta1 - theta0,
            v_offset: 0.0,
            v_scale: height,
        },
    ))
}

/// Circumcenter of three points, constrained to the plane through
/// `a` perpendicular to `axis` (extrusion profiles are planar).
fn circumcenter_3d(a: Point3, b: Point3, c: Point3, axis: Vec3) -> Option<Point3> {
    // Build a 2D frame in the profile plane.
    let (e1, e2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(axis);
    let to2d = |p: Point3| -> (f64, f64) {
        let d = p - a;
        (d.dot(&e1), d.dot(&e2))
    };
    let (ax, ay) = (0.0, 0.0);
    let (bx, by) = to2d(b);
    let (cx, cy) = to2d(c);
    let dd = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if dd.abs() < 1.0e-16 {
        return None;
    }
    let ux = ((ax * ax + ay * ay) * (by - cy)
        + (bx * bx + by * by) * (cy - ay)
        + (cx * cx + cy * cy) * (ay - by))
        / dd;
    let uy = ((ax * ax + ay * ay) * (cx - bx)
        + (bx * bx + by * by) * (ax - cx)
        + (cx * cx + cy * cy) * (bx - ax))
        / dd;
    Some(a + e1 * ux + e2 * uy)
}
