//! Cylinder × cylinder perpendicular-axes dispatch (4th-order case).
//!
//! Two cylinders with perpendicular intersecting axes meet in a
//! quartic curve. When the radii match, this is **Viviani's curve**
//! (the intersection of a sphere and a tangent cylinder, equivalently
//! the intersection of two equal perpendicular cylinders sharing an
//! origin). For unequal radii, it's a generic quartic with up to 4
//! disjoint branches.
//!
//! Unlike the analytic-output dispatchers earlier in the matrix,
//! this one produces NURBS sample-and-fit output — the curve is
//! genuinely 4th-degree algebraic and has no parameterisation in
//! Line/Arc/Ellipse/Parabola/Hyperbola form. The marching-style
//! dispatcher's value over the legacy `intersect_surfaces` fallback
//! is:
//!
//! 1. **Anchor-driven branch selection** — picks the correct one of
//!    up to 4 branches deterministically from the start anchor.
//! 2. **Closed-form sampling** — uses the implicit-equation pair to
//!    generate samples directly (no Newton iteration), which is
//!    faster and more robust than generic Newton-marching.
//! 3. **Direction conformance** — guarantees `ParameterMap::Identity`
//!    over the cyl A angular parameter `u`.
//!
//! Restriction: requires the two cylinder axes to intersect (closest
//! distance ≤ `tol.linear()`). Non-intersecting perpendicular pairs
//! fall through to the legacy path.

use crate::geometry::Point3;
use crate::geometry::Vec3;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, SurfaceSide};
use crate::geometry::edge_triple::types::{
    AnalyticPair, BranchSelector, EdgeTopology, EdgeTriple, EdgeTripleProvenance,
    EdgeTripleQuality, ParameterMap,
};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

#[allow(clippy::too_many_arguments)]
pub(super) fn try_cyl_cyl_perpendicular(
    origin_a: Point3,
    axis_a: Vec3,
    radius_a: f64,
    origin_b: Point3,
    axis_b: Vec3,
    radius_b: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    // Detection: axes perpendicular.
    if axis_a.dot(&axis_b).abs() > tol.angular() {
        return None;
    }
    // Axes must intersect (closest distance ≤ tol.linear()). For
    // perpendicular axes, the closest-points line is along
    // `axis_a × axis_b`. If `(origin_b − origin_a)` has a non-zero
    // component along that direction, the axes are skew.
    let cross = axis_a.cross(&axis_b);
    let cross_len = cross.length();
    if cross_len < tol.solver() {
        return None;
    }
    let cross_unit = cross * (1.0 / cross_len);
    let skew = (origin_b - origin_a).dot(&cross_unit).abs();
    if skew > tol.linear() {
        return None;
    }

    // Find the intersection point of the axes by projecting origin_b
    // onto axis_a's perpendicular. Equivalently: solve for `t_a`,
    // `t_b` such that `origin_a + t_a · axis_a = origin_b + t_b ·
    // axis_b`. With perpendicular axes:
    //   t_a = (origin_b − origin_a) · axis_a
    //   t_b = (origin_a − origin_b) · axis_b
    let t_a = (origin_b - origin_a).dot(&axis_a);
    let t_b = (origin_a - origin_b).dot(&axis_b);
    let intersection = origin_a + axis_a * t_a;
    // Sanity: intersection should equal origin_b + t_b · axis_b.
    let alt = origin_b + axis_b * t_b;
    if intersection.distance_to(&alt) > tol.linear() {
        return None;
    }

    // Local frame at intersection: axis_a as "z", axis_b as "x", and
    // their cross as "y" (= cross_unit, signed).
    let z_dir = axis_a;
    let x_dir = axis_b;
    let y_dir = z_dir.cross(&x_dir).normalized().unwrap_or(Vec3::zero());
    if y_dir.length() < tol.solver() {
        return None;
    }

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };

    // Local coords of an anchor.
    let local = |p: Point3| -> (f64, f64, f64) {
        let dp = p - intersection;
        (dp.dot(&x_dir), dp.dot(&y_dir), dp.dot(&z_dir))
    };
    let (sx, sy, sz) = local(start);
    let (ex, ey, ez) = local(end);

    // Cyl A (axis = z): x² + y² = r_a². So sx² + sy² = r_a², ex² + ey² = r_a².
    // Cyl B (axis = x): y² + z² = r_b². So sy² + sz² = r_b², ey² + ez² = r_b².
    // Validate anchors.
    let on_a =
        |x: f64, y: f64| -> bool { (x * x + y * y - radius_a * radius_a).abs() < tol.linear() };
    let on_b =
        |y: f64, z: f64| -> bool { (y * y + z * z - radius_b * radius_b).abs() < tol.linear() };
    if !on_a(sx, sy) || !on_b(sy, sz) {
        return Some(Err(EdgeTripleError::AnchorOffSurface {
            which: SurfaceSide::A,
            distance: ((sx * sx + sy * sy).sqrt() - radius_a)
                .abs()
                .max(((sy * sy + sz * sz).sqrt() - radius_b).abs()),
        }));
    }
    if !on_a(ex, ey) || !on_b(ey, ez) {
        return Some(Err(EdgeTripleError::AnchorOffSurface {
            which: SurfaceSide::B,
            distance: ((ex * ex + ey * ey).sqrt() - radius_a)
                .abs()
                .max(((ey * ey + ez * ez).sqrt() - radius_b).abs()),
        }));
    }

    // Parameterise by cyl A's u-angle: u = atan2(y_local, x_local).
    let u_start = sy.atan2(sx);
    let u_end_raw = ey.atan2(ex);
    let u_end = if u_end_raw < u_start {
        u_end_raw + std::f64::consts::TAU
    } else {
        u_end_raw
    };

    // Branch selection: at parameter u, x = r_a cos u, y = r_a sin u.
    // From cyl B: z² = r_b² − y². Sign of z stays consistent across
    // the curve segment between two anchors as long as we don't cross
    // z = 0. Match the start anchor's sign.
    let z_sign = sz.signum();
    if z_sign.abs() < tol.solver() {
        return None; // start anchor on z=0 — branch ambiguous
    }

    // Sample, build NURBS curve_3d in 3-D and pcurves on each cyl.
    let n_samples = 65;
    let mut pts_3d: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut pts_a: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut pts_b: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    let span = u_end - u_start;
    for i in 0..n_samples {
        let frac = i as f64 / (n_samples as f64 - 1.0);
        let u = u_start + frac * span;
        let x = radius_a * u.cos();
        let y = radius_a * u.sin();
        let z_sq = radius_b * radius_b - y * y;
        if z_sq < -tol.linear() {
            return None; // segment crosses out of the real-z region
        }
        let z = z_sq.max(0.0).sqrt() * z_sign;
        let p3 = intersection + x_dir * x + y_dir * y + z_dir * z;
        pts_3d.push(p3);
        // Cyl A pcurve: u_cyl_a, v_cyl_a. cyl_a's u = atan2(y, x) on
        // the local frame, but cyl_a's surface uses
        // orthonormal_basis(axis_a) which may differ. Compute
        // surface uv from p3.
        let (u_a, v_a) = cyl_uv(p3, origin_a, axis_a);
        pts_a.push(Point3::new(u_a, v_a, 0.0));
        let (u_b, v_b) = cyl_uv(p3, origin_b, axis_b);
        pts_b.push(Point3::new(u_b, v_b, 0.0));
        params.push(frac);
    }
    // Unwrap atan2 jumps in the cylinder u-coords.
    unwrap_u(&mut pts_a, n_samples);
    unwrap_u(&mut pts_b, n_samples);

    let mut nurbs_3d = match crate::geometry::fitting::fit_curve(&pts_3d, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };
    // G1: snap first/last CPs to the exact anchors. fit_curve uses
    // clamped boundary conditions so they're already approximately
    // equal, but FP rounding in the trigonometric sampling can
    // introduce 1-ULP differences that trip equality checks.
    if let Some(first) = nurbs_3d.control_points.first_mut() {
        *first = start;
    }
    if let Some(last) = nurbs_3d.control_points.last_mut() {
        *last = end;
    }
    let nurbs_a = match crate::geometry::fitting::fit_curve(&pts_a, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };
    let nurbs_b = match crate::geometry::fitting::fit_curve(&pts_b, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };
    let curve_3d = CurveGeometry::Nurbs(reknot_to(nurbs_3d, u_start, u_end));
    let _ = (start, end); // silence warnings about move earlier
    let pcurve_a = Pcurve::from_nurbs(reknot_to(nurbs_a, u_start, u_end));
    let pcurve_b = Pcurve::from_nurbs(reknot_to(nurbs_b, u_start, u_end));

    let cyl_a_surf = SurfaceGeometry::Cylinder {
        origin: origin_a,
        axis: axis_a,
        radius: radius_a,
    };
    let cyl_b_surf = SurfaceGeometry::Cylinder {
        origin: origin_b,
        axis: axis_b,
        radius: radius_b,
    };

    let triple = EdgeTriple {
        curve_3d,
        pcurve_a,
        pcurve_b,
        t_start: u_start,
        t_end: u_end,
        topology: EdgeTopology::Open,
        parameter_map: ParameterMap::Identity,
        provenance: EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::CylinderCylinder,
            branch: BranchSelector::OnlyBranch,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple,
        &cyl_a_surf,
        &cyl_b_surf,
        tol,
    ))
}

fn cyl_uv(p: Point3, origin: Point3, axis: Vec3) -> (f64, f64) {
    let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(axis);
    let dp = p - origin;
    let v = dp.dot(&axis);
    let perp_vec = dp - axis * v;
    let u = perp_vec.dot(&perp2).atan2(perp_vec.dot(&perp1));
    (u, v)
}

fn unwrap_u(pts: &mut [Point3], n: usize) {
    for i in 1..n {
        let du = pts[i].x - pts[i - 1].x;
        if du > std::f64::consts::PI {
            for cp in pts.iter_mut().skip(i) {
                cp.x -= std::f64::consts::TAU;
            }
        } else if du < -std::f64::consts::PI {
            for cp in pts.iter_mut().skip(i) {
                cp.x += std::f64::consts::TAU;
            }
        }
    }
}

fn reknot_to(
    mut curve: crate::geometry::nurbs_curve::NurbsCurve,
    t0: f64,
    t1: f64,
) -> crate::geometry::nurbs_curve::NurbsCurve {
    if curve.knots.is_empty() {
        return curve;
    }
    let k_lo = *curve.knots.first().unwrap();
    let k_hi = *curve.knots.last().unwrap();
    let span_old = k_hi - k_lo;
    if span_old.abs() < 1.0e-30 {
        return curve;
    }
    let span_new = t1 - t0;
    for k in &mut curve.knots {
        *k = t0 + (*k - k_lo) / span_old * span_new;
    }
    curve
}
