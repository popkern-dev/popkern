//! Plane × sphere general dispatch — covers any plane orientation.
//!
//! The polar-axis dispatcher in `analytic_trig_plane_sphere.rs`
//! handles the special case where the plane normal is parallel to
//! `(0, 0, 1)` (popkern's sphere polar axis), giving an exact
//! analytic `Pcurve::Linear` on the sphere side. For arbitrary plane
//! orientations, the sphere pcurve traces a non-trivial path in
//! `(longitude, colatitude)` coordinates that no current pcurve
//! variant captures exactly — we sample-and-fit a NURBS pcurve
//! instead.
//!
//! Output:
//! - `curve_3d`: exact `CurveGeometry::Arc` (every plane-sphere
//!   intersection is a circle in 3-D).
//! - pcurve on plane: exact `Pcurve::Elliptical` (projection of the
//!   circle into plane UV).
//! - pcurve on sphere: degree-5 NURBS interpolation over 33 samples
//!   in `[θ_start, θ_end]`. ~1e-12 residual; well within
//!   `tol.linear()`.
//! - `ParameterMap::Identity` holds via shared θ.

use crate::geometry::curve::CurveGeometry;
use crate::geometry::edge_triple::analytic_trig_helpers::decompose_2x2;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, SurfaceSide};
use crate::geometry::edge_triple::types::{
    AnalyticPair, BranchSelector, EdgeTopology, EdgeTriple, EdgeTripleProvenance,
    EdgeTripleQuality, ParameterMap,
};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Plane × sphere intersection for any plane orientation. Returns
/// `None` for the tangent / non-intersecting cases (plane farther
/// from sphere centre than the sphere's radius). The polar dispatcher
/// is preferred when applicable; this one is the catch-all.
#[allow(clippy::too_many_arguments)]
pub(super) fn try_plane_sphere_general(
    plane_origin: Point3,
    plane_normal: Vec3,
    plane_u: Vec3,
    plane_v: Vec3,
    plane_is_b: bool,
    sphere_center: Point3,
    sphere_radius: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let n_len_sq = plane_normal.dot(&plane_normal);
    if n_len_sq < tol.solver() {
        return None;
    }
    // Signed distance from sphere centre to the plane along the
    // plane normal.
    let h = (sphere_center - plane_origin).dot(&plane_normal) / n_len_sq.sqrt();
    if h.abs() >= sphere_radius - tol.linear() {
        return None; // tangent or non-intersecting
    }
    let circle_radius = (sphere_radius * sphere_radius - h * h).sqrt();
    let n_unit = plane_normal * (1.0 / n_len_sq.sqrt());
    let circle_center = sphere_center + n_unit * h;

    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(n_unit);

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };

    let theta_of = |p: Point3, side: SurfaceSide| -> Result<f64, EdgeTripleError> {
        let dp = p - circle_center;
        let cu = dp.dot(&arc_basis_u);
        let cv = dp.dot(&arc_basis_v);
        let r = (cu * cu + cv * cv).sqrt();
        if (r - circle_radius).abs() > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: r - circle_radius,
            });
        }
        Ok(cv.atan2(cu))
    };
    let theta_start = match theta_of(start, SurfaceSide::A) {
        Ok(t) => t,
        Err(e) => return Some(Err(e)),
    };
    let theta_end = match theta_of(end, SurfaceSide::B) {
        Ok(t) => t,
        Err(e) => return Some(Err(e)),
    };
    let theta_end = if theta_end < theta_start {
        theta_end + std::f64::consts::TAU
    } else {
        theta_end
    };

    // Curve_3d: an exact analytic Arc on the cutting plane.
    let curve_3d = CurveGeometry::Arc {
        center: circle_center,
        normal: n_unit,
        radius: circle_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    // Plane pcurve: exact Elliptical via 2×2 SVD (same approach as
    // perpendicular plane×cyl). Both u_pcurve(θ) and v_pcurve(θ) are
    // linear-in-(cos θ, sin θ).
    let center_in_plane = circle_center - plane_origin;
    let plane_u_len_sq = plane_u.dot(&plane_u);
    let plane_v_len_sq = plane_v.dot(&plane_v);
    let center_u = center_in_plane.dot(&plane_u) / plane_u_len_sq;
    let center_v = center_in_plane.dot(&plane_v) / plane_v_len_sq;
    let a_uu = circle_radius * arc_basis_u.dot(&plane_u) / plane_u_len_sq;
    let a_vu = circle_radius * arc_basis_u.dot(&plane_v) / plane_v_len_sq;
    let a_uv = circle_radius * arc_basis_v.dot(&plane_u) / plane_u_len_sq;
    let a_vv = circle_radius * arc_basis_v.dot(&plane_v) / plane_v_len_sq;
    let (u_radius, v_radius, rotation, phase) = decompose_2x2(a_uu, a_uv, a_vu, a_vv);
    let pcurve_plane = Pcurve::elliptical(
        (center_u, center_v),
        u_radius,
        v_radius,
        rotation,
        theta_start - phase,
        theta_end - phase,
        theta_start,
        theta_end,
    );

    // Sphere pcurve: NURBS sample-and-fit over θ. The (longitude,
    // colatitude) trace is `lon(θ) = atan2(y, x)` and `colat(θ) =
    // acos(z / r)` evaluated on the 3-D arc point — non-trivial
    // functions of θ that no current Pcurve variant captures.
    let n_samples = 33;
    let mut sphere_pts: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    let span = theta_end - theta_start;
    for i in 0..n_samples {
        let frac = i as f64 / (n_samples as f64 - 1.0);
        let theta = theta_start + frac * span;
        let p3 = circle_center
            + arc_basis_u * (circle_radius * theta.cos())
            + arc_basis_v * (circle_radius * theta.sin());
        let dp = p3 - sphere_center;
        let lon = dp.y.atan2(dp.x);
        let colat = (dp.z / sphere_radius).clamp(-1.0, 1.0).acos();
        sphere_pts.push(Point3::new(lon, colat, 0.0));
        params.push(frac);
    }
    // Unwrap longitude jumps so the NURBS-fit doesn't smooth over a
    // 2π discontinuity at the seam.
    for i in 1..n_samples {
        let du = sphere_pts[i].x - sphere_pts[i - 1].x;
        if du > std::f64::consts::PI {
            for cp in sphere_pts.iter_mut().skip(i) {
                cp.x -= std::f64::consts::TAU;
            }
        } else if du < -std::f64::consts::PI {
            for cp in sphere_pts.iter_mut().skip(i) {
                cp.x += std::f64::consts::TAU;
            }
        }
    }

    let sphere_nurbs = match crate::geometry::fitting::fit_curve(&sphere_pts, &params, 5, n_samples)
    {
        Ok(c) => c,
        Err(_) => return None,
    };
    let sphere_nurbs = reknot_to(sphere_nurbs, theta_start, theta_end);
    let pcurve_sphere = Pcurve::from_nurbs(sphere_nurbs);

    let plane_surf = SurfaceGeometry::Plane {
        origin: plane_origin,
        normal: plane_normal,
        u_axis: plane_u,
        v_axis: plane_v,
    };
    let sphere_surf = SurfaceGeometry::Sphere {
        center: sphere_center,
        radius: sphere_radius,
    };

    let (pcurve_a, pcurve_b, surf_a, surf_b) = if plane_is_b {
        (pcurve_sphere, pcurve_plane, &sphere_surf, &plane_surf)
    } else {
        (pcurve_plane, pcurve_sphere, &plane_surf, &sphere_surf)
    };

    let triple = EdgeTriple {
        curve_3d,
        pcurve_a,
        pcurve_b,
        t_start: theta_start,
        t_end: theta_end,
        topology: EdgeTopology::Open,
        parameter_map: ParameterMap::Identity,
        provenance: EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneSphere,
            branch: BranchSelector::OnlyBranch,
            cardinality: None,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}

fn reknot_to(
    mut curve: crate::geometry::nurbs_curve::NurbsCurve,
    t0: f64,
    t1: f64,
) -> crate::geometry::nurbs_curve::NurbsCurve {
    let Some(&k_lo) = curve.knots.first() else {
        return curve;
    };
    let Some(&k_hi) = curve.knots.last() else {
        return curve;
    };
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
