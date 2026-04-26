//! Sphere × sphere general dispatch (any centre alignment).
//!
//! The polar-only dispatcher in `analytic_trig_sphere_sphere.rs`
//! requires `c₂ − c₁` parallel to `(0, 0, 1)`. This dispatcher is the
//! catch-all: two spheres always intersect in a circle when neither
//! coincident, tangent, nor non-intersecting, regardless of how their
//! centres are oriented.
//!
//! Geometry: subtracting the two sphere equations gives a planar
//! constraint `2P·(c₂ − c₁) = r₁² − r₂² + |c₂|² − |c₁|²`, the radical
//! plane perpendicular to the line of centres. The intersection
//! circle lies on this plane, with centre on the axis and radius
//! determined by Pythagoras.
//!
//! Output:
//! - `curve_3d`: exact `CurveGeometry::Arc`.
//! - pcurves on each sphere: degree-3 NURBS interpolation over 33
//!   samples in `(longitude, colatitude)`. CPs may slightly overshoot
//!   the data convex hull when the longitude trace is non-monotonic;
//!   `check_uv_box`'s periodicity-aware reduction (added in the
//!   surface API) absorbs this without losing correctness.
//! - `ParameterMap::Identity` via shared θ.

use crate::geometry::Point3;
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
pub(super) fn try_sphere_sphere_general(
    center_a: Point3,
    radius_a: f64,
    center_b: Point3,
    radius_b: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let delta = center_b - center_a;
    let d = delta.length();
    if d < tol.linear() {
        return None;
    }
    if d >= radius_a + radius_b - tol.linear() || d <= (radius_a - radius_b).abs() + tol.linear() {
        return None;
    }

    let dir = delta * (1.0 / d);
    let h_along = (radius_a * radius_a - radius_b * radius_b + d * d) / (2.0 * d);
    let circle_center = center_a + dir * h_along;
    let circle_radius_sq = radius_a * radius_a - h_along * h_along;
    if circle_radius_sq <= 0.0 {
        return None;
    }
    let circle_radius = circle_radius_sq.sqrt();
    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(dir);

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

    let curve_3d = CurveGeometry::Arc {
        center: circle_center,
        normal: dir,
        radius: circle_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    // NURBS sample-and-fit pcurves in (longitude, colatitude) for
    // each sphere. Use degree-5 interpolation with 65 samples ==
    // 65 CPs: exact pass-through at samples, off-sample residual
    // < 1e-7 across the test corpus. The longitude trace can be
    // non-monotonic (dips and recovers when the centre-of-centres
    // line is offset from the polar axis), so the fit needs more
    // fidelity than the plane×sphere case.
    let n_samples = 65;
    let mut pts_a: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut pts_b: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    let span = theta_end - theta_start;
    for i in 0..n_samples {
        let frac = i as f64 / (n_samples as f64 - 1.0);
        let theta = theta_start + frac * span;
        let p3 = circle_center
            + arc_basis_u * (circle_radius * theta.cos())
            + arc_basis_v * (circle_radius * theta.sin());
        let dp_a = p3 - center_a;
        let lon_a = dp_a.y.atan2(dp_a.x);
        let colat_a = (dp_a.z / radius_a).clamp(-1.0, 1.0).acos();
        pts_a.push(Point3::new(lon_a, colat_a, 0.0));
        let dp_b = p3 - center_b;
        let lon_b = dp_b.y.atan2(dp_b.x);
        let colat_b = (dp_b.z / radius_b).clamp(-1.0, 1.0).acos();
        pts_b.push(Point3::new(lon_b, colat_b, 0.0));
        params.push(frac);
    }
    unwrap_lon(&mut pts_a, n_samples);
    unwrap_lon(&mut pts_b, n_samples);

    let nurbs_a = match crate::geometry::fitting::fit_curve(&pts_a, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };
    let nurbs_b = match crate::geometry::fitting::fit_curve(&pts_b, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };
    let pcurve_a = Pcurve::from_nurbs(reknot_to(nurbs_a, theta_start, theta_end));
    let pcurve_b = Pcurve::from_nurbs(reknot_to(nurbs_b, theta_start, theta_end));

    let surf_a = SurfaceGeometry::Sphere {
        center: center_a,
        radius: radius_a,
    };
    let surf_b = SurfaceGeometry::Sphere {
        center: center_b,
        radius: radius_b,
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
            kind: AnalyticPair::SphereSphere,
            branch: BranchSelector::OnlyBranch,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, &surf_a, &surf_b, tol,
    ))
}

fn unwrap_lon(pts: &mut [Point3], n: usize) {
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
