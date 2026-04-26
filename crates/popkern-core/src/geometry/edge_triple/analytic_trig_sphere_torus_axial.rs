//! Sphere × torus axial dispatch — octic factors into circles.
//!
//! When the sphere centre lies on the torus axis, the algebraic
//! intersection (generically a degree-8 curve from a degree-2
//! sphere × degree-4 torus) **factors into 0, 1, or 2 horizontal
//! circles**, each at a fixed axial height. This is one of the
//! cleanest factorisations in the quartic-pair space.
//!
//! Math: write `P` on the torus as
//!
//!   P = T_c + (R + r·cos v)·(cos u·p1 + sin u·p2) + r·sin v·a
//!
//! and the sphere centre as `S_c = T_c + h·a` (axis-collinear).
//! Then `P − S_c` has axial component `r·sin v − h` and radial
//! component `R + r·cos v`, so
//!
//!   |P − S_c|² = (R + r·cos v)² + (r·sin v − h)²
//!
//! depends only on `v`. Setting equal to `ρ²` gives the constraint
//!
//!   R·cos v − h·sin v = (ρ² − R² − r² − h²) / (2r)
//!
//! which is a one-amplitude phase-shifted cosine: real solutions iff
//! the right-hand side fits within `±√(R² + h²)`. When it does,
//! `v_± = −ψ ± acos(rhs / √(R²+h²))` where `ψ = atan2(h, R)`. Each
//! `v_±` defines an axial-height circle of radius `R + r·cos v_±`.
//!
//! Output: exact `CurveGeometry::Arc` plus `Pcurve::Linear` on the
//! torus side. Sphere pcurve is `Pcurve::Linear` when the torus
//! axis is parallel to the sphere's polar axis `(0, 0, 1)`; otherwise
//! NURBS sample-and-fit.

use crate::geometry::Point3;
use crate::geometry::Vec3;
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
use crate::tolerance::Tolerance;

#[allow(clippy::too_many_arguments)]
pub(super) fn try_sphere_torus_axial(
    sphere_center: Point3,
    sphere_radius: f64,
    sphere_is_b: bool,
    torus_center: Point3,
    torus_axis: Vec3,
    major_radius: f64,
    minor_radius: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    // Detection: sphere centre on torus axis.
    let dp = sphere_center - torus_center;
    let h = dp.dot(&torus_axis);
    let perp = dp - torus_axis * h;
    if perp.length() > tol.linear() {
        return None;
    }

    // Solve R·cos v − h·sin v = C / (2r) for v, with
    // C = ρ² − R² − r² − h².
    let big_r = major_radius;
    let r = minor_radius;
    let rho = sphere_radius;
    let c = rho * rho - big_r * big_r - r * r - h * h;
    let amp = (big_r * big_r + h * h).sqrt();
    let rhs = c / (2.0 * r * amp);
    if !(-1.0..=1.0).contains(&rhs) {
        return None; // no real intersection (sphere too big or too small)
    }
    let psi = h.atan2(big_r); // tan ψ = h / R
    let arc_acos = rhs.acos();
    let v_plus = -psi + arc_acos;
    let v_minus = -psi - arc_acos;

    // Pick branch by start-anchor's v.
    let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(torus_axis);
    let solve_v = |p: Point3| -> f64 {
        let dp = p - torus_center;
        let v_axis = dp.dot(&torus_axis);
        let perp_vec = dp - torus_axis * v_axis;
        let tube_r = perp_vec.length();
        let sin_v = (v_axis / r).clamp(-1.0, 1.0);
        let cos_v = ((tube_r - big_r) / r).clamp(-1.0, 1.0);
        sin_v.atan2(cos_v)
    };
    let dist_v = |a: f64, b: f64| -> f64 {
        let d = (a - b).abs();
        d.min(std::f64::consts::TAU - d)
    };
    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };
    let v_start_truth = solve_v(start);
    let v_branch = if dist_v(v_plus, v_start_truth) <= dist_v(v_minus, v_start_truth) {
        v_plus
    } else {
        v_minus
    };
    let circle_radius = big_r + r * v_branch.cos();
    if circle_radius <= tol.linear() {
        return None; // degenerate (zero-radius circle)
    }
    let axial_height = r * v_branch.sin();
    let circle_center = torus_center + torus_axis * axial_height;

    // Curve is a circle in the plane perpendicular to torus axis
    // at the chosen axial height.
    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(torus_axis);
    let theta_of = |p: Point3, side: SurfaceSide| -> Result<f64, EdgeTripleError> {
        let d = p - circle_center;
        let cu = d.dot(&arc_basis_u);
        let cv = d.dot(&arc_basis_v);
        let r_real = (cu * cu + cv * cv).sqrt();
        if (r_real - circle_radius).abs() > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: r_real - circle_radius,
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
        normal: torus_axis,
        radius: circle_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    // Torus pcurve: Linear (constant v_branch, linear u with frame
    // phase offset). Same machinery as `try_plane_torus_perpendicular`.
    let torus_u_offset = arc_basis_u.dot(&perp2).atan2(arc_basis_u.dot(&perp1));
    let v_branch_canon = if v_branch < 0.0 {
        v_branch + std::f64::consts::TAU
    } else {
        v_branch
    };
    let pcurve_torus = Pcurve::linear(
        (theta_start + torus_u_offset, v_branch_canon),
        (theta_end + torus_u_offset, v_branch_canon),
        theta_start,
        theta_end,
    );

    // Sphere pcurve: when the torus axis is parallel to the sphere's
    // polar axis (0, 0, 1), the colatitude is constant and the
    // longitude is linear in θ — Linear pcurve. Otherwise, sample-
    // and-fit a NURBS pcurve.
    let polar = Vec3::new(0.0, 0.0, 1.0);
    let pcurve_sphere = if torus_axis.dot(&polar).abs() >= 1.0 - tol.angular() {
        // Linear case. dp.z constant = circle_center.z − sphere_center.z.
        let z_offset = (circle_center.z - sphere_center.z) / sphere_radius;
        let colat = z_offset.clamp(-1.0, 1.0).acos();
        let global_x = Vec3::new(1.0, 0.0, 0.0);
        let global_y = Vec3::new(0.0, 1.0, 0.0);
        let sphere_u_offset = arc_basis_u.dot(&global_y).atan2(arc_basis_u.dot(&global_x));
        Pcurve::linear(
            (theta_start + sphere_u_offset, colat),
            (theta_end + sphere_u_offset, colat),
            theta_start,
            theta_end,
        )
    } else {
        // NURBS sample-and-fit fallback for off-polar torus axis.
        let n_samples = 65;
        let mut pts: Vec<Point3> = Vec::with_capacity(n_samples);
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
            pts.push(Point3::new(lon, colat, 0.0));
            params.push(frac);
        }
        // Unwrap longitude.
        for i in 1..n_samples {
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
        let nurbs = match crate::geometry::fitting::fit_curve(&pts, &params, 5, n_samples) {
            Ok(c) => c,
            Err(_) => return None,
        };
        Pcurve::from_nurbs(reknot_to(nurbs, theta_start, theta_end))
    };

    let _ = decompose_2x2; // currently unused — sphere pcurve uses Linear/NURBS only

    let sphere_surf = SurfaceGeometry::Sphere {
        center: sphere_center,
        radius: sphere_radius,
    };
    let torus_surf = SurfaceGeometry::Torus {
        center: torus_center,
        axis: torus_axis,
        major_radius,
        minor_radius,
    };

    let (pcurve_a, pcurve_b, surf_a, surf_b) = if sphere_is_b {
        (pcurve_torus, pcurve_sphere, &torus_surf, &sphere_surf)
    } else {
        (pcurve_sphere, pcurve_torus, &sphere_surf, &torus_surf)
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
            kind: AnalyticPair::SphereTorus,
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
