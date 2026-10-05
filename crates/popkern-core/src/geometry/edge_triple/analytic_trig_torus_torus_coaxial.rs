//! Torus × torus coaxial dispatch — octic factors into circles.
//!
//! When two tori share an axis (parallel, with one centre lying on
//! the other's axis), the algebraic intersection (degree-4 × degree-4
//! = degree-8 generally) collapses to ≤ 2 horizontal circles by the
//! axisymmetric-collapse argument.
//!
//! Math: in the meridian half-plane `(X, Y)` (where `X` = perp
//! distance from shared axis, `Y` = signed axial distance from
//! `T_a`'s centre along `axis_a`), each torus tube traces a circle:
//!
//!   torus A:  (X − R_a)² + Y² = r_a²        centre (R_a, 0)
//!   torus B:  (X − R_b)² + (Y − h)² = r_b²  centre (R_b, h)
//!
//! where `h = (T_b − T_a) · axis_a` (signed axial offset). Two-circle
//! intersection in the plane gives 0, 1, or 2 `(X, Y)` solutions,
//! each defining a horizontal circle of radius `X` at axial height
//! `Y` from `T_a`'s centre.
//!
//! Output: exact `CurveGeometry::Arc` + `Pcurve::Linear` on both
//! surfaces (constant v on each, linear u).
//!
//! Axis-sign handling: when `axis_b · axis_a = -1`, the v parameter
//! on torus B reads `Y_b = -(Y - h)` (axial offset along `axis_b`).
//! Frame phase is captured per surface via `orthonormal_basis(axis)`,
//! matching the cone × torus axial pattern.

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
use popkern_math::PortableMath;

#[allow(clippy::too_many_arguments)]
pub(super) fn try_torus_torus_coaxial(
    center_a: Point3,
    axis_a: Vec3,
    big_r_a: f64,
    r_a: f64,
    center_b: Point3,
    axis_b: Vec3,
    big_r_b: f64,
    r_b: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    // Detection: axes parallel.
    let axis_align = axis_a.dot(&axis_b);
    if axis_align.abs() < 1.0 - tol.angular() {
        return None;
    }
    // T_b's centre lies on T_a's axis.
    let dp = center_b - center_a;
    let perp = dp - axis_a * dp.dot(&axis_a);
    if perp.length() > tol.linear() {
        return None;
    }

    let h_axis = dp.dot(&axis_a); // signed
    // Two-circle intersection in (X, Y) meridian plane.
    // Circle A: centre (R_a, 0), radius r_a.
    // Circle B: centre (R_b, h_axis), radius r_b.
    let dx = big_r_b - big_r_a;
    let dy = h_axis;
    let dist_sq = dx * dx + dy * dy;
    let dist = dist_sq.sqrt();
    // Coincident-centre degenerate: same tubes (no curve) or nested
    // tubes (no intersection) — defer to NURBS marching either way.
    if dist < tol.linear() {
        return None;
    }
    let r_sum = r_a + r_b;
    let r_dif = (r_a - r_b).abs();
    if dist > r_sum + tol.linear() || dist < r_dif - tol.linear() {
        return None; // tubes disjoint or one inside the other
    }
    // Standard two-circle intersection formula.
    //   a = (d² + r_a² − r_b²) / (2d)
    //   h_off = √(r_a² − a²)
    let a = (dist_sq + r_a * r_a - r_b * r_b) / (2.0 * dist);
    let disc = r_a * r_a - a * a;
    let h_off = disc.max(0.0).sqrt();
    let ux = dx / dist;
    let uy = dy / dist;
    // Perpendicular (rotate +90° in (X, Y) plane).
    let nx = -uy;
    let ny = ux;
    let mid_x = big_r_a + a * ux;
    let mid_y = a * uy;
    let sol_plus = (mid_x + h_off * nx, mid_y + h_off * ny);
    let sol_minus = (mid_x - h_off * nx, mid_y - h_off * ny);

    // Branch selection by anchor's (X, Y) projection.
    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };
    let project_xy = |p: Point3| -> (f64, f64) {
        let dp = p - center_a;
        let y = dp.dot(&axis_a);
        let perp = dp - axis_a * y;
        let x = perp.length();
        (x, y)
    };
    let (xs, ys) = project_xy(start);
    let dist_sq_to = |sol: (f64, f64)| -> f64 {
        let dx = sol.0 - xs;
        let dy = sol.1 - ys;
        dx * dx + dy * dy
    };
    let (sol_x, sol_y) = if dist_sq_to(sol_plus) <= dist_sq_to(sol_minus) {
        sol_plus
    } else {
        sol_minus
    };
    let circle_radius = sol_x;
    if circle_radius <= tol.linear() {
        return None;
    }
    let circle_center = center_a + axis_a * sol_y;

    // v on each torus: cos v = (X − R) / r, sin v = Y_local / r.
    let sign_b = axis_align.signum();
    let cos_v_a = ((sol_x - big_r_a) / r_a).clamp(-1.0, 1.0);
    let sin_v_a = (sol_y / r_a).clamp(-1.0, 1.0);
    let v_a = sin_v_a.patan2(cos_v_a);
    let y_b_local = sign_b * (sol_y - h_axis);
    let cos_v_b = ((sol_x - big_r_b) / r_b).clamp(-1.0, 1.0);
    let sin_v_b = (y_b_local / r_b).clamp(-1.0, 1.0);
    let v_b = sin_v_b.patan2(cos_v_b);
    let canon = |v: f64| -> f64 {
        if v < 0.0 {
            v + std::f64::consts::TAU
        } else {
            v
        }
    };
    let v_a_canon = canon(v_a);
    let v_b_canon = canon(v_b);

    let (perp1_a, perp2_a) = crate::geometry::analytic_surfaces::orthonormal_basis(axis_a);
    let (perp1_b, perp2_b) = crate::geometry::analytic_surfaces::orthonormal_basis(axis_b);
    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(axis_a);
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
        Ok(cv.patan2(cu))
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
        normal: axis_a,
        radius: circle_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    let a_u_offset = arc_basis_u.dot(&perp2_a).patan2(arc_basis_u.dot(&perp1_a));
    let pcurve_a = Pcurve::linear(
        (theta_start + a_u_offset, v_a_canon),
        (theta_end + a_u_offset, v_a_canon),
        theta_start,
        theta_end,
    );
    let b_u_offset = arc_basis_u.dot(&perp2_b).patan2(arc_basis_u.dot(&perp1_b));
    let pcurve_b = Pcurve::linear(
        (theta_start + b_u_offset, v_b_canon),
        (theta_end + b_u_offset, v_b_canon),
        theta_start,
        theta_end,
    );

    let surf_a = SurfaceGeometry::Torus {
        center: center_a,
        axis: axis_a,
        major_radius: big_r_a,
        minor_radius: r_a,
    };
    let surf_b = SurfaceGeometry::Torus {
        center: center_b,
        axis: axis_b,
        major_radius: big_r_b,
        minor_radius: r_b,
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
            kind: AnalyticPair::TorusTorus,
            branch: BranchSelector::OnlyBranch,
            cardinality: None,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, &surf_a, &surf_b, tol,
    ))
}
