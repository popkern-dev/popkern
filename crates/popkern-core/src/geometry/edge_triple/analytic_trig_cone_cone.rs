//! Cone × cone coaxial dispatch.
//!
//! When two cones share an axis (parallel directions, both apexes
//! lying on the shared axis) and have different half-angles, the
//! surfaces meet in exactly one circle at the height where their
//! radii match:
//!
//!   v_a · tan(θ_a) = v_b · tan(θ_b),    v_a − v_b = d
//!
//! Solving:
//!
//!   v_a = d · tan(θ_b) / (tan(θ_b) − tan(θ_a))
//!
//! When the half-angles match: surfaces either coincide (`d = 0`)
//! or are parallel-and-non-intersecting; both return `None`.
//!
//! Output:
//! - `curve_3d`: exact `CurveGeometry::Arc`.
//! - both pcurves: `Pcurve::Linear` (linear u with frame-phase
//!   offset, constant v at the cone's intersection height).
//! - `ParameterMap::Identity` via shared θ.

use crate::geometry::curve::CurveGeometry;
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

#[allow(clippy::too_many_arguments)]
pub(super) fn try_cone_cone_coaxial(
    apex_a: Point3,
    axis_a: Vec3,
    half_angle_a: f64,
    apex_b: Point3,
    axis_b: Vec3,
    half_angle_b: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    // Axes must be parallel (allow opposite signs).
    if axis_a.dot(&axis_b).abs() < 1.0 - tol.angular() {
        return None;
    }
    // Apex_b must lie on cone_a's axis.
    let dp = apex_b - apex_a;
    let along = dp.dot(&axis_a);
    let perp = dp - axis_a * along;
    if perp.length() > tol.linear() {
        return None;
    }
    // Use cone_a's axis as the canonical axis direction.
    let axis = axis_a;
    let axis_b_sign = axis_b.dot(&axis_a).signum();
    let _ = axis_b_sign;

    let t_a = half_angle_a.tan();
    let t_b = half_angle_b.tan();
    if (t_a - t_b).abs() < tol.solver() {
        return None; // matching half-angles ⇒ coincident or parallel
    }

    // Solve v_a from v_a · t_a = v_b · t_b and v_a − v_b = d.
    let d = along; // signed apex offset, apex_a → apex_b along cone_a axis
    let v_a = d * t_b / (t_b - t_a);
    let v_b = v_a - d;
    if v_a <= tol.linear() || v_b <= tol.linear() {
        return None; // intersection on negative-v half (other nappe)
    }
    let circle_radius = v_a * t_a;
    if circle_radius <= tol.linear() {
        return None;
    }
    let circle_center = apex_a + axis * v_a;

    let (perp_a_a, perp_a_b) = crate::geometry::analytic_surfaces::orthonormal_basis(axis_a);
    let (perp_b_a, perp_b_b) = crate::geometry::analytic_surfaces::orthonormal_basis(axis_b);
    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(axis);
    let cone_a_u_offset = arc_basis_u.dot(&perp_a_b).atan2(arc_basis_u.dot(&perp_a_a));
    let cone_b_u_offset = arc_basis_u.dot(&perp_b_b).atan2(arc_basis_u.dot(&perp_b_a));

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
        normal: axis,
        radius: circle_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    let pcurve_cone_a = Pcurve::linear(
        (theta_start + cone_a_u_offset, v_a),
        (theta_end + cone_a_u_offset, v_a),
        theta_start,
        theta_end,
    );
    let pcurve_cone_b = Pcurve::linear(
        (theta_start + cone_b_u_offset, v_b),
        (theta_end + cone_b_u_offset, v_b),
        theta_start,
        theta_end,
    );

    let cone_a_surf = SurfaceGeometry::Cone {
        apex: apex_a,
        axis: axis_a,
        half_angle: half_angle_a,
    };
    let cone_b_surf = SurfaceGeometry::Cone {
        apex: apex_b,
        axis: axis_b,
        half_angle: half_angle_b,
    };

    let triple = EdgeTriple {
        curve_3d,
        pcurve_a: pcurve_cone_a,
        pcurve_b: pcurve_cone_b,
        t_start: theta_start,
        t_end: theta_end,
        topology: EdgeTopology::Open,
        parameter_map: ParameterMap::Identity,
        provenance: EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::ConeCone,
            branch: BranchSelector::OnlyBranch,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple,
        &cone_a_surf,
        &cone_b_surf,
        tol,
    ))
}
