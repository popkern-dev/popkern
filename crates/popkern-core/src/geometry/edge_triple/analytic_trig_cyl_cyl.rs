//! Cylinder × cylinder analytic dispatch.
//!
//! When the two cylinder axes are parallel, the cross-sections are
//! two coplanar circles. Their intersection is 0, 1, or 2 points
//! that, lifted along the common axis, form 0, 1, or 2 lines in 3-D.
//!
//! Each line is parallel to the common axis, with `u` constant on
//! each cylinder (the angular position of the intersection point).
//! The pcurve on each cyl is a vertical line in `(u, v)` space.
//!
//! Skew or intersecting-axes cases are out of scope — those produce
//! ellipses (intersecting) or generic 4th-order curves (skew), and
//! fall through to the marching solver.

use crate::geometry::curve::CurveGeometry;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, EdgeHint, SurfaceSide};
use crate::geometry::edge_triple::types::{
    AnalyticPair, BranchSelector, EdgeTopology, EdgeTriple, EdgeTripleProvenance,
    EdgeTripleQuality, ParameterMap,
};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

#[allow(clippy::too_many_arguments)]
pub(super) fn try_cyl_cyl_parallel(
    cyl_a_origin: Point3,
    cyl_a_axis: Vec3,
    cyl_a_radius: f64,
    cyl_b_origin: Point3,
    cyl_b_axis: Vec3,
    cyl_b_radius: f64,
    anchors: &EdgeAnchors,
    hint: Option<&EdgeHint>,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    // Detect parallel axes.
    let dot = cyl_a_axis.dot(&cyl_b_axis);
    if dot.abs() < 1.0 - tol.angular() {
        return None;
    }
    let same_direction = dot > 0.0;
    let common_axis = cyl_a_axis;

    // Project cyl_b origin onto the cross-section plane through
    // cyl_a origin (perpendicular to common_axis).
    let off = cyl_b_origin - cyl_a_origin;
    let v_offset = off.dot(&common_axis);
    let radial = off - common_axis * v_offset;
    let d = radial.length();
    if d < tol.linear() {
        return None; // co-linear axes — degenerate
    }
    if d >= cyl_a_radius + cyl_b_radius - tol.linear() {
        return None; // separated
    }
    if d <= (cyl_a_radius - cyl_b_radius).abs() + tol.linear() {
        return None; // one inside the other
    }

    // Two intersection points in the cross-section. Set up
    // 2-D coordinates in the cross-section plane: x along
    // `radial / d` (toward cyl_b), y perpendicular.
    let x_basis = radial * (1.0 / d);
    let y_basis = common_axis.cross(&x_basis);
    let a = (d * d + cyl_a_radius * cyl_a_radius - cyl_b_radius * cyl_b_radius) / (2.0 * d);
    let h = (cyl_a_radius * cyl_a_radius - a * a).max(0.0).sqrt();
    if h < tol.linear() {
        return None;
    }

    // Two lines, each parallel to common_axis.
    // Branch 0: y = +h. Branch 1: y = -h.
    let pick_top = match (anchors, hint) {
        (_, Some(EdgeHint::BranchIndex(i))) => *i == 0,
        (EdgeAnchors::Open { start, .. }, _) => (*start - cyl_a_origin).dot(&y_basis) >= 0.0,
        _ => true,
    };
    let y_sign = if pick_top { 1.0 } else { -1.0 };
    let line_origin_2d = (a, y_sign * h);
    let line_origin = cyl_a_origin + x_basis * line_origin_2d.0 + y_basis * line_origin_2d.1;

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };

    // Verify anchors lie on the line within tolerance.
    let anchor_check = |p: Point3| -> Option<f64> {
        let off = p - line_origin;
        let proj = off.dot(&common_axis);
        let perp = off - common_axis * proj;
        if perp.length() > tol.linear() {
            return None;
        }
        Some(proj)
    };
    let s_start = match anchor_check(start) {
        Some(s) => s,
        None => {
            return Some(Err(EdgeTripleError::AnchorOffSurface {
                which: SurfaceSide::A,
                distance: (start - line_origin).length(),
            }));
        }
    };
    let s_end = match anchor_check(end) {
        Some(s) => s,
        None => {
            return Some(Err(EdgeTripleError::AnchorOffSurface {
                which: SurfaceSide::B,
                distance: (end - line_origin).length(),
            }));
        }
    };
    let (s_lo, s_hi) = if s_start <= s_end {
        (s_start, s_end)
    } else {
        (s_end, s_start)
    };

    // Compute the intersection point's u-coordinate on each cyl
    // (the angular position around its axis).
    let (cyl_a_basis_u, cyl_a_basis_v) =
        crate::geometry::analytic_surfaces::orthonormal_basis(cyl_a_axis);
    let u_a = (line_origin - cyl_a_origin)
        .dot(&cyl_a_basis_v)
        .patan2((line_origin - cyl_a_origin).dot(&cyl_a_basis_u));
    // For cyl B, account for direction sign: if cyl_b_axis points
    // opposite to cyl_a_axis, the v-axis direction flips for cyl B.
    let cyl_b_axis_oriented = if same_direction {
        cyl_b_axis
    } else {
        cyl_b_axis * -1.0
    };
    let (cyl_b_basis_u, cyl_b_basis_v) =
        crate::geometry::analytic_surfaces::orthonormal_basis(cyl_b_axis_oriented);
    let off_b = line_origin - cyl_b_origin;
    let u_b = off_b.dot(&cyl_b_basis_v).patan2(off_b.dot(&cyl_b_basis_u));

    // Pcurve on each cyl: u = constant, v varies linearly along
    // the line. Since the line is parallel to common_axis, v on
    // cyl A is just `s` (the line parameter). On cyl B, account for
    // direction: if same, v = s + offset; if anti-parallel, v = -s + offset.
    let v_a_start = s_lo;
    let v_a_end = s_hi;
    let v_offset_b = (line_origin - cyl_b_origin).dot(&cyl_b_axis);
    let v_b_start = if same_direction {
        v_offset_b + s_lo
    } else {
        v_offset_b - s_lo
    };
    let v_b_end = if same_direction {
        v_offset_b + s_hi
    } else {
        v_offset_b - s_hi
    };

    let pcurve_a = Pcurve::linear((u_a, v_a_start), (u_a, v_a_end), s_lo, s_hi);
    let pcurve_b = Pcurve::linear((u_b, v_b_start), (u_b, v_b_end), s_lo, s_hi);

    let curve_3d = CurveGeometry::Line {
        origin: line_origin,
        direction: common_axis,
    };

    let cyl_a = SurfaceGeometry::Cylinder {
        origin: cyl_a_origin,
        axis: cyl_a_axis,
        radius: cyl_a_radius,
    };
    let cyl_b = SurfaceGeometry::Cylinder {
        origin: cyl_b_origin,
        axis: cyl_b_axis,
        radius: cyl_b_radius,
    };

    let triple = EdgeTriple {
        curve_3d,
        pcurve_a,
        pcurve_b,
        t_start: s_lo,
        t_end: s_hi,
        topology: EdgeTopology::Open,
        parameter_map: ParameterMap::Identity,
        provenance: EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::CylinderCylinder,
            branch: if pick_top {
                BranchSelector::Indexed(0)
            } else {
                BranchSelector::Indexed(1)
            },
            cardinality: None,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, &cyl_a, &cyl_b, tol,
    ))
}
