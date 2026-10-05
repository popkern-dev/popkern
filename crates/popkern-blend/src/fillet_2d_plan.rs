//! Pre-surgery analysis for `fillet_corner_2d`.
//!
//! This module owns the curve-agnostic geometry that computes
//! where the fillet's tangent points and arc center land on the
//! two adjacent edges of a sheet-body corner. Splitting it out of
//! `fillet_2d.rs` keeps each file under the 500-line structural
//! cap and reflects the natural boundary between "analysis"
//! (compute a `CornerPlan`) and "surgery" (consume it via
//! `apply_corner_plan`).
//!
//! ## Dispatch
//!
//! [`corner_plan`] dispatches on the curve types of the two
//! adjacent edges:
//!
//! - `(Line, Line)` — classic 2D fillet math: `setback = r /
//!   tan(θ/2)`, arc centre offset `= r / sin(θ/2)` along the
//!   interior bisector.
//! - `(Line, Arc)` / `(Arc, Line)` — offset-curve intersection
//!   via [`solve_line_arc_fillet`]. The fillet centre sits on an
//!   offset-parallel of the line (distance `r` on the interior
//!   side) and on a concentric offset of the arc (radius `R ± r`
//!   depending on internal vs external tangency), which yields a
//!   scalar quadratic in the line parameter.
//! - `(Arc, Arc)` and other combinations — rejected (future
//!   work).

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};
use popkern_math::PortableMath;

use crate::corner_2d_helpers::{
    circle_basis_like_core, edge_arc_length, setback_from_corner, tangent_at_corner_outgoing,
};
use crate::fillet_2d_solvers::{solve_arc_arc_fillet, solve_generic_fillet, solve_line_arc_fillet};

/// Pre-surgery analysis captured before we touch the mutable BRep.
/// Collecting everything up front lets us bail out cleanly on any
/// input validation failure without having to roll back partial
/// state.
#[derive(Debug, Clone, Copy)]
pub(crate) struct CornerPlan {
    /// Outgoing half-edge at the corner whose edge will carry the
    /// first tangent point `t1`. Walked out of the corner vertex
    /// toward the far end of `edge_a`.
    pub(crate) he_out_a: HalfEdgeHandle,
    /// Outgoing half-edge at the corner toward the far end of
    /// `edge_b`.
    pub(crate) he_out_b: HalfEdgeHandle,
    /// Tangent point on `edge_a`.
    pub(crate) t1: Point3,
    /// Tangent point on `edge_b`.
    pub(crate) t2: Point3,
    /// Arc centre for the fillet's new edge.
    pub(crate) arc_center: Point3,
    /// Unit arc-plane normal — `dir_a × dir_b` normalised, where
    /// `dir_a`, `dir_b` are the outgoing unit tangents at the corner.
    pub(crate) arc_normal: Vec3,
    /// Arc start angle at `t1` in the `circle_basis(arc_normal)` frame.
    pub(crate) start_angle: f64,
    /// Arc end angle at `t2` in the same frame.
    pub(crate) end_angle: f64,
    /// Fillet radius (copied in for convenience).
    pub(crate) radius: f64,
}

pub(crate) fn corner_plan(
    brep: &BRep,
    vertex: VertexHandle,
    radius: f64,
) -> Result<CornerPlan, String> {
    // Collect all outgoing half-edges at the corner. A sheet-body
    // corner has two: one on the front face's outer loop, one on
    // the backside face's loop. Each points to a different edge.
    let mut outgoing: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.origin == vertex)
        .map(|(h, _)| h)
        .collect();
    outgoing.sort();

    if outgoing.len() != 2 {
        return Err(format!(
            "corner vertex must have exactly 2 outgoing half-edges, got {}",
            outgoing.len()
        ));
    }

    let corner_pos = brep.vertices[vertex].position;

    let edge_a = brep.half_edges[outgoing[0]].edge;
    let edge_b = brep.half_edges[outgoing[1]].edge;
    if edge_a == edge_b {
        return Err("both outgoing half-edges belong to the same edge (degenerate loop)".into());
    }

    // Unit tangent directions at the corner, outward along each
    // adjacent edge. `tangent_at_corner_outgoing` handles both
    // `Line` and `Arc` curves; other curve kinds are rejected here.
    let dir_a = tangent_at_corner_outgoing(brep, outgoing[0])?;
    let dir_b = tangent_at_corner_outgoing(brep, outgoing[1])?;

    let len_a = edge_arc_length(brep, edge_a)?;
    let len_b = edge_arc_length(brep, edge_b)?;
    if len_a < 1e-12 || len_b < 1e-12 {
        return Err("adjacent edge has near-zero length".into());
    }

    let cos_angle = dir_a.dot(&dir_b).clamp(-1.0, 1.0);
    let angle = cos_angle.pacos();
    if angle < 1.0e-4 {
        return Err("edges are nearly colinear (no meaningful corner)".into());
    }
    if (std::f64::consts::PI - angle).abs() < 1.0e-4 {
        return Err("edges form a near-180° corner (degenerate)".into());
    }

    // Interior bisector and fillet plane normal are curve-type-
    // agnostic — they only depend on the tangent directions at the
    // corner, which `tangent_at_corner_outgoing` already resolved.
    let bisector = (dir_a + dir_b)
        .normalized()
        .ok_or_else(|| "bisector degenerates — edges near-parallel".to_string())?;
    let raw_normal = dir_a.cross(&dir_b);
    let arc_normal = raw_normal
        .normalized()
        .ok_or_else(|| "arc normal degenerate".to_string())?;

    // Dispatch on curve types to compute the tangent points and
    // fillet center. (Line, Line) uses the classic straight-line
    // formula; (Line, Arc) / (Arc, Line) solves offset-line-
    // offset-circle intersection. `tangent_at_corner_outgoing`
    // above has already rejected edges without curve geometry, so
    // this re-unwrap pathway is defensive.
    let curve_a = brep.edges[edge_a]
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge_a:?} has no curve geometry"))?;
    let curve_b = brep.edges[edge_b]
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge_b:?} has no curve geometry"))?;

    // If either adjacent curve is a NURBS, dispatch to the
    // generic Newton solver — offset curves of NURBS aren't
    // NURBS in general, so closed-form intersection doesn't
    // apply. Line-Line, Line-Arc, Arc-Line, and Arc-Arc keep
    // their specialised closed-form paths below.
    let involves_nurbs = matches!(curve_a.geometry, CurveGeometry::Nurbs(_))
        || matches!(curve_b.geometry, CurveGeometry::Nurbs(_));

    let (t1, t2, arc_center) = if involves_nurbs {
        solve_generic_fillet(
            &curve_a.geometry,
            &curve_b.geometry,
            corner_pos,
            dir_a,
            dir_b,
            bisector,
            arc_normal,
            radius,
            angle,
        )?
    } else {
        match (&curve_a.geometry, &curve_b.geometry) {
            (CurveGeometry::Line { .. }, CurveGeometry::Line { .. }) => {
                let half_angle = angle * 0.5;
                let tan_half = half_angle.ptan();
                let sin_half = half_angle.psin();
                let setback = radius / tan_half;
                let offset = radius / sin_half;
                if setback >= len_a || setback >= len_b {
                    return Err(format!(
                        "radius {radius:.6} too large: setback {setback:.6} ≥ adjacent edge length ({len_a:.6} / {len_b:.6})"
                    ));
                }
                let t1 = corner_pos + dir_a * setback;
                let t2 = corner_pos + dir_b * setback;
                let arc_center = corner_pos + bisector * offset;
                (t1, t2, arc_center)
            }
            (
                CurveGeometry::Line { .. },
                CurveGeometry::Arc {
                    center: ac,
                    radius: ar,
                    ..
                },
            ) => {
                let (line_tp, arc_tp, fc) =
                    solve_line_arc_fillet(corner_pos, dir_a, *ac, *ar, bisector, radius)?;
                (line_tp, arc_tp, fc)
            }
            (
                CurveGeometry::Arc {
                    center: ac,
                    radius: ar,
                    ..
                },
                CurveGeometry::Line { .. },
            ) => {
                let (line_tp, arc_tp, fc) =
                    solve_line_arc_fillet(corner_pos, dir_b, *ac, *ar, bisector, radius)?;
                (arc_tp, line_tp, fc)
            }
            (
                CurveGeometry::Arc {
                    center: ca,
                    radius: ra,
                    ..
                },
                CurveGeometry::Arc {
                    center: cb,
                    radius: rb,
                    ..
                },
            ) => {
                solve_arc_arc_fillet(corner_pos, *ca, *ra, *cb, *rb, bisector, arc_normal, radius)?
            }
            (other_a, other_b) => {
                return Err(format!(
                    "unsupported curve types for 2D fillet: {other_a:?} / {other_b:?}"
                ));
            }
        }
    };

    // Universal setback checks: the computed tangent points must
    // lie ahead of the corner (positive setback) and within the
    // respective edge's arc-length budget. For the Line-Line branch
    // above this is a redundancy check; for Line-Arc / Arc-Line it's
    // the primary bounds check on the arc side.
    let setback_a = setback_from_corner(brep, outgoing[0], t1)?;
    let setback_b = setback_from_corner(brep, outgoing[1], t2)?;
    if setback_a < -1e-9 || setback_b < -1e-9 {
        return Err(format!(
            "fillet tangent point lands behind corner (setbacks {setback_a:.6}, {setback_b:.6})"
        ));
    }
    if setback_a >= len_a || setback_b >= len_b {
        return Err(format!(
            "radius {radius:.6} too large: setbacks {setback_a:.6} / {setback_b:.6} ≥ adjacent edge lengths ({len_a:.6} / {len_b:.6})"
        ));
    }

    // `CurveGeometry::Arc::point_at` internally uses the
    // `circle_basis(normal)` frame, which picks `u` by crossing
    // the normal with +x or +y (whichever isn't parallel) and
    // derives `v = normal × u`. We must compute our start and end
    // angles IN THAT FRAME so the evaluated arc endpoints actually
    // land on `t1` and `t2`. Replicating `circle_basis` locally
    // here avoids coupling to a private helper.
    let (u_basis, v_basis) = circle_basis_like_core(arc_normal);
    let t1_rel = t1 - arc_center;
    let t2_rel = t2 - arc_center;
    let start_angle = t1_rel.dot(&v_basis).patan2(t1_rel.dot(&u_basis));
    let end_angle = t2_rel.dot(&v_basis).patan2(t2_rel.dot(&u_basis));

    Ok(CornerPlan {
        he_out_a: outgoing[0],
        he_out_b: outgoing[1],
        t1,
        t2,
        arc_center,
        arc_normal,
        start_angle,
        end_angle,
        radius,
    })
}
