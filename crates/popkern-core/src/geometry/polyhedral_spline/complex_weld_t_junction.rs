//! Boundary T-junctions between surface complex components (PK-1119).
//!
//! Neighbouring components refined differently meet along a seam where
//! one side has a vertex the other's boundary edge runs straight past:
//! a coarse edge P→Q facing a finer chain Q→M→P. The weld pairs
//! boundary edges by their endpoints, so such a seam stays open and the
//! components touch only at P and Q, a sheet joined through a vertex
//! that re-reads as separate parts. Before pairing, each boundary edge
//! is split at every boundary vertex of another component lying on its
//! curve inside it, so both sides of the seam carry the same vertices
//! and every edge finds its partner.
//!
//! The test is geometric and local: the vertex within the weld's own C0
//! tolerance of the edge's curve and farther than it from both of the
//! edge's ends. A vertex beside an edge but off its curve is no
//! junction and is left alone; nothing moves, the edge is split at the
//! curve point nearest the vertex.

use std::collections::HashMap;

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::curve_inversion::invert_curve_at_point;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, edge_segment};

/// Split `brep`'s boundary edges at the T-junctions where a boundary
/// vertex of another component lies on them, within `tolerance`; each
/// new edge joins its parent's component in `edge_component`. Returns
/// how many splits were made.
///
/// # Errors
///
/// The message of a split [`BRep::split_edge_at_point`] declines.
pub(super) fn split_boundary_t_junctions(
    brep: &mut BRep,
    edge_component: &mut HashMap<EdgeHandle, usize>,
    tolerance: f64,
) -> Result<usize, String> {
    let tol = Tolerance::default().with_linear(tolerance);
    let mut splits = 0;
    // Each split adds one boundary vertex and resolves one junction:
    // the boundary vertex count bounds the passes.
    let limit = 2 * brep.free_edges().len() + 1;
    while splits < limit {
        let Some((edge, point)) = next_junction(brep, edge_component, &tol) else {
            break;
        };
        let component = edge_component.get(&edge).copied().unwrap_or_default();
        let split = brep
            .split_edge_at_point(edge, point, &tol)
            .map_err(|e| format!("boundary T-junction split declined: {e}"))?;
        edge_component.insert(split.edge_a, component);
        edge_component.insert(split.edge_b, component);
        splits += 1;
    }
    Ok(splits)
}

/// The first boundary edge with another component's boundary vertex
/// inside it on its curve, and that vertex's position.
fn next_junction(
    brep: &BRep,
    edge_component: &HashMap<EdgeHandle, usize>,
    tol: &Tolerance,
) -> Option<(EdgeHandle, Point3)> {
    let band = tol.linear();
    let boundary: Vec<(EdgeHandle, usize, Point3, Point3)> = brep
        .free_edges()
        .into_iter()
        .filter_map(|e| {
            let (start, end) = edge_segment(brep, e)?;
            Some((e, edge_component.get(&e).copied()?, start, end))
        })
        .collect();
    for &(edge, component, start, end) in &boundary {
        let Some(curve) = brep.edges.get(edge).and_then(|e| e.curve.as_ref()) else {
            continue;
        };
        let lo = curve.t_start.min(curve.t_end);
        let hi = curve.t_start.max(curve.t_end);
        let ends = boundary
            .iter()
            .filter(|b| b.1 != component)
            .flat_map(|b| [b.2, b.3]);
        for p in ends {
            if p.distance_to(&start) <= band || p.distance_to(&end) <= band {
                continue;
            }
            let Ok(t) = invert_curve_at_point(curve, p, tol) else {
                continue;
            };
            if t <= lo || t >= hi || curve.geometry.point_at(t).distance_to(&p) > band {
                continue;
            }
            return Some((edge, p));
        }
    }
    None
}
