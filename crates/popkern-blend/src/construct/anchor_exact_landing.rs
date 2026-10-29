//! Exact anchor landing for chamfer bevels (PK-994, PK-1008): a plane
//! intersection lands on existing topology only within the linear
//! band. Split from `helpers_anchor_insert.rs` (line budget).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

/// Land an EXACT anchor (a chamfer bevel's plane intersection) on
/// existing topology: reuse a vertex, or split a straight edge, only
/// when the position lies on it within the linear band. The loose
/// gates of [`insert_anchor_at_position`] scale with the curve sagitta
/// they were built for (`linear × 1e4`, an absolute 1e-2 at the default
/// band) and captured a pocket rim's reflex rail end on the block's
/// OUTER rim 4e-3 away at 1/100 scale (PK-994).
pub(super) fn land_on_boundary_exactly(
    brep: &mut BRep,
    position: Point3,
    tolerance: &Tolerance,
    edge_splits: &mut Vec<(EdgeHandle, EdgeHandle, EdgeHandle)>,
) -> Result<VertexHandle, String> {
    let gate = (tolerance.linear() * 10.0).max(1e-12);
    if let Some((vh, _)) = brep
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&position) <= gate)
    {
        return Ok(vh);
    }
    let mut hit: Option<(EdgeHandle, f64, Point3)> = None;
    for (eh, edge) in brep.edges.iter() {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let (a, b) = (
            brep.vertices[he.origin].position,
            brep.vertices[twin.origin].position,
        );
        let ab = b - a;
        let len_sq = ab.dot(&ab);
        if len_sq <= gate * gate {
            continue;
        }
        let t = (position - a).dot(&ab) / len_sq;
        let foot = a + ab * t;
        // Interior of a straight edge: the chord is the edge.
        let straight = edge.curve.as_ref().is_none_or(|c| {
            let m = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
            (m - a).cross(&ab).length() <= gate * len_sq.sqrt()
        });
        if t > 0.0 && t < 1.0 && straight && foot.distance_to(&position) <= gate {
            hit = Some((eh, t, foot));
            break;
        }
    }
    let (eh, t, foot) = hit.ok_or_else(|| "exact anchor is not on the boundary".to_string())?;
    // Split at the chord fraction (the edge is its chord): a sub-edge
    // can still carry its parent's curve span, and splitting by curve
    // inversion against that span placed a cone's rail end 1.7 away
    // from its foot. The new vertex must sit at the foot.
    let mut candidate = brep.clone();
    let split = candidate
        .split_edge(eh, t)
        .map_err(|e| format!("exact anchor split failed: {e}"))?;
    let landed = candidate
        .vertices
        .get(split.vertex)
        .map(|v| v.position.distance_to(&foot))
        .unwrap_or(f64::INFINITY);
    if landed > gate {
        return Err(format!(
            "exact anchor split landed {landed:.3e} from its foot"
        ));
    }
    *brep = candidate;
    edge_splits.push((eh, split.edge_a, split.edge_b));
    Ok(split.vertex)
}
