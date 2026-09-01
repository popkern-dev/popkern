//! Vertex snap with adjacent-chord re-anchoring — split from
//! `edge_reconcile_fit.rs` (500-line ratchet).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle};

#[allow(dead_code)]
fn pt_zero() -> Point3 {
    Point3::new(0.0, 0.0, 0.0)
}

/// Move `vertex` to `target` and rebuild adjacent edges' endpoint
/// geometry: chord Lines are rebuilt between their (moved) vertex
/// positions; NURBS/polyline curves get their nearest endpoint
/// control point snapped when it sits within the same class of
/// distance. `except` (the edge being re-trimmed) is left alone.
pub(crate) fn snap_vertex_and_adjacent_chords(
    brep: &mut BRep,
    vertex: popkern_core::topology::VertexHandle,
    target: popkern_core::geometry::Point3,
    except: EdgeHandle,
) {
    let Some(v) = brep.vertices.get_mut(vertex) else {
        return;
    };
    let old_pos = v.position;
    if old_pos.distance_to(&target) < 1e-12 {
        return;
    }
    v.position = target;

    // Collect edges incident to the vertex.
    let incident: Vec<EdgeHandle> = brep
        .edges
        .iter()
        .filter(|(eh, e)| {
            *eh != except && {
                let he = e.half_edge;
                brep.half_edges.get(he).is_some_and(|h| {
                    h.origin == vertex
                        || brep
                            .half_edges
                            .get(h.twin)
                            .is_some_and(|t| t.origin == vertex)
                })
            }
        })
        .map(|(eh, _)| eh)
        .collect();
    for eh in incident {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let Some(curve) = edge.curve.clone() else {
            continue;
        };
        let he = edge.half_edge;
        let (va, vb) = {
            let h = &brep.half_edges[he];
            (h.origin, brep.half_edges[h.twin].origin)
        };
        let pa = brep.vertices[va].position;
        let pb = brep.vertices[vb].position;
        match &curve.geometry {
            popkern_core::geometry::CurveGeometry::Line { .. } => {
                // Rebuild the line between the current vertex
                // positions.
                let dir = pb - pa;
                let len = dir.length();
                if len > 1e-12
                    && let Some(e) = brep.edges.get_mut(eh)
                {
                    e.curve = Some(popkern_core::geometry::CurveOnEdge {
                        geometry: popkern_core::geometry::CurveGeometry::Line {
                            origin: pa,
                            direction: dir * (1.0 / len),
                        },
                        t_start: 0.0,
                        t_end: len,
                    });
                }
            }
            popkern_core::geometry::CurveGeometry::Nurbs(n) => {
                // Snap the endpoint control point nearest the moved
                // vertex.
                let ceval: &dyn CurveEval = &curve.geometry;
                let start = ceval.point_at(curve.t_start);
                let end = ceval.point_at(curve.t_end);
                let mut n2 = n.clone();
                let snap_first = start.distance_to(&old_pos) <= end.distance_to(&old_pos);
                let moved = if snap_first {
                    n2.control_points.first_mut()
                } else {
                    n2.control_points.last_mut()
                };
                // The CP belongs to this endpoint when it sits much
                // closer to the old vertex than the move distance
                // (exact-match 1e-6 silently left junction-adjacent
                // curve starts 2e-2 behind a 5e-2 settle move —
                // measured on the boss smooth-junction settles).
                let move_len = old_pos.distance_to(&target);
                if let Some(cp) = moved
                    && cp.distance_to(&old_pos) < (move_len * 0.5).max(1e-6)
                {
                    *cp = target;
                    if let Some(e) = brep.edges.get_mut(eh) {
                        e.curve = Some(popkern_core::geometry::CurveOnEdge {
                            geometry: popkern_core::geometry::CurveGeometry::Nurbs(n2),
                            t_start: curve.t_start,
                            t_end: curve.t_end,
                        });
                    }
                }
            }
            _ => {
                // Analytic curves (arcs etc.) keep their geometry;
                // a moved vertex within the snap class stays within
                // downstream association gates.
            }
        }
    }
}
