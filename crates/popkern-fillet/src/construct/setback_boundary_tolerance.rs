//! Entity-tolerance recording for the setback fan's accepted
//! boundary inexactness. Split from `setback_fan_seam_refresh.rs`
//! (500-line structural cap); the contract rationale lives on the
//! function.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Record the fan's accepted boundary inexactness as per-edge
/// tolerances (the entity-tolerance contract check 9 already honors:
/// "the edge's own tolerance when set, else the model tolerance").
///
/// The setback cap's rim rows are FITTED onto the live sweeps under
/// the cap-outer G0 budget (5e-3); the residual that fit legitimately
/// leaves (measured 2.1e-5..2.8e-3 on the blunt fixtures) is real,
/// bounded, and accepted — but with nothing recorded, the serving
/// curve-surface validation fails the body at the default band and a
/// Part-level caller cannot receive the fan at all. Each boundary
/// edge of the cap faces gets `edge.tolerance` stamped with the
/// measured worst curve-vs-adjacent-surface distance (1.25x
/// headroom, global nearest-point projection). Never lowers an
/// existing tolerance; sub-1e-6 residuals are not recorded. This
/// pass writes TOLERANCE FIELDS ONLY — it must not feed the
/// tessellator's fallback decisions (an earlier variant that widened
/// the mesher's off-surface gate with these values repartitioned the
/// boundary-first/uniform-grid face split and cracked previously
/// watertight fixtures: hexagon6_blunt 0 -> 482 open edges —
/// falsified, reverted).
pub(super) fn record_boundary_tolerances(brep: &mut BRep, cap_faces: &[FaceHandle]) {
    use popkern_core::geometry::surface::SurfaceEval;
    let mut edges: Vec<EdgeHandle> = Vec::new();
    for &face in cap_faces {
        let Some(face_data) = brep.faces.get(face) else {
            continue;
        };
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut he = start;
        for _ in 0..64 {
            let Some(hd) = brep.half_edges.get(he) else {
                break;
            };
            if !edges.contains(&hd.edge) {
                edges.push(hd.edge);
            }
            he = hd.next;
            if he == start {
                break;
            }
        }
    }
    // One vertex-ring out: the rim refit's endpoint-coherence repair
    // pulls the incident arc/corner-path curves' ends onto the moved
    // fan-corner vertices, leaving an accepted ~1e-3 deviation on
    // edges OUTSIDE the cap loops (measured: the blunt tower's
    // serving gate declined at 1.0040e-3 on such an edge while the
    // fan itself was watertight). Record on those too.
    let ring_vertices: std::collections::BTreeSet<popkern_core::topology::VertexHandle> = edges
        .iter()
        .filter_map(|&eh| {
            let edge = brep.edges.get(eh)?;
            let he = brep.half_edges.get(edge.half_edge)?;
            let twin = brep.half_edges.get(he.twin)?;
            Some([he.origin, twin.origin])
        })
        .flatten()
        .collect();
    for (eh, edge) in brep.edges.iter() {
        if edges.contains(&eh) {
            continue;
        }
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        if ring_vertices.contains(&he.origin) || ring_vertices.contains(&twin.origin) {
            edges.push(eh);
        }
    }
    for eh in edges {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let ce: &dyn CurveEval = &curve.geometry;
        let samples: Vec<Point3> = (0..=16)
            .map(|i| ce.point_at(curve.t_start + i as f64 / 16.0 * (curve.t_end - curve.t_start)))
            .collect();
        let mut residual = 0.0_f64;
        let Some(he_data) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        for half in [edge.half_edge, he_data.twin] {
            let Some(face) = brep
                .half_edges
                .get(half)
                .and_then(|hd| brep.loops.get(hd.loop_))
                .map(|l| l.face)
            else {
                continue;
            };
            let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            for p in &samples {
                let Some((u, v)) = popkern_core::geometry::uv_projection::project_point_to_uv(
                    &surface.geometry,
                    p,
                ) else {
                    continue;
                };
                residual = residual.max(surface.geometry.point_at(u, v).distance_to(p));
            }
        }
        if residual > 1.0e-6 {
            let recorded = residual * 1.25;
            let edge = &mut brep.edges[eh];
            if edge.tolerance.map(|t| recorded > t).unwrap_or(true) {
                edge.tolerance = Some(recorded);
            }
            // Entity tolerances must stay inside the body's
            // imprecision budget or validate() fails the commit with
            // ExceedsMaxBound (measured: the whole surgery rolled
            // back on octagon8_mixed and the G2 vertex-blend gate,
            // shipping fallbacks). Widening the budget for recorded
            // inexactness is the sanctioned pattern (see
            // brep_edge_surface_compatibility_helpers).
            if recorded > brep.max_imprecision {
                brep.max_imprecision = recorded * 1.01;
            }
        }
    }
}
