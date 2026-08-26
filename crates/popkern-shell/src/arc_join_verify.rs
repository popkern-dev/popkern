//! Arc-join verification and rim-edge curve repair: the geometric
//! postconditions every join face must pass before serving, and the
//! offset-rail rebuild for free-rim edges whose Stage-2 recompute
//! declined (PK-326 boundary corners).

use std::collections::HashMap;

use crate::arc_join_geometry::offset_rail;
use crate::error::ShellError;
use crate::skin_map::{JoinRecord, SkinMap};
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Rebuild a FREE-RIM edge's curve as the offset image of its current
/// (stale clone) curve on the live side's offset surface. Returns
/// `false` when the edge is not a rim edge or the rebuild fails.
pub(crate) fn rebuild_rim_edge_curve(
    skin: &mut BRep,
    eh: EdgeHandle,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    skin_map: &SkinMap,
    tolerance: &Tolerance,
) -> bool {
    let Some(edge) = skin.edges.get(eh) else {
        return false;
    };
    let Some(curve) = edge.curve.clone() else {
        return false;
    };
    let Some(he) = skin.half_edges.get(edge.half_edge) else {
        return false;
    };
    let Some(twin) = skin.half_edges.get(he.twin) else {
        return false;
    };
    // Exactly one live side qualifies this as a rim edge.
    let live = match (
        skin.loops.get(he.loop_).map(|l| l.face),
        skin.loops.get(twin.loop_).map(|l| l.face),
    ) {
        (Some(f), None) | (None, Some(f)) => f,
        _ => return false,
    };
    let reverse_face = skin_map.reverse_face_map();
    let Some(outer_face) = reverse_face.get(&live) else {
        return false;
    };
    let Some(offset) = offset_surfaces.get(outer_face) else {
        return false;
    };
    let Some((rail, _dev)) = offset_rail(&curve, offset, tolerance) else {
        return false;
    };
    if let Some(edge) = skin.edges.get_mut(eh) {
        edge.curve = Some(rail);
        true
    } else {
        false
    }
}

/// Geometric postconditions for every join face: boundary curves exist,
/// end at their vertices, and lie on the join surface.
pub(crate) fn verify_join_faces(
    skin: &BRep,
    records: &[JoinRecord],
    rail_deviation: f64,
    tolerance: &Tolerance,
) -> Result<(), ShellError> {
    let gate = (tolerance.linear() * 100.0).max(2.0 * rail_deviation + tolerance.linear());
    for record in records {
        let face = record.face;
        let Some(surface) = skin.faces.get(face).and_then(|f| f.surface.as_ref()) else {
            return Err(join_verify_error(face, "join face has no surface"));
        };
        let loop_ = skin.faces[face].outer_loop;
        let start = skin.loops[loop_].half_edge;
        let mut cursor = start;
        let mut steps = 0usize;
        loop {
            steps += 1;
            if steps > skin.half_edges.len() {
                return Err(join_verify_error(
                    face,
                    "join face boundary walk does not close (broken next cycle)",
                ));
            }
            let he = &skin.half_edges[cursor];
            let twin = &skin.half_edges[he.twin];
            let edge = &skin.edges[he.edge];
            let Some(curve) = edge.curve.as_ref() else {
                return Err(join_verify_error(face, "join boundary edge has no curve"));
            };
            let a = skin.vertices[he.origin].position;
            let b = skin.vertices[twin.origin].position;
            let c0 = curve.geometry.point_at(curve.t_start);
            let c1 = curve.geometry.point_at(curve.t_end);
            let fwd = a.distance_to(&c0).max(b.distance_to(&c1));
            let rev = a.distance_to(&c1).max(b.distance_to(&c0));
            if fwd.min(rev) > gate {
                return Err(join_verify_error(
                    face,
                    &format!(
                        "join boundary curve misses its vertices by {:.3e} (gate {gate:.3e})",
                        fwd.min(rev)
                    ),
                ));
            }
            for frac in [0.25, 0.5, 0.75] {
                let t = curve.t_start + frac * (curve.t_end - curve.t_start);
                let p = curve.geometry.point_at(t);
                let Some(d) = crate::curved_recompute::surface_distance(&surface.geometry, &p)
                else {
                    return Err(join_verify_error(
                        face,
                        "join surface distance query failed",
                    ));
                };
                if d > gate {
                    return Err(join_verify_error(
                        face,
                        &format!(
                            "join boundary curve leaves the join surface by {d:.3e} \
                             (gate {gate:.3e})"
                        ),
                    ));
                }
            }
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
    }
    Ok(())
}

pub(crate) fn join_verify_error(face: FaceHandle, mechanism: &str) -> ShellError {
    ShellError::Topology(format!(
        "arc-join geometric postcondition failed on join face {face:?}: {mechanism}"
    ))
}
