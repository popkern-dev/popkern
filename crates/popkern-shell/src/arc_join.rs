//! Arc-join construction: rolling-ball strips along gap edges, sphere
//! window caps at gap corners (design §5 Stage 3a, PK-325).
//!
//! All topology edits go through Euler operators (each transactional
//! with its own candidate-verify-commit): per gap edge two `mef`s build
//! a bigon strip on each side and one `kef` merges them into a
//! rail-bounded strip; per gap corner `split_vertex` peels each face
//! run to its offset corner (the connectors ARE the sphere arcs),
//! `mef`s close the detours into cap triangles and `kef`s dissolve the
//! interior connectors, leaving one k-gon sphere window. Closed gap
//! edges (cylinder rims, lens rims, bore rims) are `split_edge`-split
//! once and join the general path (their k=2 vertices need one split
//! and no cap).

use std::collections::{BTreeMap, BTreeSet, HashMap};

use crate::arc_join_geometry::{great_arc, offset_rail, strip_surface, surface_geo_normal};
use crate::curved_recompute::project_to_surface;
use crate::error::{JoinDecline, ShellError};
use crate::join_classification::{GapEdge, JoinPlan};
use crate::skin_map::{JoinKind, JoinRecord, SkinMap};
use crate::types::JoinStyle;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Construct all arc joins on one skin. Returns the join provenance
/// records; declines transactionally (the caller's input body is
/// untouched — the skin is a pipeline-internal clone).
///
/// # Errors
///
/// [`ShellError::JoinFailed`] when any gap edge or corner cannot be
/// served (spindle torus, antipodal arcs, projection ambiguity);
/// [`ShellError::Core`] on Euler-operator postcondition failures.
pub(crate) fn construct_arc_joins(
    skin: &mut BRep,
    plan: &JoinPlan,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    skin_map: &SkinMap,
    outer: &BRep,
    tolerance: &Tolerance,
) -> Result<Vec<JoinRecord>, ShellError> {
    if plan.effective != JoinStyle::Arc || plan.gap_edges.is_empty() {
        return Ok(Vec::new());
    }
    let reverse_vertex = skin_map.reverse_vertex_map();
    let reverse_face = skin_map.reverse_face_map();

    // Phase 0: split closed gap edges so every gap edge has two
    // distinct endpoints and every loop has at least two half-edges.
    let mut gap_edges: BTreeMap<EdgeHandle, GapEdge> = plan.gap_edges.clone();
    let closed: Vec<EdgeHandle> = gap_edges
        .keys()
        .copied()
        .filter(|eh| edge_is_closed(skin, *eh))
        .collect();
    for eh in closed {
        let gap = gap_edges[&eh].clone();
        let split = skin.split_edge(eh, 0.5)?;
        gap_edges.insert(split.edge_a, gap.clone());
        gap_edges.insert(split.edge_b, gap);
    }

    // Phase 1: strips.
    let mut declines: Vec<JoinDecline> = Vec::new();
    let mut records: Vec<JoinRecord> = Vec::new();
    let mut strip_sources: BTreeMap<FaceHandle, EdgeHandle> = BTreeMap::new();
    let mut corner_vertices: BTreeSet<VertexHandle> = BTreeSet::new();
    let mut rail_deviation = 0.0f64;

    for (eh, gap) in &gap_edges {
        let Some(clone_curve) = skin.edges.get(*eh).and_then(|e| e.curve.clone()) else {
            declines.push(decline(
                gap,
                "gap edge carries no curve to anchor the canal on",
            ));
            continue;
        };
        let off1 = &offset_surfaces[&gap.outer_face1];
        let off2 = &offset_surfaces[&gap.outer_face2];
        let Some((rail1, dev1)) = offset_rail(&clone_curve, off1, tolerance) else {
            declines.push(decline(
                gap,
                "rail projection onto the first offset surface failed",
            ));
            continue;
        };
        let Some((rail2, dev2)) = offset_rail(&clone_curve, off2, tolerance) else {
            declines.push(decline(
                gap,
                "rail projection onto the second offset surface failed",
            ));
            continue;
        };
        rail_deviation = rail_deviation.max(dev1).max(dev2);
        let surface = match strip_surface(&clone_curve, gap.ball_radius, off1, off2, &rail1, &rail2)
        {
            Ok(s) => s,
            Err(mechanism) => {
                declines.push(decline(gap, &mechanism));
                continue;
            }
        };

        // Surgery: bigon strip on each side, then merge across the
        // original edge.
        let h1 = skin.edges[*eh].half_edge;
        let h1_next = skin.half_edges[h1].next;
        let res1 = skin.mef(h1, h1_next)?;
        let h2 = skin.half_edges[h1].twin;
        let h2_next = skin.half_edges[h2].next;
        let res2 = skin.mef(h2, h2_next)?;
        skin.kef(h1)?;

        skin.edges[res1.edge].curve = Some(rail1);
        skin.edges[res2.edge].curve = Some(rail2);

        // The strip face survives as res1.face; measure its orientation
        // flag against the first face's flag-honest normal at the rail
        // contact (the pair is G1-tangent there).
        let mid_t = 0.5 * (clone_curve.t_start + clone_curve.t_end);
        let p_mid = clone_curve.geometry.point_at(mid_t);
        let orientation =
            strip_orientation(skin, res1.face, &surface, off1, &p_mid, gap).unwrap_or(true);
        skin.faces[res1.face].surface = Some(SurfaceOnFace {
            geometry: surface,
            orientation,
        });

        let v_a = skin.half_edges[skin.edges[res1.edge].half_edge].origin;
        let v_b = skin.half_edges[skin.half_edges[skin.edges[res1.edge].half_edge].twin].origin;
        corner_vertices.insert(v_a);
        corner_vertices.insert(v_b);
        strip_sources.insert(res1.face, gap.outer_edge);
        records.push(JoinRecord {
            kind: JoinKind::PipeFill {
                source_edge: gap.outer_edge,
            },
            face: res1.face,
        });
    }

    if !declines.is_empty() {
        return Err(ShellError::JoinFailed {
            style: JoinStyle::Arc,
            declines,
        });
    }

    // Phase 2: corners.
    let strip_faces: BTreeSet<FaceHandle> = strip_sources.keys().copied().collect();
    for v in corner_vertices {
        match crate::arc_join_corner::construct_corner(
            skin,
            v,
            &strip_faces,
            offset_surfaces,
            &reverse_vertex,
            &reverse_face,
            outer,
            tolerance,
        ) {
            Ok(Some(cap_face)) => {
                let source_vertex = reverse_vertex.get(&v).copied();
                records.push(JoinRecord {
                    kind: JoinKind::SphereFill { source_vertex },
                    face: cap_face,
                });
            }
            Ok(None) => {}
            Err(mechanism) => {
                return Err(ShellError::JoinFailed {
                    style: JoinStyle::Arc,
                    declines: vec![JoinDecline {
                        edge: reverse_vertex
                            .get(&v)
                            .and_then(|ov| first_edge_at(outer, *ov))
                            .unwrap_or_default(),
                        mechanism: format!("corner construction failed: {mechanism}"),
                        arc_alternative: false,
                    }],
                });
            }
        }
    }

    // Phase 3: geometric postconditions — every join-face boundary
    // curve must end at its vertices and lie on the join surface.
    verify_join_faces(skin, &records, rail_deviation, tolerance)?;

    Ok(records)
}

fn decline(gap: &GapEdge, mechanism: &str) -> JoinDecline {
    JoinDecline {
        edge: gap.outer_edge,
        mechanism: mechanism.to_string(),
        arc_alternative: false,
    }
}

fn edge_is_closed(skin: &BRep, eh: EdgeHandle) -> bool {
    let Some(edge) = skin.edges.get(eh) else {
        return false;
    };
    let Some(he) = skin.half_edges.get(edge.half_edge) else {
        return false;
    };
    let Some(twin) = skin.half_edges.get(he.twin) else {
        return false;
    };
    he.origin == twin.origin
}

fn first_edge_at(outer: &BRep, vh: VertexHandle) -> Option<EdgeHandle> {
    let anchor = outer.vertices.get(vh)?.outgoing;
    outer.half_edges.get(anchor).map(|he| he.edge)
}

/// Orientation flag for a strip: its flag-honest normal must agree with
/// the adjacent offset face's flag-honest normal at the (G1-tangent)
/// rail contact.
fn strip_orientation(
    skin: &BRep,
    strip_face: FaceHandle,
    strip_surface: &SurfaceGeometry,
    off1: &SurfaceGeometry,
    p_mid: &Point3,
    gap: &GapEdge,
) -> Option<bool> {
    let q1 = project_to_surface(off1, p_mid)?;
    let n_strip_geo = surface_geo_normal(strip_surface, &q1)?;
    let n_face_geo = surface_geo_normal(off1, &q1)?;
    // Find the skin face carrying off1's geometry (the neighbor across
    // rail1) to read its orientation flag: it is the face of the loop
    // holding the twin of the strip's rail1 half-edge.
    let flag = neighbor_flag(skin, strip_face, gap)?;
    let n_face_flag = if flag { n_face_geo } else { n_face_geo * -1.0 };
    Some(n_strip_geo.dot(&n_face_flag) > 0.0)
}

/// The orientation flag of the offset face on side 1 of the strip: walk
/// the strip's loop and read the first neighbor whose face is not the
/// strip itself.
fn neighbor_flag(skin: &BRep, strip_face: FaceHandle, _gap: &GapEdge) -> Option<bool> {
    let loop_ = skin.faces.get(strip_face)?.outer_loop;
    let start = skin.loops.get(loop_)?.half_edge;
    let mut cursor = start;
    for _ in 0..skin.half_edges.len() {
        let he = skin.half_edges.get(cursor)?;
        let twin = skin.half_edges.get(he.twin)?;
        let nb_face = skin.loops.get(twin.loop_)?.face;
        if nb_face != strip_face {
            return skin
                .faces
                .get(nb_face)?
                .surface
                .as_ref()
                .map(|s| s.orientation);
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    None
}

/// Geometric postconditions for every join face: boundary curves exist,
/// end at their vertices, and lie on the join surface.
fn verify_join_faces(
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
        loop {
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

fn join_verify_error(face: FaceHandle, mechanism: &str) -> ShellError {
    ShellError::Topology(format!(
        "arc-join geometric postcondition failed on join face {face:?}: {mechanism}"
    ))
}

/// Shared helper for corner construction: the offset corner position of
/// an original point on a face's offset surface.
pub(crate) fn offset_corner_position(
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    outer_face: FaceHandle,
    p_orig: &Point3,
) -> Option<(Point3, f64)> {
    let off = offset_surfaces.get(&outer_face)?;
    let q = project_to_surface(off, p_orig)?;
    let r = q.distance_to(p_orig);
    Some((q, r))
}

/// Author the great-circle arc curve for a cap-boundary edge from its
/// current vertex positions.
pub(crate) fn author_corner_arc(
    skin: &mut BRep,
    edge: EdgeHandle,
    center: Point3,
    radius: f64,
    tolerance: &Tolerance,
) -> Result<(), String> {
    let he = skin.edges[edge].half_edge;
    let a = skin.vertices[skin.half_edges[he].origin].position;
    let twin = skin.half_edges[he].twin;
    let b = skin.vertices[skin.half_edges[twin].origin].position;
    let arc: CurveOnEdge = great_arc(center, radius, a, b, tolerance)
        .ok_or_else(|| "corner arc endpoints are antipodal or degenerate".to_string())?;
    skin.edges[edge].curve = Some(arc);
    Ok(())
}

/// Measured orientation for a sphere cap: agree with the adjacent strip
/// (tangent across the shared arc).
pub(crate) fn cap_orientation(
    skin: &BRep,
    cap_face: FaceHandle,
    sphere: &SurfaceGeometry,
) -> Option<bool> {
    let loop_ = skin.faces.get(cap_face)?.outer_loop;
    let start = skin.loops.get(loop_)?.half_edge;
    let mut cursor = start;
    for _ in 0..skin.half_edges.len() {
        let he = skin.half_edges.get(cursor)?;
        let twin = skin.half_edges.get(he.twin)?;
        let nb_face = skin.loops.get(twin.loop_)?.face;
        if nb_face != cap_face
            && let Some(nb) = skin.faces.get(nb_face).and_then(|f| f.surface.as_ref())
            && let Some(curve) = skin.edges.get(he.edge).and_then(|e| e.curve.as_ref())
        {
            let mid = curve.geometry.point_at(0.5 * (curve.t_start + curve.t_end));
            let n_cap = surface_geo_normal(sphere, &mid)?;
            let n_nb_geo = surface_geo_normal(&nb.geometry, &mid)?;
            let n_nb = if nb.orientation {
                n_nb_geo
            } else {
                n_nb_geo * -1.0
            };
            return Some(n_cap.dot(&n_nb) > 0.0);
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    None
}

#[cfg(test)]
#[path = "arc_join_tests.rs"]
mod tests;
