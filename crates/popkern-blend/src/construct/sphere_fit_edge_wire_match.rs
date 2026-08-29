//! Boundary-segment matching for sphere-fit edge wiring.

use super::sphere_fit_edge_wire_helpers::{collect_loop_half_edges, project_point_on_curve};
use super::*;

pub(super) type WireJob = (HalfEdgeHandle, usize, bool, Point3, Point3);

pub(super) fn build_wire_worklist(
    brep: &BRep,
    sub_faces: &[FaceHandle],
    semantic_cap_indices: &[Option<usize>],
    segment_curves: &[&NurbsCurve],
    segment_endpoints: &[[Point3; 2]],
    tol: f64,
    sphere_fit_wide_tol: bool,
) -> Vec<WireJob> {
    // Assignment scan: for each cap sub-face's outer-loop
    // half-edge, find the best matching segment (proximity, with
    // the wider tolerance). Build a worklist of
    // `(he, seg_idx, reversed, target_p_origin, target_p_dest)`.
    // We do NOT mutate during the scan so the matching sees
    // pre-relocation positions consistently.
    //
    // When injection is active, allow multiple
    // sub-edges to share a single segment (`sphere_fit_wide_tol`). The
    // cap face's outer loop after sliver-absorb topology surgery
    // can have 5+ sub-edges along a single logical boundary
    // (FilletEnd[k] / FaceEdge[k]); the original "consume each
    // segment once" policy left short sub-edges unwired, holding
    // worst-case `cap_outer_g0` at ~10 cm. Letting multiple
    // sub-edges share a segment as their 3D curve writes the
    // analytic boundary to ALL of them — endpoint mismatch
    // tolerated since downstream consumers (validation, pcurve
    // refit) operate on the curve geometry, not the vertex-to-
    // curve-endpoint identity.
    let allow_multi_sub_edge = sphere_fit_wide_tol;
    let trace_match = popkern_core::trace::trace_enabled("wire-match");
    // When injection is active, restrict
    // matching to each sub-face's semantic cap index. cap[j]
    // gets row 0 = FilletEnd[j] (segment 2j) and col 9 =
    // FaceEdge[j] (segment 2j+1). Twin face surface kind tells
    // us which boundary the edge sits on; combined with the cap
    // identity, the segment index is fully determined.
    //
    // Without this constraint, geometric-proximity matching
    // routes cap[1]'s cylinder edge to FilletEnd[2] (the wrong
    // fillet) when the rolling-ball footprint of cap[1] sits
    // closer in 3D to FilletEnd[2]'s analytic arc than to
    // FilletEnd[1]'s — empirically observed on pyramid4_mixed
    // because the rolling-ball z (0.755) and analytic-arc z
    // (0.818) gap (~6 cm) is comparable to the inter-fillet XY
    // separation, so all four FilletEnd arcs are roughly
    // equidistant. Rewritten exact-boundary loops can also rotate
    // vector order away from boundary order, so callers provide
    // cylinder/sphere identity when available.
    let constrained_by_cap_index = sphere_fit_wide_tol;
    let mut consumed = vec![false; segment_curves.len()];
    let mut worklist: Vec<WireJob> = Vec::new();
    // Cap-face set for cap-vs-cap edge
    // detection. The cap face's surface at wire time is NurbsSurface
    // (Gregory assignment happens AFTER wire in the pipeline), so
    // we can't filter by surface kind. Identify cap faces by
    // membership in `sub_faces` instead.
    let cap_face_set: std::collections::BTreeSet<FaceHandle> = sub_faces.iter().copied().collect();

    for (cap_k, &sub_face) in sub_faces.iter().enumerate() {
        let Some(face_data) = brep.faces.get(sub_face) else {
            continue;
        };
        let outer_loop = face_data.outer_loop;
        let Some(loop_data) = brep.loops.get(outer_loop) else {
            continue;
        };
        let start_he = loop_data.half_edge;

        let hes = collect_loop_half_edges(brep, start_he);

        for he_handle in hes {
            let Some(he) = brep.half_edges.get(he_handle) else {
                continue;
            };
            let v_origin = he.origin;
            let next_he = he.next;
            let Some(next_he_data) = brep.half_edges.get(next_he) else {
                continue;
            };
            let v_dest = next_he_data.origin;

            // Skip cap-vs-cap edges (twin
            // face is in `sub_faces`, i.e., another cap face).
            // They aren't FilletEnd or FaceEdge boundaries —
            // they're inner spokes / shared cap edges between
            // adjacent per-fillet caps. Wiring them would write
            // the wrong analytic curve and could corrupt cap-cap
            // continuity. Only wire cap-vs-non-cap edges.
            //
            // Earlier surface-based filter (Sphere | Gregory) was
            // ineffective because the cap face's surface at wire
            // time is still NurbsSurface (Gregory replaces it
            // *after* the wire step in the pipeline). Membership
            // in `sub_faces` is the topology-stable signal.
            let twin_he = he.twin;
            let twin_face = brep
                .half_edges
                .get(twin_he)
                .and_then(|twin| brep.loops.get(twin.loop_))
                .map(|l| l.face);
            let twin_is_cap = twin_face
                .map(|f| cap_face_set.contains(&f))
                .unwrap_or(false);
            if twin_is_cap {
                if trace_match {
                    eprintln!(
                        "[wire-match] he={he_handle:?} edge={:?} SKIPPED (twin is cap face)",
                        he.edge
                    );
                }
                continue;
            }

            let Some(origin_vertex) = brep.vertices.get(v_origin) else {
                continue;
            };
            let Some(dest_vertex) = brep.vertices.get(v_dest) else {
                continue;
            };
            let p_origin = origin_vertex.position;
            let p_dest = dest_vertex.position;

            // Matching strategy:
            //
            // * Default (`allow_multi_sub_edge = false`): endpoint-sum
            //   matching. Each segment matches at most one sub-edge,
            //   chosen by minimum sum of endpoint distances. Works
            //   well when each cap face's outer loop has exactly one
            //   sub-edge per segment (cube_uniform).
            //
            // * `allow_multi_sub_edge = true` (Stage 4.2):
            //   midpoint-projection matching. For each sub-edge,
            //   compute its midpoint and project onto each segment;
            //   match to the segment with minimum projection distance.
            //   Robust against multi-sub-edge topologies where a
            //   single FilletEnd/FaceEdge has been split into 3+
            //   sub-edges by sliver-absorb. Endpoint-sum fails for
            //   short sub-edges in segment interior (both endpoints
            //   far from segment endpoints) — projection works
            //   because sub-edges geometrically lie on (or near) the
            //   parent fillet cylinder, whose tangent circle = the
            //   FilletEnd segment.
            let mut best: Option<(usize, f64, bool)> = None;
            if constrained_by_cap_index {
                // Route by semantic cap index
                // + twin surface kind. sub_face owns segments
                // {2j, 2j+1}; twin face kind picks parity.
                let semantic_cap_k = semantic_cap_indices
                    .get(cap_k)
                    .and_then(|cap| *cap)
                    .unwrap_or(cap_k);
                let twin_face_data = twin_face.and_then(|f| brep.faces.get(f));
                let twin_kind = twin_face_data
                    .and_then(|f| f.surface.as_ref())
                    .map(|sof| &sof.geometry);
                let target_seg = match twin_kind {
                    Some(SurfaceGeometry::Cylinder { .. }) => Some(2 * semantic_cap_k),
                    Some(SurfaceGeometry::Plane { .. }) => Some(2 * semantic_cap_k + 1),
                    _ => None,
                };
                if let Some(seg_i) = target_seg
                    && seg_i < segment_curves.len()
                {
                    let segment = segment_curves[seg_i];
                    let p_mid = Point3::new(
                        (p_origin.x + p_dest.x) * 0.5,
                        (p_origin.y + p_dest.y) * 0.5,
                        (p_origin.z + p_dest.z) * 0.5,
                    );
                    let t_mid = project_point_on_curve(segment, p_mid);
                    let proj_mid = segment.point_at(t_mid);
                    let d = proj_mid.distance_to(&p_mid);
                    let t_origin = project_point_on_curve(segment, p_origin);
                    let t_dest = project_point_on_curve(segment, p_dest);
                    let reversed = t_origin > t_dest;
                    best = Some((seg_i, d, reversed));
                } else {
                    let p_mid = Point3::new(
                        (p_origin.x + p_dest.x) * 0.5,
                        (p_origin.y + p_dest.y) * 0.5,
                        (p_origin.z + p_dest.z) * 0.5,
                    );
                    for seg_i in [2 * semantic_cap_k, 2 * semantic_cap_k + 1] {
                        let Some(segment) = segment_curves.get(seg_i) else {
                            continue;
                        };
                        let t_mid = project_point_on_curve(segment, p_mid);
                        let proj_mid = segment.point_at(t_mid);
                        let d = proj_mid.distance_to(&p_mid);
                        let t_origin = project_point_on_curve(segment, p_origin);
                        let t_dest = project_point_on_curve(segment, p_dest);
                        let reversed = t_origin > t_dest;
                        if best.map(|(_, bd, _)| d < bd).unwrap_or(true) {
                            best = Some((seg_i, d, reversed));
                        }
                    }
                }
            } else if allow_multi_sub_edge {
                let p_mid = Point3::new(
                    (p_origin.x + p_dest.x) * 0.5,
                    (p_origin.y + p_dest.y) * 0.5,
                    (p_origin.z + p_dest.z) * 0.5,
                );
                for (i, segment) in segment_curves.iter().enumerate() {
                    let t_mid = project_point_on_curve(segment, p_mid);
                    let proj_mid = segment.point_at(t_mid);
                    let d = proj_mid.distance_to(&p_mid);
                    // Determine direction by projecting both endpoints.
                    let t_origin = project_point_on_curve(segment, p_origin);
                    let t_dest = project_point_on_curve(segment, p_dest);
                    let reversed = t_origin > t_dest;
                    if best.map(|(_, bd, _)| d < bd).unwrap_or(true) {
                        best = Some((i, d, reversed));
                    }
                }
            } else {
                for (i, eps) in segment_endpoints.iter().enumerate() {
                    if consumed[i] {
                        continue;
                    }
                    let d_fwd = p_origin.distance_to(&eps[0]) + p_dest.distance_to(&eps[1]);
                    let d_rev = p_origin.distance_to(&eps[1]) + p_dest.distance_to(&eps[0]);
                    let (d, reversed) = if d_fwd <= d_rev {
                        (d_fwd, false)
                    } else {
                        (d_rev, true)
                    };
                    if best.map(|(_, bd, _)| d < bd).unwrap_or(true) {
                        best = Some((i, d, reversed));
                    }
                }
            }
            if trace_match {
                let chord = p_origin.distance_to(&p_dest);
                let twin_kind = brep
                    .faces
                    .get(twin_face.unwrap_or(sub_face))
                    .and_then(|f| f.surface.as_ref())
                    .map(|sof| match &sof.geometry {
                        SurfaceGeometry::Cylinder { .. } => "Cyl",
                        SurfaceGeometry::Plane { .. } => "Plane",
                        SurfaceGeometry::Sphere { .. } => "Sphere",
                        SurfaceGeometry::Gregory(_) => "Gregory",
                        _ => "Other",
                    })
                    .unwrap_or("none");
                let cap_idx = sub_faces
                    .iter()
                    .position(|f| *f == sub_face)
                    .unwrap_or(usize::MAX);
                eprintln!(
                    "[wire-match] cap[{cap_idx}]={sub_face:?} he={he_handle:?} edge={:?} twin={twin_kind} chord={chord:.3e} v_origin=({:.3},{:.3},{:.3}) v_dest=({:.3},{:.3},{:.3}) best_seg={:?} best_d={:?}",
                    he.edge,
                    p_origin.x,
                    p_origin.y,
                    p_origin.z,
                    p_dest.x,
                    p_dest.y,
                    p_dest.z,
                    best.map(|(i, _, _)| i),
                    best.map(|(_, d, _)| d),
                );
            }
            let Some((seg_idx, d, reversed)) = best else {
                continue;
            };
            if d > tol {
                if trace_match {
                    eprintln!("[wire-match]   REJECTED seg={seg_idx} d={d:.3e} > tol={tol:.3e}");
                }
                continue;
            }
            if trace_match {
                eprintln!("[wire-match]   ACCEPTED seg={seg_idx} d={d:.3e} reversed={reversed}");
            }
            if !allow_multi_sub_edge {
                consumed[seg_idx] = true;
            }

            let eps = segment_endpoints[seg_idx];
            let (target_p_origin, target_p_dest) = if reversed {
                (eps[1], eps[0])
            } else {
                (eps[0], eps[1])
            };
            worklist.push((he_handle, seg_idx, reversed, target_p_origin, target_p_dest));
        }
    }

    worklist
}
