//! Multi-chord face splits for Phase B (Round 5r+).
//!
//! When a single original face carries multiple fillet trim chords
//! whose 3D lines intersect in the face interior (e.g. all four
//! chords on a cube face filleted on every edge), a naive
//! sequential `split_face` approach fails on the second chord
//! because its anchor pair has been separated onto two different
//! sub-faces by the first chord's split.
//!
//! This module handles that case with a lazy "insert crossing
//! vertex, then split in sub-chord segments" algorithm:
//!
//! 1. Precompute all pairwise segment intersections among the
//!    trim chords on this face.
//! 2. Apply trims in sequence. For each trim, before its own
//!    `split_face`, walk its precomputed intersection list and
//!    — for every *earlier* trim whose chord crossed this one —
//!    locate the matching sub-edge of the earlier trim in the
//!    current BRep and split it at the crossing point. The new
//!    vertex becomes an intermediate anchor on the current trim.
//! 3. Apply `split_face` on each sub-chord segment in turn. With
//!    the intermediate vertices now on the earlier trims' sub-face
//!    boundaries, each sub-chord has both endpoints on one
//!    sub-face and the split succeeds.
//! 4. Pick the LONGEST sub-chord's consumed face as the "true"
//!    strip for this trim. On a cube's inset-square decomposition
//!    the longest sub-chord is the middle chord between the two
//!    inset corners, and its consumed half is the rectangular
//!    fillet strip. Shorter outer sub-chords produce corner
//!    regions that become sphere vertex blends later.

use super::helpers;
use super::multi_chord_span_geometry::{
    find_sub_edge_containing_point, intersect_segments, remap_split, trim_subcurve_between,
};
use super::phase_split::{FaceTrim, TrimSplitResult, find_face_containing_vertices};
use crate::error::BlendFailure;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
/// Outcome of applying all trims on one original face via
/// `apply_multi_chord_for_face`.
pub(super) struct MultiChordOutcome {
    pub trim_results: Vec<TrimSplitResult>,
    pub failures: Vec<(EdgeHandle, BlendFailure)>,
    /// All new descendant edges produced by `split_edge` calls
    /// inside the multi-chord algorithm (splitting an earlier
    /// trim's sub-chord at a crossing point). The caller folds
    /// these into its fillet-edge descendants map, so Phase C
    /// sees them as valid descendants when identifying consumed
    /// faces.
    pub extra_edge_splits: Vec<(EdgeHandle, EdgeHandle, EdgeHandle)>,
    /// Trims resolved through the PK-141 region path — their
    /// consumed faces are merged regions with no split ancestry, so
    /// the caller's longest-descendant re-pick must skip them.
    pub region_resolved: Vec<(usize, u8)>,
}
/// Apply all of a face's trims, honouring chord-chord
/// intersections. This is the Round 5r+ replacement for the
/// per-face inner loop in Phase B.
///
/// `fillet_edges_by_trim[k]` is the canonical fillet edge for
/// `trims[k]`. `descendants` is the per-fillet-edge descendants
/// map — used here in read-only fashion to identify the strip via
/// `helpers::identify_consumed_face_any`.
pub(super) fn apply_multi_chord_for_face(
    brep: &mut BRep,
    original_face: FaceHandle,
    trims: &[FaceTrim],
    fillet_edges_by_trim: &[EdgeHandle],
    descendants: &std::collections::HashMap<EdgeHandle, Vec<EdgeHandle>>,
    tolerance: &Tolerance,
) -> MultiChordOutcome {
    let mut trim_results: Vec<TrimSplitResult> = Vec::new();
    let mut failures: Vec<(EdgeHandle, BlendFailure)> = Vec::new();
    let mut extra_edge_splits: Vec<(EdgeHandle, EdgeHandle, EdgeHandle)> = Vec::new();
    let n = trims.len();
    if n == 0 {
        return MultiChordOutcome {
            trim_results,
            failures,
            extra_edge_splits,
            region_resolved: Vec::new(),
        };
    }
    // Precompute pairwise chord intersections.
    let segments: Vec<(Point3, Point3)> = trims
        .iter()
        .map(|t| {
            (
                brep.vertices[t.anchors.0].position,
                brep.vertices[t.anchors.1].position,
            )
        })
        .collect();
    let mut intersections: Vec<Vec<(f64, Point3, usize)>> = vec![Vec::new(); n];
    let cross_gate = tolerance.linear().max(1e-9) * 1.0e4;
    for i in 0..n {
        for j in (i + 1)..n {
            let (a, b) = segments[i];
            let (c, d) = segments[j];
            if let Some((s, t, pt)) = intersect_segments(a, b, c, d) {
                intersections[i].push((s, pt, j));
                intersections[j].push((t, pt, i));
            } else if let Some((s, t, pt)) = super::multi_chord_span_geometry::intersect_trim_curves(
                &trims[i].trim_3d,
                &trims[j].trim_3d,
                cross_gate,
            ) {
                // Curved-support crossings: the anchor chords pass
                // inside the surface and never intersect; the sampled
                // trim curves carry the true crossing (PK-141).
                intersections[i].push((s, pt, j));
                intersections[j].push((t, pt, i));
            }
        }
    }
    for list in intersections.iter_mut() {
        list.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }
    let mut per_trim_sub_edges: Vec<Vec<EdgeHandle>> = vec![Vec::new(); n];
    let mut sub_faces: Vec<FaceHandle> = vec![original_face];
    let mut deferred: Vec<super::multi_chord_region::DeferredTrim> = Vec::new();
    for i in 0..n {
        let trim = &trims[i];
        let fillet_edge = fillet_edges_by_trim[i];
        let fillet_desc: &[EdgeHandle] = descendants
            .get(&fillet_edge)
            .map(|v| v.as_slice())
            .unwrap_or(std::slice::from_ref(&fillet_edge));
        // Step 1: for every precomputed intersection with an
        // earlier trim `j < i`, split the earlier trim's current
        // sub-edge at the crossing point to create an
        // intermediate vertex.
        let mut intermediate_vertices: Vec<VertexHandle> = Vec::new();
        for &(_, pt, j) in &intersections[i] {
            if j >= i {
                continue;
            }
            let Some((eh, t_on)) = find_sub_edge_containing_point(brep, &per_trim_sub_edges[j], pt)
            else {
                continue;
            };
            let t_clamped = t_on.clamp(0.001, 0.999);
            match brep.split_edge(eh, t_clamped) {
                Ok(split) => {
                    extra_edge_splits.push((eh, split.edge_a, split.edge_b));
                    remap_split(&mut per_trim_sub_edges, eh, split.edge_a, split.edge_b);
                    intermediate_vertices.push(split.vertex);
                }
                Err(_) => {
                    // Fall through — the sub-chord below will
                    // still be attempted without this crossing
                    // vertex and likely fail cleanly.
                }
            }
        }
        // Step 2: build the ordered vertex sequence for this
        // trim and apply one `split_face` per sub-chord.
        let mut endpoints: Vec<VertexHandle> = Vec::with_capacity(intermediate_vertices.len() + 2);
        endpoints.push(trim.anchors.0);
        endpoints.extend(intermediate_vertices.iter().copied());
        endpoints.push(trim.anchors.1);
        let mut produced_edges: Vec<EdgeHandle> = Vec::new();
        // Candidates: (length, consumed, pcurve, trim_3d).
        // We pick the LONGEST sub-chord's consumed as the true
        // strip for this fillet edge; shorter outer sub-chords
        // produce corner regions.
        let mut strip_candidates: Vec<(f64, FaceHandle, NurbsCurve, NurbsCurve)> = Vec::new();
        let mut interior_candidates: Vec<(f64, NurbsCurve, NurbsCurve)> = Vec::new();
        let mut sub_chord_failed = false;
        for k in 0..endpoints.len() - 1 {
            let mut start_v = endpoints[k];
            let mut end_v = endpoints[k + 1];
            let mut target = find_face_containing_vertices(brep, &sub_faces, start_v, end_v);
            // Chain-holder transplant (PK-141): a trim's cut belongs
            // on the face carrying its own fillet's chain. When both
            // anchors are valence-1 spur tips whose loop face does
            // NOT hold a descendant while another sub-face does (an
            // earlier trim's split left the zero-width spur
            // appendages on the wrong side — measured: the -y
            // vertical wall trim's split_face returned a 2-sided
            // banana inside the kept piece, the whole unbounded-
            // region cascade), transplant both anchors onto the
            // holder and cut there.
            let holder = sub_faces.iter().copied().find(|&f| {
                super::merge::face_loop_contains_any_edge(brep, f, fillet_desc).unwrap_or(false)
            });
            let target_holds = target
                .map(|t| {
                    super::merge::face_loop_contains_any_edge(brep, t, fillet_desc).unwrap_or(false)
                })
                .unwrap_or(false);
            if let Some(h) = holder
                && target != Some(h)
                && !target_holds
                && let Some((nv_start, nv_end)) =
                    super::multi_chord_anchor_transplant::transplant_spur_anchor_pair(
                        brep, start_v, end_v, h, tolerance,
                    )
            {
                if std::env::var_os("PK_MULTICHORD_TRACE").is_some() {
                    eprintln!(
                        "[mc-target] eg{} side{} sub {k}: transplanted anchors \
                         {start_v:?}/{end_v:?} to chain holder {h:?}",
                        trim.edge_idx, trim.side
                    );
                }
                start_v = nv_start;
                end_v = nv_end;
                endpoints[k] = nv_start;
                endpoints[k + 1] = nv_end;
                target = Some(h);
            }
            // PK-141 spur re-home (see rehome_spur_tip).
            if target.is_none() {
                if let Some((nv, tf)) = super::multi_chord_span_geometry::rehome_spur_tip(
                    brep, end_v, start_v, &sub_faces,
                ) {
                    end_v = nv;
                    endpoints[k + 1] = nv;
                    target = Some(tf);
                } else if let Some((nv, tf)) = super::multi_chord_span_geometry::rehome_spur_tip(
                    brep, start_v, end_v, &sub_faces,
                ) {
                    start_v = nv;
                    endpoints[k] = nv;
                    target = Some(tf);
                }
            }
            let sub_midpoint = {
                let p_a = brep.vertices[start_v].position;
                let p_b = brep.vertices[end_v].position;
                Some(Point3::new(
                    (p_a.x + p_b.x) * 0.5,
                    (p_a.y + p_b.y) * 0.5,
                    (p_a.z + p_b.z) * 0.5,
                ))
            };
            let Some(target) = target else {
                if std::env::var_os("PK_MULTICHORD_TRACE").is_some() {
                    let pa = brep.vertices.get(start_v).map(|v| v.position);
                    let pb = brep.vertices.get(end_v).map(|v| v.position);
                    let faces_of = |v: popkern_core::topology::VertexHandle| -> Vec<popkern_core::topology::FaceHandle> {
                        let mut out = std::collections::BTreeSet::new();
                        for (_, he) in brep.half_edges.iter() {
                            if he.origin == v
                                && let Some(l) = brep.loops.get(he.loop_)
                            {
                                out.insert(l.face);
                            }
                        }
                        out.into_iter().collect()
                    };
                    eprintln!(
                        "[mc] sub-chord {start_v:?}@{pa:?} faces {:?} -> {end_v:?}@{pb:?} faces {:?}; sub_faces {sub_faces:?}",
                        faces_of(start_v),
                        faces_of(end_v)
                    );
                }
                failures.push((
                    fillet_edge,
                    BlendFailure::InvalidTopology {
                        reason: format!(
                            "multi-chord sub-chord anchors not found on any sub-face of {original_face:?}"
                        ),
                        location: sub_midpoint,
                    },
                ));
                sub_chord_failed = true;
                break;
            };
            let p_a = brep.vertices[start_v].position;
            let p_b = brep.vertices[end_v].position;
            // #136 slice 6 / PK-139: the sub-chord's 3-D trim is the
            // SUB-SPAN of the fillet's actual trim curve, not the
            // straight chord between the anchors — a by-points
            // variable law's rail bulges between its breakpoints and
            // the chord line cut straight through it (measured: rail
            // Line at offset 0.05 under a 0.12 bulge apex, 40 welded
            // opens). Straight trims still canonicalize back to
            // analytic lines inside `split_face`. Falls back to the
            // chord when the anchors do not invert onto the trim
            // (corner sub-chords off the trim span).
            let t3d_raw = trim_subcurve_between(&trim.trim_3d, &p_a, &p_b, tolerance);
            // Vertex-onto-trim snap: the sliced trim's ends are exact
            // geometry; anchor vertices carry perp-walk noise. Move
            // the vertices (and their incident chords) onto the trim
            // ends so split_face's endpoint matching AND the curve
            // stay machine-exact.
            if let Some(t) = t3d_raw.as_ref() {
                use popkern_core::geometry::CurveEval;
                let (tt0, tt1) = t.parameter_range();
                let e0 = t.point_at(tt0);
                let e1 = t.point_at(tt1);
                let (v_first, v_last) = if e0.distance_to(&p_a) <= e0.distance_to(&p_b) {
                    (start_v, end_v)
                } else {
                    (end_v, start_v)
                };
                for (vh, target) in [(v_first, e0), (v_last, e1)] {
                    let d = brep
                        .vertices
                        .get(vh)
                        .map(|v| v.position.distance_to(&target))
                        .unwrap_or(0.0);
                    if d > tolerance.linear() && d < 0.1 {
                        if let Some(v) = brep.vertices.get_mut(vh) {
                            v.position = target;
                        }
                        crate::construct::topology_repair::snap_incident_edge_curve_endpoints_robust(
                            brep, vh, target,
                        );
                    }
                }
            }
            let t3d = t3d_raw
                .or_else(|| {
                    // WHOLE-TRIM identity: when the anchors ARE the
                    // trim's own endpoints, the sub-span is the trim
                    // itself — the subcurve extraction can fail on a
                    // near-full-span inversion, and the chord fallback
                    // then cuts across territory the trim routes
                    // around (measured: the notched flange's
                    // interrupted-rim trim spans ~331°; the chord
                    // crossed the notch void and the split minted a
                    // phantom top-face fragment inside it).
                    use popkern_core::geometry::curve::CurveEval;
                    let (t0, t1) = trim.trim_3d.parameter_range();
                    let s0 = trim.trim_3d.point_at(t0);
                    let s1 = trim.trim_3d.point_at(t1);
                    // Terminal-construction noise between Phase-A
                    // anchors and the trim's own row ends measures
                    // ~2e-2 on the notched flange; gate relative to
                    // the trim's span.
                    let span = s0
                        .distance_to(&s1)
                        .max(trim.trim_3d.point_at(0.5 * (t0 + t1)).distance_to(&s0));
                    let tol_d = (tolerance.linear() * 100.0).max(span * 0.02);
                    let forward = s0.distance_to(&p_a) < tol_d && s1.distance_to(&p_b) < tol_d;
                    let reverse = s0.distance_to(&p_b) < tol_d && s1.distance_to(&p_a) < tol_d;
                    (forward || reverse).then(|| {
                        // Snap the trim's endpoint CPs onto the anchor
                        // vertices (terminal noise ~2e-2); split_face
                        // requires exact endpoint-vertex identity, and
                        // the end-of-pipeline reconcile re-trims the
                        // authored edge back onto the exact analytic
                        // intersection.
                        let mut t = trim.trim_3d.clone();
                        let (want_first, want_last) = if forward { (p_a, p_b) } else { (p_b, p_a) };
                        if let Some(first) = t.control_points.first_mut() {
                            *first = want_first;
                        }
                        if let Some(last) = t.control_points.last_mut() {
                            *last = want_last;
                        }
                        t
                    })
                })
                .unwrap_or_else(|| helpers::line_nurbs(&p_a, &p_b));
            let pc = match helpers::pcurve_for_face_along(brep, target, &t3d) {
                Ok(pc) => pc,
                Err(e) => {
                    failures.push((
                        fillet_edge,
                        BlendFailure::InvalidTopology {
                            reason: e,
                            location: sub_midpoint,
                        },
                    ));
                    sub_chord_failed = true;
                    break;
                }
            };
            let split = match brep.split_face(target, &pc, &t3d, tolerance) {
                Ok(s) => s,
                Err(e) => {
                    failures.push((
                        fillet_edge,
                        BlendFailure::InvalidTopology {
                            reason: format!("face split failed: {e}"),
                            location: sub_midpoint,
                        },
                    ));
                    sub_chord_failed = true;
                    break;
                }
            };
            produced_edges.extend(split.new_edges.iter().copied());
            // identify_consumed_face_any matches any descendant.
            // MID-SPAN descendants take priority: corner/seam
            // remnant segments also carry descendant status, and
            // identification through a remnant picked the rim-side
            // piece (measured: 0.7 rim deviation).
            let mid_span_desc = super::multi_chord_region_seeding::mid_span_descendants(
                brep,
                fillet_desc,
                &trim.trim_3d,
            );
            let identify_desc: &[EdgeHandle] = if mid_span_desc.is_empty() {
                fillet_desc
            } else {
                &mid_span_desc
            };
            let (kept, consumed) = match helpers::identify_consumed_face_any(
                brep,
                split.face_a,
                split.face_b,
                identify_desc,
            ) {
                Ok(kc) => kc,
                Err(_) => {
                    if std::env::var_os("PK_MULTICHORD_TRACE").is_some() {
                        eprintln!(
                            "[mc] sub-chord {k} on {target:?}: neither half ({:?}/{:?}) carries a descendant of {fillet_edge:?} (holders: {:?})",
                            split.face_a,
                            split.face_b,
                            super::multi_chord_span_geometry::descendant_holders(brep, fillet_desc)
                        );
                    }
                    // Neither half contains any descendant —
                    // this sub-chord is interior; keep its curves
                    // for the region fallback (PK-141).
                    interior_candidates.push((p_a.distance_to(&p_b), pc.clone(), t3d.clone()));
                    if let Some(pos) = sub_faces.iter().position(|&f| f == target) {
                        sub_faces[pos] = split.face_a;
                        sub_faces.push(split.face_b);
                    }
                    continue;
                }
            };
            if let Some(pos) = sub_faces.iter().position(|&f| f == target) {
                sub_faces[pos] = kept;
                sub_faces.push(consumed);
            }
            let length = p_a.distance_to(&p_b);
            strip_candidates.push((length, consumed, pc, t3d));
        }
        if sub_chord_failed {
            per_trim_sub_edges[i] = produced_edges;
            continue;
        }
        per_trim_sub_edges[i] = produced_edges.clone();
        // Pick the LONGEST sub-chord's consumed face as the true
        // strip for this fillet edge.
        strip_candidates.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
        let Some((_, consumed, pcurve, trim_3d)) = strip_candidates.into_iter().next() else {
            // PK-141: defer to the region fallback.
            interior_candidates
                .sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
            deferred.push(super::multi_chord_region::DeferredTrim {
                edge_idx: trim.edge_idx,
                side: trim.side,
                fillet_edge,
                candidate: interior_candidates.into_iter().next(),
                produced_edges,
            });
            continue;
        };
        trim_results.push(TrimSplitResult {
            edge_idx: trim.edge_idx,
            side: trim.side,
            consumed,
            trim_edges: produced_edges,
            pcurve,
            trim_3d,
        });
    }
    super::multi_chord_region::reroute_all_when_deferred(
        &mut deferred,
        &mut trim_results,
        trims,
        fillet_edges_by_trim,
    );
    let all_trim_edges: Vec<EdgeHandle> = per_trim_sub_edges.iter().flatten().copied().collect();
    let pre_deferred_len = trim_results.len();
    super::multi_chord_region::process_deferred(
        brep,
        deferred,
        trims,
        descendants,
        &sub_faces,
        &all_trim_edges,
        &mut trim_results,
        &mut failures,
        tolerance,
    );
    let region_resolved: Vec<(usize, u8)> = trim_results[pre_deferred_len..]
        .iter()
        .map(|r| (r.edge_idx, r.side))
        .collect();
    MultiChordOutcome {
        trim_results,
        failures,
        extra_edge_splits,
        region_resolved,
    }
}
