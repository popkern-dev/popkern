//! Cellular boolean: the primary boolean pipeline using cell decomposition.
//!
//! Reuses the intersection graph infrastructure (phases 1-5 of graph-merge)
//! to produce face fragments, then classifies fragments via point-in-solid
//! to determine cell membership, selects cells by operation predicate, and
//! derives fragment selections from the cell selection.
//!
//! The assembly step reuses graph-merge's `assemble_result`, so the only
//! difference from graph-merge is HOW fragments are selected — cellular
//! selection is order-independent because it classifies volumetric cells
//! rather than individual fragments.

use super::fragment_interior::{fragment_interior_sample, planar_outward_normal};

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::{BoolDiagnostics, BoolOperation, PointClassification};
use crate::intersect::IntersectionResult;
use crate::merge::fragment_classify::FragmentSelection;
use crate::merge::graph::{IntersectionGraph, Operand};
use crate::merge::graph_build::insert_operand_topology;
use crate::merge::graph_intersect::insert_intersection_curves;
use crate::merge::graph_split::{
    EndpointClass, classify_endpoint, remove_collinear_duplicates, split_edges_at_shared_nodes,
    split_graph_edge, split_shared_edges_at_nodes,
};
use crate::merge::radial::sort_all_radial;
use crate::merge::walk::{Fragment, walk_all_fragments_ext};

/// Execute a boolean operation using the cellular decomposition approach.
///
/// Phases 1-5 (graph build, edge-face splits, intersection curves, radial
/// sort, fragment walk) are identical to graph-merge. Phase 6 replaces
/// per-fragment ray casting with cell-based classification:
///
/// 1. For each fragment, classify a sample point against both bodies
/// 2. Derive a cell membership (in A, in B, in both, in neither)
/// 3. Select cells by the operation predicate
/// 4. Keep fragments whose cell is selected
///
/// This is order-independent because the classification depends only on
/// geometry (point-in-solid), not on iteration order.
pub(crate) fn cellular_graph_merge(
    original_a: &BRep,
    original_b: &BRep,
    intersection_result: &IntersectionResult,
    op: BoolOperation,
    tolerance: &Tolerance,
) -> Result<BRep, BoolDiagnostics> {
    let snap_tol = tolerance.linear() * 2.0;
    let mut graph = IntersectionGraph::new();

    // Phase 1: Insert original topology from both operands.
    let (_vmap_a, _hemap_a) = insert_operand_topology(&mut graph, original_a, Operand::A, snap_tol);
    let (_vmap_b, _hemap_b) = insert_operand_topology(&mut graph, original_b, Operand::B, snap_tol);

    // Phase 2: Process edge-face hits.
    if std::env::var_os("PK_BOOL_TRACE").is_some() {
        eprintln!(
            "  [cellular] {} edge-face hits",
            intersection_result.edge_face_hits.len()
        );
        for hit in &intersection_result.edge_face_hits {
            eprintln!(
                "    hit at ({:.3},{:.3},{:.3})",
                hit.point.x, hit.point.y, hit.point.z
            );
        }
    }
    let mut hit_nodes: Vec<usize> = Vec::new();
    for hit in &intersection_result.edge_face_hits {
        let endpoint_class = classify_endpoint(&graph, &hit.point, snap_tol);
        if let EndpointClass::OnEdge(he_idx, t_param) = endpoint_class {
            hit_nodes.push(split_graph_edge(&mut graph, he_idx, hit.point, t_param));
        }
    }

    // Phase 3: Insert intersection curves as shared edges, then
    // split them at edge-face hit nodes they pass through (the loop
    // walk must be able to route through those junctions).
    // Closed intersection curves must CLOSE at an edge-face hit
    // point when one lies on them: the loop walk enters/bounces off
    // an intersection circle only at its closure node (§7.4), and
    // the analytic cylinder-hole cases have always closed at the
    // cylinder seam — the crossing — by construction. A circle from
    // an offset walk start (or a swept boss whose seam is elsewhere)
    // must be rotated to that convention.
    let rotated: Vec<crate::intersect::IntersectionRecord> = intersection_result
        .records
        .iter()
        .map(|r| rotate_closed_record_to_hits(r, &intersection_result.edge_face_hits, tolerance))
        .collect();
    insert_intersection_curves(&mut graph, &rotated, snap_tol);
    let _ = &hit_nodes;
    // Split SHARED edges at record×record T-junctions: one record's
    // endpoint frequently lands mid-span of another (a dovetail
    // rail's flank trace ends ON the joint-shoulder trace), and the
    // walk cannot route the junction without a node-split there
    // (measured: the rail's floor strip never became a fragment and
    // the carve went non-manifold). Every node is a candidate — the
    // splitter skips nodes not on a shared edge's interior.
    let all_nodes: Vec<usize> = (0..graph.nodes.len()).collect();
    split_shared_edges_at_nodes(&mut graph, &all_nodes, snap_tol);
    if std::env::var_os("PK_BOOL_TRACE").is_some() {
        let mut i = 0;
        while i < graph.half_edges.len() {
            let he = &graph.half_edges[i];
            if i < he.twin && matches!(he.source, crate::merge::graph::EdgeSource::Shared { .. }) {
                let a = graph.nodes[he.origin].position;
                let b = graph.nodes[graph.half_edges[he.twin].origin].position;
                eprintln!(
                    "  [graph] shared he {i}: ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) curve={}",
                    a.x,
                    a.y,
                    a.z,
                    b.x,
                    b.y,
                    b.z,
                    he.curve.is_some()
                );
            }
            i += 2;
        }
    }

    // Phase 3.1–3.2: Handle coincident face boundaries; also split
    // graph edges at nodes created by edge-face hits — an
    // intersection curve passing THROUGH such a node (e.g. a seam
    // edge crossing a boss/plate intersection circle) must get a
    // node there or the loop walk cannot route through it.
    if !intersection_result.edge_face_hits.is_empty()
        || intersection_result
            .coincidence
            .face_overlaps
            .iter()
            .any(|o| !o.full_overlap)
    {
        split_edges_at_shared_nodes(&mut graph, snap_tol);
        if intersection_result
            .coincidence
            .face_overlaps
            .iter()
            .any(|o| !o.full_overlap)
        {
            remove_collinear_duplicates(&mut graph);
        }
    }

    // Phase 4: Sort radial edges at all nodes.
    sort_all_radial(&mut graph);

    // Phase 5: Walk loops to find face fragments.
    let overlap_pairs: Vec<_> = intersection_result
        .coincidence
        .face_overlaps
        .iter()
        .map(|o| (o.face_a, o.face_b))
        .collect();
    let mut walk_fallbacks: Vec<(crate::merge::graph::FaceRef, String)> = Vec::new();
    let fragments = walk_all_fragments_ext(
        &graph,
        original_a,
        original_b,
        !intersection_result.coincidence.face_overlaps.is_empty(),
        true,
        &overlap_pairs,
        &mut walk_fallbacks,
    );

    // Phase 6 (CELLULAR): Classify fragments by cell membership.
    let selections = classify_fragments_cellular(
        &fragments,
        &graph,
        original_a,
        original_b,
        op,
        tolerance,
        &intersection_result.coincidence.face_overlaps,
    )?;

    // Phase 7: Order-independent assembly from selected fragments.
    let kept: Vec<usize> = (0..fragments.len())
        .filter(|&i| selections[i].keep)
        .collect();

    let result = super::assemble::cellular_assemble(
        &graph,
        &fragments,
        &selections,
        original_a,
        original_b,
        *tolerance,
    )
    .map_err(|mut diag| {
        // ChartSubdivisionReport promotion (PK-106): a failed cellular
        // boolean carries, per face that could not chart-walk, WHY it
        // took the legacy fallback — the substrate defects are usually
        // the trailhead for the assembly failure being reported.
        for (face_ref, reason) in &walk_fallbacks {
            diag.add(crate::error::DiagnosticEntry {
                severity: crate::error::Severity::Warning,
                phase: crate::error::Stage::Merge,
                message: format!("chart walk fell back to legacy on {face_ref:?}: {reason}"),
                entities: vec![crate::error::EntityRef::Face(face_ref.face)],
                measured_value: None,
                location: None,
                remediation: None,
            });
        }
        diag
    })?;

    // Build provenance map (available via boolean_with_provenance).
    let _provenance =
        super::assemble_cleanup::build_face_provenance(&result, &kept, &fragments, &selections);

    Ok(result)
}

/// Classify fragments using cellular logic: determine which volumetric
/// cell each fragment borders, then select based on cell membership.
///
/// For each fragment:
/// 1. Sample a point on the fragment
/// 2. Test it against BOTH bodies via point-in-solid
/// 3. Compute cell membership: bits 0=in_A, 1=in_B
/// 4. Apply operation predicate to the cell membership
#[allow(clippy::too_many_arguments)]
fn classify_fragments_cellular(
    fragments: &[Fragment],
    graph: &IntersectionGraph,
    brep_a: &BRep,
    brep_b: &BRep,
    op: BoolOperation,
    tolerance: &Tolerance,
    face_overlaps: &[crate::coincidence::types::FaceOverlap],
) -> Result<Vec<FragmentSelection>, BoolDiagnostics> {
    use crate::adaptive_precision::robust_point_in_solid_final;
    use crate::intersect::compute_face_aabbs;

    let aabbs_a = compute_face_aabbs(brep_a);
    let aabbs_b = compute_face_aabbs(brep_b);

    // Faces participating in a COINCIDENT overlap pair: fragments of
    // these lie in the other body's boundary plane, so a raw sample
    // has ill-defined ray parity (measured: two mirror-image cap
    // fragments classified Inside on one side and Outside on the
    // other — the kept duplicate skin went non-manifold). They are
    // routed to the offset-pair classifier, whose (in, out) verdict
    // maps to OnBoundary and the operation's keep-one-skin rules.
    let overlap_a: std::collections::HashSet<_> =
        face_overlaps.iter().map(|fo| fo.face_a).collect();
    let overlap_b: std::collections::HashSet<_> =
        face_overlaps.iter().map(|fo| fo.face_b).collect();

    // Probe offset for side/coorientation sampling. It must sit
    // decades ABOVE the linear tolerance band (so boundary parity
    // noise cannot flip a verdict) yet decades BELOW feature scale
    // (so the probe point stays within the local cell). The old
    // `linear × 1e3` scaling honored only the first bound: at a
    // coarse linear tolerance of 1e-4 the ±16× retry probes stepped
    // 1.6 — clear outside a unit-sized model — and every
    // coincident-skin verdict degraded to parity noise (measured:
    // flush cube∪cube at 1e-4 kept both overlap skins, χ = −10,
    // 16 open edges). Cap by the combined model diagonal; the floor
    // stays above the tolerance band even for tiny models.
    let model_diag = {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for brep in [brep_a, brep_b] {
            for (_, v) in &brep.vertices {
                for (i, c) in [v.position.x, v.position.y, v.position.z]
                    .into_iter()
                    .enumerate()
                {
                    lo[i] = lo[i].min(c);
                    hi[i] = hi[i].max(c);
                }
            }
        }
        let d = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    };
    let probe_eps = if model_diag.is_finite() && model_diag > 0.0 {
        (tolerance.linear() * 1.0e3)
            .max(1.0e-3)
            .min(model_diag * 2.0e-3)
            .max(tolerance.linear() * 8.0)
    } else {
        (tolerance.linear() * 1.0e3).max(1.0e-3)
    };
    // The coarser retry scale (corner-graze disambiguation) obeys the
    // same upper bound: an unclamped ×16 retry at linear=1e-2 probed
    // 1.28 — half the model away — read the far side of the solid,
    // and the trust-the-coarser-scale rule overrode a CORRECT
    // fine-scale verdict (measured: cube∪sphere at 1e-2 dropped the
    // sphere wall, χ = 8). When the clamp collapses the ladder the
    // retry degenerates to the fine scale, which is simply a no-op.
    let probe_eps_coarse = if model_diag.is_finite() && model_diag > 0.0 {
        (probe_eps * 16.0).min(model_diag * 5.0e-2).max(probe_eps)
    } else {
        probe_eps * 16.0
    };

    let mut selections = Vec::with_capacity(fragments.len());

    for fragment in fragments {
        // DEGENERATE loops (a lone non-closed half-edge, or a
        // he/twin bigon) enclose no area — boundary-graze residue.
        // Keeping one produced a dangling face in the assembly.
        let degenerate = match fragment.outer_loop.len() {
            1 => {
                let he = &graph.half_edges[fragment.outer_loop[0]];
                if he.origin != graph.half_edges[he.twin].origin {
                    true
                } else {
                    // A lone SELF-LOOP is a real fragment boundary when
                    // it is a closed intersection curve (a circle cut
                    // through a face). A ZERO-EXTENT self-loop is a
                    // tangency-point record (a surface grazing another
                    // at one point) — keeping it assembled a point-face
                    // that broke the solid's Euler count. Cull by
                    // sampled curve extent against the tolerance.
                    let extent = he
                        .curve
                        .as_ref()
                        .or_else(|| graph.half_edges[he.twin].curve.as_ref())
                        .map(|c| {
                            use popkern_core::geometry::curve::CurveEval;
                            let origin = graph.nodes[he.origin].position;
                            (1..8)
                                .map(|i| {
                                    let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 8.0;
                                    c.geometry.point_at(t).distance_to(&origin)
                                })
                                .fold(0.0, f64::max)
                        })
                        .unwrap_or(0.0);
                    extent < tolerance.linear() * 10.0
                }
            }
            2 => {
                let a = fragment.outer_loop[0];
                if graph.half_edges[a].twin == fragment.outer_loop[1] {
                    // A he/twin bigon is graze residue only for RECORD
                    // pairs. A BOUNDARY pair walked as one loop is a
                    // periodic seam spine — the outer loop of a full
                    // periodic chart (a full sphere's two meridian
                    // copies bound the whole wall face, with the
                    // imprint loops as its holes).
                    matches!(
                        graph.half_edges[a].source,
                        crate::merge::graph::EdgeSource::Shared { .. }
                    )
                } else {
                    // An OPEN two-half-edge walk with (sampled) zero
                    // area is residue of a seeded walk that swallowed
                    // a tangency self-loop and dead-ended — e.g. the
                    // unconsumed direction of a boundary-collinear
                    // corner record plus the zero-extent point-contact
                    // record at its far end. It can never assemble
                    // (measured: kept as OnBoundary, it became an
                    // unclosed 20 mm loop and the whole flush-arm
                    // chain fell to graph-merge). Real two-arc bigons
                    // (a split circle) have genuine sampled area and
                    // stay.
                    let start = graph.half_edges[fragment.outer_loop[0]].origin;
                    let last = fragment.outer_loop[1];
                    let end = graph.half_edges[graph.half_edges[last].twin].origin;
                    let open = start != end;
                    let area =
                        crate::merge::walk::walk_geometry::loop_newell(graph, &fragment.outer_loop)
                            .length();
                    open && area < (tolerance.linear() * 10.0).powi(2)
                }
            }
            _ => false,
        };
        if degenerate {
            if std::env::var_os("PK_CELL_TRACE").is_some() {
                eprintln!(
                    "  [sel] {:?} outer={} DEGENERATE — dropped",
                    fragment.face,
                    fragment.outer_loop.len()
                );
            }
            selections.push(FragmentSelection {
                keep: false,
                flip: false,
            });
            continue;
        }
        // Sample a point on this fragment. Fragments with inner
        // loops or closed-rim outers need the interior-aware
        // sampler: the edge-midpoint shortcut degenerates to the
        // SEAM VERTEX on a closed rim (maximally unstable for
        // point-in-solid — measured: the bottom cap of a cylinder
        // with a star socket cut through it was dropped and the
        // tool's outside face kept), and a plain centroid of a
        // holed cap lands INSIDE THE HOLE.
        let frag_surface = match fragment.face.operand {
            Operand::A => brep_a.faces.get(fragment.face.face),
            Operand::B => brep_b.faces.get(fragment.face.face),
        }
        .and_then(|f| f.surface.as_ref())
        .map(|s| (s.geometry.clone(), s.orientation));
        // The edge-midpoint shortcut is the tuned default. On a
        // coincident-overlap face it can land ON the other body's
        // boundary EDGE (a flush cutter's flank midpoint sits on the
        // tab's corner line), where the offset-pair probes straddle
        // an edge instead of a face and the coorientation verdict
        // inverts. Resample interior-aware ONLY in that case —
        // blanket interior sampling regressed the flush union, whose
        // arc-bounded overlap fragments sample less reliably through
        // the scanline than through the midpoint shortcut.
        let on_overlap_face_early = match fragment.face.operand {
            Operand::A => overlap_a.contains(&fragment.face.face),
            Operand::B => overlap_b.contains(&fragment.face.face),
        };
        let other_for_edges = match fragment.face.operand {
            Operand::A => brep_b,
            Operand::B => brep_a,
        };
        let shortcut = sample_fragment_point(graph, &fragment.outer_loop);
        let eps_probe = probe_eps;
        // Second trigger: the shortcut point must admit a DECISIVE
        // own-body side probe (one side in, one side out). A shortcut
        // landing on the fragment's own crease edge is indecisive at
        // any offset scale — and every probe made from it, against
        // either body, is noise (measured: a flush cutter's flank
        // fragment probed from the tab's corner edge and the
        // relief-mouth membrane was kept).
        let calibratable_at = |p: &popkern_core::geometry::Point3| -> bool {
            let own_brep_ref = match fragment.face.operand {
                Operand::A => brep_a,
                Operand::B => brep_b,
            };
            let own_aabbs_ref = match fragment.face.operand {
                Operand::A => &aabbs_a,
                Operand::B => &aabbs_b,
            };
            let Some(n0) = planar_outward_normal(own_brep_ref, fragment.face.face) else {
                return true; // no planar normal — legacy path, no opinion
            };
            for e in [eps_probe, probe_eps_coarse] {
                let a = robust_point_in_solid_final(
                    &(*p + n0 * e),
                    own_brep_ref,
                    own_aabbs_ref,
                    tolerance,
                );
                let b = robust_point_in_solid_final(
                    &(*p - n0 * e),
                    own_brep_ref,
                    own_aabbs_ref,
                    tolerance,
                );
                if matches!(
                    (a, b),
                    (PointClassification::Inside, PointClassification::Outside)
                        | (PointClassification::Outside, PointClassification::Inside)
                ) {
                    return true;
                }
            }
            false
        };
        let sample = if on_overlap_face_early {
            // Both candidate samplers can fail on coincident-overlap
            // fragments: the midpoint shortcut lands on crease edges
            // and the other body's corners; the scanline interior
            // sampler mis-rows on arc-heavy charts. Pick whichever
            // admits a DECISIVE own-body side probe — the tuned
            // shortcut first unless it sits in a known hazard — and
            // keep the shortcut when neither calibrates (legacy
            // behavior, orientation stays unknown downstream).
            let interior = frag_surface.as_ref().and_then(|(surf, ori)| {
                fragment_interior_sample(graph, fragment, surf, *ori, true)
            });
            let hazard = point_near_any_edge(other_for_edges, &shortcut, eps_probe * 10.0);
            // HOLED fragments must sample interior-aware regardless
            // of hazard: an all-shared ring fragment's shortcut falls
            // through to the loop-vertex CENTROID, which for an
            // annular region lands in the hole — a point decisively
            // calibratable against the own body yet not ON the
            // fragment at all (measured: a moat cut's tube-material
            // ring on the flush block bottom sampled the tube's void,
            // classified Outside, and the subtract kept the ring as a
            // duplicate skin).
            let ordered = if hazard || !fragment.inner_loops.is_empty() {
                [interior, Some(shortcut)]
            } else {
                [Some(shortcut), interior]
            };
            ordered
                .into_iter()
                .flatten()
                .find(|p| calibratable_at(p))
                .unwrap_or(shortcut)
        } else if let Some(p) = frag_surface
            .as_ref()
            .and_then(|(surf, ori)| fragment_interior_sample(graph, fragment, surf, *ori, false))
        {
            p
        } else {
            // A fragment CUT by the intersection (Shared edge in its
            // loop) can have every surviving boundary edge hugging
            // the other body while only its interior is clear of it
            // (measured on the fanned flush arm: the swoop web's cut
            // curve lies on the plate underside, its flush side edges
            // on the arm walls, and both ends embed in root and
            // plate — which boundary edge the walk enumerated first
            // decided whether the web survived). Vote the non-shared
            // boundary midpoints: unanimity keeps the tuned shortcut
            // (identical to legacy behavior); disagreement means the
            // boundary straddles the other body, and the one sample
            // that is decisive by construction is a point stepped OFF
            // the cut onto the fragment's own side.
            let other = match fragment.face.operand {
                Operand::A => (brep_b, &aabbs_b),
                Operand::B => (brep_a, &aabbs_a),
            };
            let cut = fragment.outer_loop.iter().any(|&h| {
                matches!(
                    graph.half_edges[h].source,
                    crate::merge::graph::EdgeSource::Shared { .. }
                )
            });
            let mids = if cut {
                non_shared_chord_midpoints(graph, &fragment.outer_loop, 6)
            } else {
                Vec::new()
            };
            let unanimous = mids.len() < 2 || {
                let first = robust_point_in_solid_final(&mids[0], other.0, other.1, tolerance);
                matches!(
                    first,
                    PointClassification::Inside | PointClassification::Outside
                ) && mids[1..]
                    .iter()
                    .all(|p| robust_point_in_solid_final(p, other.0, other.1, tolerance) == first)
            };
            if unanimous {
                // The shortcut itself can land ON the other body's
                // boundary when the fragment is a thin CUT sliver —
                // a bore-cap piece beyond a breached wall has only
                // TWO boundary edges (the cut chord + a rim arc), so
                // the disagreement vote above can never fire
                // (mids.len() < 2) and the first-edge shortcut sits
                // exactly on the cutting plane, a parity coin flip
                // (measured: the sliver classified Inside and the
                // subtract kept a flap outside the body). When the
                // shortcut is hazard-positioned and a non-shared
                // boundary midpoint exists, that midpoint is the
                // decisive sample — it lies on the fragment's OWN
                // original boundary, clear of the cut.
                let hazard = cut
                    && !mids.is_empty()
                    && point_near_any_edge(other.0, &shortcut, eps_probe * 10.0);
                if hazard { mids[0] } else { shortcut }
            } else {
                frag_surface
                    .as_ref()
                    .and_then(|surf| {
                        super::fragment_sampling::fragment_cut_offset_sample(
                            graph, fragment, &surf.0,
                        )
                    })
                    .unwrap_or(shortcut)
            }
        };

        // Classify against the OPPOSITE body only.
        let own_operand = fragment.face.operand;
        let (other_brep, other_aabbs) = match own_operand {
            Operand::A => (brep_b, &aabbs_b),
            Operand::B => (brep_a, &aabbs_a),
        };

        // OFFSET-PAIR classification: a fragment can lie exactly on
        // (or extend) the other solid's trim geometry — a socket
        // tool's wall below the part shares its column with the
        // freshly cut hole outline, so every on-face sample grazes
        // the other body's trim and single-point parity is a coin
        // flip (measured: below-plate tool-wall fragments classified
        // Inside and kept, duplicating pocket walls while the true
        // cap was dropped). Classify the pair (sample ± ε·normal)
        // instead: agreement is decisive; disagreement means the
        // fragment sits ON the other boundary.
        // Offset-pair classification is scoped to fragments whose
        // geometry can coincide with the other body's trim curves
        // (closed-rim outers and holed fragments): applying it to
        // every fragment regressed plain cube-notch subtractions
        // (fixtures that pass with single-point classification).
        let on_overlap_face = match fragment.face.operand {
            Operand::A => overlap_a.contains(&fragment.face.face),
            Operand::B => overlap_b.contains(&fragment.face.face),
        };
        let needs_offset_pair = on_overlap_face
            || !fragment.inner_loops.is_empty()
            || fragment.outer_loop.iter().any(|&h| {
                let he = &graph.half_edges[h];
                he.origin == graph.half_edges[he.twin].origin
            });
        // `cooriented`: for OnBoundary fragments with a MATERIAL-TRUE
        // outward normal available, whether the other body's outward
        // normal at the contact agrees with this fragment's. Contact
        // orientation decides keep-one-skin semantics per operation:
        // a COORIENTED subtract contact means the removed region sits
        // exactly behind the fragment (the relief mouth is open — the
        // piece must go), while an ANTI-oriented contact is a tool
        // pressed against the wall from outside (the wall stays).
        // `None` = orientation unknown (no planar outward normal);
        // legacy rules apply.
        let mut cooriented: Option<bool> = None;
        let class = if !needs_offset_pair {
            robust_point_in_solid_final(&sample, other_brep, other_aabbs, tolerance)
        } else {
            use popkern_core::geometry::surface::SurfaceEval;
            let own_brep = match fragment.face.operand {
                Operand::A => brep_a,
                Operand::B => brep_b,
            };
            let frag_face = own_brep.faces.get(fragment.face.face);
            // Prefer the material-true outward normal (the stored
            // outer-loop Newell vector points INTO the material in
            // this kernel; negate it). The chart normal's sign is
            // arbitrary relative to material, which is fine for the
            // agree/disagree verdict but useless for coorientation.
            let outward = planar_outward_normal(own_brep, fragment.face.face);
            let normal = outward.or_else(|| {
                frag_face
                    .and_then(|f| f.surface.as_ref())
                    .and_then(|sof| {
                        popkern_core::geometry::uv_projection::project_point_to_uv(
                            &sof.geometry,
                            &sample,
                        )
                        .map(|(u, v)| sof.geometry.normal_at(u, v))
                    })
                    .and_then(|n| n.normalized())
            });
            match normal {
                Some(n_raw) => {
                    let eps = probe_eps;
                    // SELF-CALIBRATE the normal against the fragment's
                    // OWN body: the loop-Newell sign convention is not
                    // uniform across construction paths (measured: a
                    // unit-cube face's Newell points into the material,
                    // an extruded prism's points out of it). Probe both
                    // sides against the own solid; outward is the side
                    // that reads Outside. An ambiguous self-probe keeps
                    // orientation unknown.
                    let (own_brep_ref, own_aabbs) = match fragment.face.operand {
                        Operand::A => (brep_a, &aabbs_a),
                        Operand::B => (brep_b, &aabbs_b),
                    };
                    let calibrate = |e: f64| {
                        let s_plus = robust_point_in_solid_final(
                            &(sample + n_raw * e),
                            own_brep_ref,
                            own_aabbs,
                            tolerance,
                        );
                        let s_minus = robust_point_in_solid_final(
                            &(sample - n_raw * e),
                            own_brep_ref,
                            own_aabbs,
                            tolerance,
                        );
                        match (s_plus, s_minus) {
                            (PointClassification::Outside, PointClassification::Inside) => {
                                Some(n_raw)
                            }
                            (PointClassification::Inside, PointClassification::Outside) => {
                                Some(n_raw * -1.0)
                            }
                            _ => None,
                        }
                    };
                    // Same corner-noise hazard as the classification
                    // pair: retry an indecisive calibration at the
                    // coarser offset.
                    let n_outward = calibrate(eps).or_else(|| calibrate(probe_eps_coarse));
                    let n = n_outward.unwrap_or(n_raw);
                    let pair_at = |e: f64| {
                        let p_plus = sample + n * e;
                        let p_minus = sample - n * e;
                        let v_plus = robust_point_in_solid_final(
                            &p_plus,
                            other_brep,
                            other_aabbs,
                            tolerance,
                        );
                        let v_minus = robust_point_in_solid_final(
                            &p_minus,
                            other_brep,
                            other_aabbs,
                            tolerance,
                        );
                        // PK_PAIR_TRACE: per-probe forensics for the
                        // offset-pair classifier — prints the exact
                        // probe points and their verdicts (PK_SEL_TRACE
                        // shows only the resolved pair). This is how
                        // the out-of-corridor ×16 retry was caught.
                        if std::env::var_os("PK_PAIR_TRACE").is_some() {
                            eprintln!(
                                "        [pp] e={e:.3e} +({:.3},{:.3},{:.3})={v_plus:?} -({:.3},{:.3},{:.3})={v_minus:?}",
                                p_plus.x, p_plus.y, p_plus.z, p_minus.x, p_minus.y, p_minus.z
                            );
                        }
                        (v_plus, v_minus)
                    };
                    // Probe at the tuned offset AND at 16× it: 1 µm
                    // probes graze corner-adjacent faces when the
                    // sample sits on the other body's edge or corner
                    // fold, and the ray parity inverts the verdict
                    // (measured: a flush cutter's flank fragment read
                    // anti-oriented at its corner sample and the
                    // relief-mouth membrane was kept). When the two
                    // scales disagree, trust the coarser one — still
                    // orders of magnitude below feature size.
                    let (c_plus_f, c_minus_f) = pair_at(eps);
                    let (c_plus, c_minus) = {
                        let (c_plus_c, c_minus_c) = pair_at(probe_eps_coarse);
                        if (c_plus_f, c_minus_f) != (c_plus_c, c_minus_c) {
                            (c_plus_c, c_minus_c)
                        } else {
                            (c_plus_f, c_minus_f)
                        }
                    };
                    if std::env::var_os("PK_SEL_TRACE").is_some() {
                        eprintln!(
                            "      [pair] {:?} outer={} n=({:.2},{:.2},{:.2}) calibrated={} eps={eps:.2e} c_plus={c_plus:?} c_minus={c_minus:?}",
                            fragment.face,
                            fragment.outer_loop.len(),
                            n.x,
                            n.y,
                            n.z,
                            n_outward.is_some()
                        );
                    }
                    match (c_plus, c_minus) {
                        (PointClassification::Inside, PointClassification::Inside) => {
                            PointClassification::Inside
                        }
                        (PointClassification::Outside, PointClassification::Outside) => {
                            PointClassification::Outside
                        }
                        (c_p, c_m) => {
                            if n_outward.is_some() {
                                // Other body's material on the −n side
                                // (behind this fragment's skin) means
                                // its outward normal ≈ +n: cooriented.
                                cooriented = match (c_p, c_m) {
                                    (PointClassification::Outside, PointClassification::Inside) => {
                                        Some(true)
                                    }
                                    (PointClassification::Inside, PointClassification::Outside) => {
                                        Some(false)
                                    }
                                    _ => None,
                                };
                            }
                            PointClassification::OnBoundary
                        }
                    }
                }
                None => robust_point_in_solid_final(&sample, other_brep, other_aabbs, tolerance),
            }
        };

        // Selection rules equivalent to graph-merge but expressed via
        // cellular logic: keep fragment if the cell on its outer side
        // is not selected by the operation predicate.
        if std::env::var_os("PK_CELL_TRACE").is_some() {
            eprintln!(
                "  [sel] {:?} outer={} sample=({:.2},{:.2},{:.2}) class={class:?} co={cooriented:?}",
                fragment.face,
                fragment.outer_loop.len(),
                sample.x,
                sample.y,
                sample.z
            );
        }
        // Keep-one-skin rules with contact ORIENTATION where known:
        // COORIENTED contact (both outward normals agree) means the
        // two boundaries coincide as the same skin — union/intersect
        // keep exactly one copy (operand A's), and subtract drops A's
        // piece (the removed region sits directly behind it, opening
        // the mouth). ANTI-oriented contact means the bodies face
        // each other across the surface — the union of the two
        // materials makes it interior (drop), while a subtract
        // leaves the wall intact (keep). Unknown orientation keeps
        // the legacy cooriented assumption.
        let (keep, flip) = match op {
            BoolOperation::Union => (
                class == PointClassification::Outside
                    || (own_operand == Operand::A
                        && class == PointClassification::OnBoundary
                        && cooriented != Some(false)),
                false,
            ),
            BoolOperation::Intersection => (
                class == PointClassification::Inside
                    || (own_operand == Operand::A
                        && class == PointClassification::OnBoundary
                        && cooriented != Some(false)),
                false,
            ),
            BoolOperation::Subtraction => match own_operand {
                Operand::A => (
                    class == PointClassification::Outside
                        || (class == PointClassification::OnBoundary && cooriented != Some(true)),
                    false,
                ),
                Operand::B => (class == PointClassification::Inside, true),
            },
        };

        selections.push(FragmentSelection { keep, flip });
    }

    Ok(selections)
}

/// True when the point lies within `band` of any edge of the body —
/// sampled 16-station chordal distance per edge (exact for lines).
fn point_near_any_edge(brep: &BRep, p: &popkern_core::geometry::Point3, band: f64) -> bool {
    use popkern_core::geometry::curve::CurveEval as _;
    for (_, edge) in brep.edges.iter() {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let a = brep.vertices[he.origin].position;
        let b = brep
            .half_edges
            .get(he.twin)
            .map(|t| brep.vertices[t.origin].position)
            .unwrap_or(a);
        if let Some(c) = edge.curve.as_ref() {
            let mut prev = c.geometry.point_at(c.t_start);
            for i in 1..=16 {
                let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 16.0;
                let q = c.geometry.point_at(t);
                if seg_dist(&prev, &q, p) < band {
                    return true;
                }
                prev = q;
            }
        } else if seg_dist(&a, &b, p) < band {
            return true;
        }
    }
    false
}

fn seg_dist(
    a: &popkern_core::geometry::Point3,
    b: &popkern_core::geometry::Point3,
    p: &popkern_core::geometry::Point3,
) -> f64 {
    let ab = *b - *a;
    let l2 = ab.dot(&ab);
    if l2 < 1e-24 {
        return p.distance_to(a);
    }
    let t = ((*p - *a).dot(&ab) / l2).clamp(0.0, 1.0);
    (*a + ab * t).distance_to(p)
}

/// Chord midpoints of the loop's non-Shared half-edges, in loop
/// order, capped at `cap`. The first entry equals the tuned
/// edge-midpoint shortcut, so a unanimous vote reproduces legacy
/// sampling exactly.
/// Midpoint of a graph half-edge, ON THE EDGE'S CURVE when it has
/// one (chord midpoint otherwise). A curved edge's chord midpoint is
/// not a point of the fragment at all — a shallow rim-arc sliver cut
/// off by a breached wall has its arc CHORD lying exactly IN the
/// cutting plane, so the chord midpoint sampled the parity coin-flip
/// point on the other body's surface while the arc's belly sits
/// decisively outside it (measured: the sliver classified Inside and
/// a subtract kept a flap beyond the body). Straight edges are
/// unchanged (chord == curve).
fn he_curve_midpoint(graph: &IntersectionGraph, he_idx: usize) -> popkern_core::geometry::Point3 {
    use popkern_core::geometry::curve::CurveEval as _;
    let he = &graph.half_edges[he_idx];
    if let Some(c) = he
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[he.twin].curve.as_ref())
    {
        return c.geometry.point_at(0.5 * (c.t_start + c.t_end));
    }
    let a = graph.nodes[he.origin].position;
    let b = graph.nodes[graph.half_edges[he.twin].origin].position;
    popkern_core::geometry::Point3::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5)
}

fn non_shared_chord_midpoints(
    graph: &IntersectionGraph,
    loop_hes: &[usize],
    cap: usize,
) -> Vec<popkern_core::geometry::Point3> {
    use crate::merge::graph::EdgeSource;
    let mut out = Vec::new();
    for &he_idx in loop_hes {
        if out.len() >= cap {
            break;
        }
        let he = &graph.half_edges[he_idx];
        if matches!(he.source, EdgeSource::Shared { .. }) {
            continue;
        }
        out.push(he_curve_midpoint(graph, he_idx));
    }
    out
}

fn sample_fragment_point(
    graph: &IntersectionGraph,
    loop_hes: &[usize],
) -> popkern_core::geometry::Point3 {
    use crate::merge::graph::EdgeSource;

    // Prefer midpoint of a non-shared edge (away from intersection curve)
    for &he_idx in loop_hes {
        let he = &graph.half_edges[he_idx];
        match &he.source {
            EdgeSource::OriginalA(_)
            | EdgeSource::OriginalB(_)
            | EdgeSource::SplitA { .. }
            | EdgeSource::SplitB { .. } => {
                return he_curve_midpoint(graph, he_idx);
            }
            EdgeSource::Shared { .. } => continue,
        }
    }

    // Fallback: centroid of all loop vertices
    if loop_hes.is_empty() {
        return popkern_core::geometry::Point3::new(0.0, 0.0, 0.0);
    }
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sz = 0.0;
    for &he_idx in loop_hes {
        let p = graph.nodes[graph.half_edges[he_idx].origin].position;
        sx += p.x;
        sy += p.y;
        sz += p.z;
    }
    let n = loop_hes.len() as f64;
    popkern_core::geometry::Point3::new(sx / n, sy / n, sz / n)
}

/// Rotate a CLOSED intersection record's parameterization so its
/// closure point coincides with an edge-face hit lying on the curve.
/// Open records and closed records with no hit on them pass through
/// unchanged.
fn rotate_closed_record_to_hits(
    record: &crate::intersect::IntersectionRecord,
    hits: &[crate::intersect_edge_face::EdgeFaceHit],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> crate::intersect::IntersectionRecord {
    use popkern_core::geometry::curve::CurveEval as _;
    let c = &record.curve.curve_3d;
    let (t0, t1) = c.parameter_range();
    // Rotation is only meaningful for a genuinely CLOSED record. The
    // record's own flag is authoritative — the clip stage keeps it
    // honest (a multi-window clip emits open arcs flagged open). The
    // geometric gap check remains as a torn-closure guard, NOT as the
    // closure test itself: a pure distance gate at containment scale
    // (linear × 100) swallowed open corner arcs whose tip gap was
    // 0.17 and "rotated" them into gap-bridging garbage curves
    // (measured: cube∪sphere grew 12 arcs running through out-of-face
    // regions and the assembly went non-manifold).
    let closed = record.curve.is_closed
        && c.point_at(t0).distance_to(&c.point_at(t1)) < tolerance.linear() * 100.0;
    if !closed {
        return record.clone();
    }
    if std::env::var_os("PK_BOOL_TRACE").is_some() {
        let gap = c.point_at(t0).distance_to(&c.point_at(t1));
        eprintln!(
            "  [rotate] closed record considered for hit rotation: gap={gap:.3e} t[{t0:.3},{t1:.3}]"
        );
    }
    for hit in hits {
        // Nearest parameter on the curve to the hit point.
        let mut best = (f64::INFINITY, t0);
        for k in 0..=256 {
            let t = t0 + (t1 - t0) * k as f64 / 256.0;
            let d = c.point_at(t).distance_to(&hit.point);
            if d < best.0 {
                best = (d, t);
            }
        }
        let (mut lo, mut hi) = (
            (best.1 - (t1 - t0) / 256.0).max(t0),
            (best.1 + (t1 - t0) / 256.0).min(t1),
        );
        for _ in 0..48 {
            let m1 = lo + (hi - lo) / 3.0;
            let m2 = hi - (hi - lo) / 3.0;
            if c.point_at(m1).distance_to(&hit.point) < c.point_at(m2).distance_to(&hit.point) {
                hi = m2;
            } else {
                lo = m1;
            }
        }
        let t_hit = 0.5 * (lo + hi);
        let d_hit = c.point_at(t_hit).distance_to(&hit.point);
        if std::env::var_os("PK_BOOL_TRACE").is_some() {
            eprintln!(
                "  [rotate] hit ({:.3},{:.3},{:.3}) nearest t={t_hit:.4} d={d_hit:.3e}",
                hit.point.x, hit.point.y, hit.point.z
            );
        }
        if d_hit > tolerance.linear() * 100.0 {
            continue;
        }
        // Interior rotation only.
        if (t_hit - t0).abs() < (t1 - t0) * 1e-4 || (t1 - t_hit).abs() < (t1 - t0) * 1e-4 {
            continue;
        }
        let mut out = record.clone();
        let rotate = |curve: &popkern_core::geometry::NurbsCurve|
            -> Option<popkern_core::geometry::NurbsCurve> {
            let (a, b) = curve.split_at(t_hit).ok()?;
            // New curve = b then a, reparameterized to [0, 1] with
            // the join at the original closure (multiplicity-degree
            // interior knot — legal C0 join).
            let (b0, b1) = {
                let n = &b;
                (n.knots[0], *n.knots.last()?)
            };
            let (a0, a1) = {
                let n = &a;
                (n.knots[0], *n.knots.last()?)
            };
            let len_b = b1 - b0;
            let len_a = a1 - a0;
            let total = len_a + len_b;
            if total < 1e-12 {
                return None;
            }
            let mut knots: Vec<f64> = b
                .knots
                .iter()
                .map(|&k| (k - b0) / total)
                .collect();
            let join = len_b / total;
            // Skip a's leading clamp knots (degree+1 of them) — the
            // join point carries multiplicity = degree from b's tail
            // clamp minus the shared control point.
            let a_knots_tail: Vec<f64> = a
                .knots
                .iter()
                .skip(a.degree + 1)
                .map(|&k| join + (k - a0) / total)
                .collect();
            // Replace b's trailing clamp (degree+1 copies) with
            // multiplicity-degree at the join.
            for _ in 0..1 {
                knots.pop();
            }
            knots.extend(a_knots_tail);
            let mut control_points = b.control_points.clone();
            // JUNCTION CONTINUITY: the splice joins the original END
            // (b's tail) to the original START (a's head). In 3D
            // those coincide, but a UV curve on a periodic chart
            // unwraps continuously, so its two ends differ by one
            // period — splicing raw control points leaves a full-
            // period jump mid-curve and every downstream chart
            // embedding of the rotated record is garbage (measured:
            // the cross-drill pierce loop's drill-chart pcurve).
            // Translate the appended piece so the junction matches;
            // for 3D curves the delta is ~zero and this is a no-op.
            let delta = match (b.control_points.last(), a.control_points.first()) {
                (Some(bl), Some(af)) => *bl - *af,
                _ => popkern_core::geometry::Vec3::zero(),
            };
            control_points.extend(a.control_points.iter().skip(1).map(|&p| p + delta));
            let mut weights = b.weights.clone();
            weights.extend(a.weights.iter().skip(1).copied());
            Some(popkern_core::geometry::NurbsCurve {
                degree: b.degree,
                knots,
                control_points,
                weights,
            })
        };
        if let (Some(c3), Some(u1), Some(u2)) = (
            rotate(&record.curve.curve_3d),
            rotate(&record.curve.curve_uv1),
            rotate(&record.curve.curve_uv2),
        ) {
            out.curve.curve_3d = c3;
            out.curve.curve_uv1 = u1;
            out.curve.curve_uv2 = u2;
            // Weld the new closure exactly onto the hit point.
            if let Some(cp) = out.curve.curve_3d.control_points.first_mut() {
                *cp = hit.point;
            }
            if let Some(cp) = out.curve.curve_3d.control_points.last_mut() {
                *cp = hit.point;
            }
            return out;
        }
    }
    record.clone()
}
