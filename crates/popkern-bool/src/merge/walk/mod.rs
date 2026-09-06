//! Loop walking for face fragment extraction.
//!
//! Walks the intersection graph to extract closed loops of half-edges that
//! form face fragments. Each fragment consists of an outer loop (CCW) and
//! zero or more inner loops (CW holes).

mod boundary_predicates;
pub(crate) mod chart_boundary_substitute;
pub(crate) mod chart_branches;
pub(crate) mod chart_cell_fold;
pub(crate) mod chart_duplicate_cover;
pub(crate) mod chart_embedding;
pub(crate) mod chart_poles;
pub(crate) mod chart_subdivision;
pub(crate) mod chart_walk;
#[cfg(test)]
mod chart_walk_tests;
mod hole_nesting;
mod next_edge;
pub(crate) mod walk_geometry;
mod walk_loop_measures;
mod walk_trace;

use std::collections::HashSet;

use popkern_core::topology::BRep;

use super::graph::{FaceRef, IntersectionGraph};
use boundary_predicates::node_on_face_boundary;
use next_edge::next_around_face;
use walk_geometry::signed_area;

/// A face fragment — loops forming a subset of an original face.
#[derive(Debug)]
pub(crate) struct Fragment {
    /// The face this fragment belongs to.
    pub face: FaceRef,
    /// Half-edge indices forming the outer loop (CCW).
    pub outer_loop: Vec<usize>,
    /// Half-edge indices forming inner loops (CW holes).
    pub inner_loops: Vec<Vec<usize>>,
}

/// Walk all face fragments in the graph (graph-merge pipeline entry).
///
/// Chart-primary like the cellular pipeline: each face whose substrate
/// passes I1–I5 walks as a chart permutation; the rest take the 3D-Newell
/// walk per face. Graph-merge graphs lack the cellular phase's seam
/// splits, so its periodic faces fall back more often — the per-face
/// fallback absorbs exactly those (PK-108 retirement: the 3D-Newell walk
/// is fallback-only on every pipeline).
pub(crate) fn walk_all_fragments(
    graph: &IntersectionGraph,
    brep_a: &BRep,
    brep_b: &BRep,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Vec<Fragment> {
    walk_all_fragments_ext(
        graph,
        brep_a,
        brep_b,
        false,
        true,
        &[],
        &mut Vec::new(),
        tolerance,
    )
}

/// `coincident_interface`: the boolean carries coplanar face
/// overlaps — enables the disjoint-loop election override, whose
/// sampled-polygon containment test is currently only trusted on
/// that configuration (it mis-voted a through-slot ring on a long
/// tapered face — see the flush-gusset pin).
///
/// `chart_primary`: walk each face as a chart-subdivision permutation
/// when its substrate passes I1–I5, with the 3D-Newell walk as recorded
/// per-face fallback (PK-104 stage 3; the cellular pipeline sets this,
/// graph-merge keeps the 3D-Newell walk — its graphs lack seam splits).
///
/// `overlap_pairs`: the coplanar face-overlap pairs from coincidence
/// detection — the substrate's partial-cover dedup is scoped to
/// records authored by exactly these pairs (PK-106).
#[allow(clippy::too_many_arguments)]
pub(crate) fn walk_all_fragments_ext(
    graph: &IntersectionGraph,
    brep_a: &BRep,
    brep_b: &BRep,
    coincident_interface: bool,
    chart_primary: bool,
    overlap_pairs: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::topology::FaceHandle,
    )],
    walk_fallbacks: &mut Vec<(FaceRef, String)>,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Vec<Fragment> {
    let mut fragments = Vec::new();

    // Opposite-operand face AABBs (PK-438 exterior-disjoint promotion).
    let opposite_aabbs = hole_nesting::OppositeAabbs::new(brep_a, brep_b);

    for (face_ref, he_indices) in &graph.face_edges {
        walk_trace::dump_face_pre(graph, face_ref, he_indices);
        let face_hes: HashSet<usize> = he_indices.iter().copied().collect();
        let mut visited: HashSet<usize> = HashSet::new();
        let mut loops: Vec<Vec<usize>> = Vec::new();

        // PK-104 stage 3: chart-subdivision permutation walk, primary
        // for the cellular pipeline. Faces whose substrate fails the
        // I1–I5 checks fall through to the 3D-Newell walk below, with
        // the reason observable under the boolean trace topic and recorded in
        // `walk_fallbacks`.
        let mut walked_by_chart = false;
        // Chart polylines + σ + chart u-period for the role assembly
        // below (present only when the chart walk produced the orbit
        // partition).
        let mut chart_roles: Option<(Vec<chart_walk::ChartPolyline>, f64, Option<f64>)> = None;
        // Since PK-107 every face of a coincident boolean
        // chart-walks, overlap faces included: the substrate's
        // duplicate-cover dedup, gap-substitute completion, and
        // interior-graze filter give the overlap faces the same
        // fragment population the tuned arbitration expects (their
        // chart-emitted web skins classify OnBoundary with correct
        // coorientation).
        if chart_primary
            && let Some((surface, orientation)) =
                chart_subdivision::face_surface(face_ref, brep_a, brep_b)
        {
            let sub = chart_subdivision::build_face_subdivision(
                graph,
                face_ref,
                he_indices,
                surface,
                orientation,
                overlap_pairs,
                match face_ref.operand {
                    super::graph::Operand::A => brep_a,
                    super::graph::Operand::B => brep_b,
                },
                tolerance,
            );
            if sub.complete() {
                if let Some(orbits) = chart_walk::walk_chart_orbits(&sub, graph) {
                    // Chart-role assembly is scoped to faces whose
                    // stored boundary is a SEAM SPINE (every
                    // non-record half-edge's twin registered to this
                    // same face): there the 3D-Newell election
                    // degenerates BY CONSTRUCTION — two anti-parallel
                    // border copies measure |N| ≈ 0 (a full sphere's
                    // seam orbit lost the outer election to a record
                    // triangle). Ordinary faces keep the 3D-Newell
                    // election: its size/winding rules carry tuned
                    // tool-side semantics (a plane-splitter's cutting
                    // face serves BOTH directions of an interior ring
                    // as caps) that the pure role model does not
                    // reproduce (measured: unscoped roles broke four
                    // holder/rail gates).
                    let boundary_hes: Vec<usize> = he_indices
                        .iter()
                        .copied()
                        .filter(|&h| {
                            !matches!(
                                graph.half_edges[h].source,
                                super::graph::EdgeSource::Shared { .. }
                            )
                        })
                        .collect();
                    let seam_spine = !boundary_hes.is_empty()
                        && boundary_hes
                            .iter()
                            .all(|&h| face_hes.contains(&graph.half_edges[h].twin));
                    // WRAPPING-WALL extension: a closed periodic wall
                    // (rims + seam) degenerates the 3D Newell
                    // election the same way — its boundary wraps the
                    // axis, |N| ≈ 0 at wall scale, hole nesting
                    // projects onto a garbage plane, and pierce rims
                    // never nest (measured: the cross-drilled boss's
                    // wall kept 0 inner loops and sampled inside its
                    // own invisible hole). Any face whose substrate
                    // carries a seam-like embedding is such a wall;
                    // planar splitter faces (the tuned 3D-Newell
                    // semantics the scoping protects) have no seam
                    // and keep the 3D-Newell election.
                    let wrapping_wall = sub.embeddings.values().any(|e| e.seam_like);
                    // ISLAND-PAIR extension (PK-341): a face carrying
                    // BOTH directions of a closed interior ring
                    // (island splits — mutual-twin orbit pairs) routes
                    // to chart roles too: the 3D-Newell election called
                    // both copies "separate" on tool faces and its
                    // nesting rejected the concentric inner rim on
                    // winding (measured: the tube's groove walls lost
                    // their holes and the fallback served a non-closed
                    // body validate-clean). Planar splitter faces (the
                    // tuned 3D-Newell semantics the scoping protects)
                    // carry no twin-paired rings.
                    let island_pair = orbits.len() >= 2
                        && orbits.iter().enumerate().any(|(i, o)| {
                            !o.is_empty()
                                && o.iter().all(|&h| {
                                    let tw = graph.half_edges[h].twin;
                                    orbits
                                        .iter()
                                        .enumerate()
                                        .any(|(j, o2)| j != i && o2.contains(&tw))
                                })
                        });
                    if seam_spine || wrapping_wall || island_pair {
                        chart_roles = Some((
                            orbits
                                .iter()
                                .map(|o| chart_walk::orbit_polyline(&sub, o))
                                .collect(),
                            sub.sigma,
                            sub.u_period,
                        ));
                    }
                    loops = orbits;
                    walked_by_chart = true;
                    if popkern_core::trace::trace_enabled("boolean") {
                        eprintln!(
                            "  [walk] {face_ref:?}: chart orbits {} (sigma {}) filtered {:?}",
                            loops.len(),
                            sub.sigma,
                            sub.report.filtered
                        );
                    }
                } else {
                    walk_fallbacks.push((*face_ref, "chart orbit collision".to_string()));
                    if popkern_core::trace::trace_enabled("boolean") {
                        eprintln!(
                            "  [walk] {face_ref:?}: chart orbit collision — 3D-Newell fallback"
                        );
                    }
                }
            } else {
                let issues: Vec<String> = sub.report.issues.iter().map(|i| i.to_string()).collect();
                walk_fallbacks.push((*face_ref, issues.join("; ")));
                if popkern_core::trace::trace_enabled("boolean") {
                    eprintln!(
                        "  [walk] {face_ref:?}: chart substrate incomplete — 3D-Newell fallback ({})",
                        issues.join("; ")
                    );
                }
            }
        }

        // A boundary edge covered by a record is not walked (see
        // `boundary_predicates::boundary_covered_by_records`).
        if !walked_by_chart {
            boundary_predicates::skip_covered(graph, he_indices, tolerance, &mut visited);
        }
        // Start seam-like first (twin on this face: the outer "spine" that
        // connects self-loop circles), then other non-shared, then shared.
        let mut sorted_hes: Vec<usize> = he_indices.clone();
        sorted_hes.sort_by_key(|&idx| {
            let is_shared = matches!(
                graph.half_edges[idx].source,
                super::graph::EdgeSource::Shared { .. }
            );
            // Primary key: seam-like (0) → non-shared (1) → shared (2).
            // Secondary key: half-edge index (lower first) for deterministic
            // fragment ordering — ensures the walk starts from the same edge
            // regardless of face_edges insertion order.
            let priority = if is_shared {
                2u8
            } else {
                let twin = graph.half_edges[idx].twin;
                let twin_on_same_face = graph.half_edges[twin].face.as_ref() == Some(face_ref);
                if twin_on_same_face { 0 } else { 1 }
            };
            (priority, idx)
        });

        for &start_he in &sorted_hes {
            if walked_by_chart {
                break;
            }
            if visited.contains(&start_he) {
                continue;
            }
            // Skip cross-registered self-loop half-edges as starting
            // points. They appear in face_edges due to cross-registration
            // but should not form standalone loops on this face.
            let twin_he = graph.half_edges[start_he].twin;
            let is_self_loop =
                graph.half_edges[start_he].origin == graph.half_edges[twin_he].origin;
            if is_self_loop && graph.half_edges[start_he].face.as_ref() != Some(face_ref) {
                continue;
            }

            let mut current_loop = Vec::new();
            let mut he = start_he;

            loop {
                if visited.contains(&he) {
                    break;
                }
                visited.insert(he);

                current_loop.push(he);

                // Find next half-edge around this face.
                match next_around_face(
                    graph,
                    he,
                    face_ref,
                    &face_hes,
                    brep_a,
                    brep_b,
                    Some(start_he),
                    &visited,
                    &current_loop,
                ) {
                    Some(next_he) => {
                        if popkern_core::trace::trace_enabled("fragment-walk") {
                            eprintln!("      [step] {he} -> {next_he}");
                        }
                        he = next_he;
                    }
                    None => break,
                }

                if he == start_he {
                    break;
                }
            }

            // Accept loops of length 1 only for self-loop edges (closed
            // curves like circles), otherwise require at least 2 edges.
            let min_len = if current_loop.len() == 1 {
                let he = current_loop[0];
                let twin = graph.half_edges[he].twin;
                let start = graph.half_edges[he].origin;
                let end = graph.half_edges[twin].origin;
                if start == end { 1 } else { 2 }
            } else {
                2
            };
            // TOPOLOGICAL CLOSURE: a chain is a loop only when its
            // last half-edge ends at the seed's origin. A walk that
            // dead-ends (no onward candidate at a junction) used to
            // be pushed as a "loop" anyway — leftover cross-
            // registered directions at a breached-wall mouth formed
            // a 2-half-edge pseudo-loop (mouth arc + one vertical,
            // 27 mm from closing) that assembled into the result and
            // failed loop-closure validation. Dead-ended chains
            // bound nothing; their half-edges stay consumed so
            // seeding cannot livelock.
            let closes = current_loop.last().is_some_and(|&last| {
                let last_twin = graph.half_edges[last].twin;
                let end_node = graph.half_edges[last_twin].origin;
                let start_node = graph.half_edges[current_loop[0]].origin;
                // Half-split trims legitimately dead-end ON a face
                // boundary edge (the boundary was not split there in
                // this phase; downstream completes the cut — the
                // sheet splitter depends on these). Interior
                // dead-ends bound nothing.
                end_node == start_node || {
                    let last_he = current_loop[current_loop.len() - 1];
                    // Only a STRUCTURAL dead-end closes as a half-split
                    // chain: no same-face departure exists at the end
                    // node at all — the boundary genuinely was not
                    // split there. A dead-end caused by CONSUMPTION
                    // (departures exist but earlier loops claimed
                    // them) is leftover residue: accepting one wired a
                    // two-piece coincident-rim chain as an "inner
                    // loop" with a rim-chord gap and the served body
                    // failed loop closure (PK-439, flush-cap
                    // prism ∪ hub).
                    let structural = !graph.nodes[end_node]
                        .outgoing
                        .iter()
                        .any(|&o| o != last_twin && o != last_he && face_hes.contains(&o));
                    structural
                        && node_on_face_boundary(graph, end_node, &face_hes, last_he, last_twin)
                }
            });
            if current_loop.len() >= min_len && closes {
                loops.push(current_loop);
            } else if !closes && popkern_core::trace::trace_enabled("fragment-walk") {
                eprintln!(
                    "      [drop-open-chain] {} hes from he {}",
                    current_loop.len(),
                    current_loop.first().copied().unwrap_or(usize::MAX)
                );
            }
        }

        if loops.is_empty() {
            continue;
        }

        // Chart-role assembly (PK-111): orbit roles are decided by
        // σ-folded chart shoelace sign — by construction of the
        // material-left walk, not by 3D election (see
        // `assemble_chart_roles`). Falls through to the 3D-Newell
        // election only when no orbit winds positive (substrate lie).
        if let Some((polys, sigma, u_period)) = chart_roles.take()
            && let Some(asm) = chart_walk::assemble_chart_roles(
                graph,
                &loops,
                &polys,
                sigma,
                u_period,
                &chart_walk::source_inner_edge_set(face_ref, brep_a, brep_b),
            )
        {
            let dump = popkern_core::trace::trace_enabled("loop-dump-all");
            let trace = popkern_core::trace::trace_enabled("boolean");
            if trace {
                let (oi, holes) = &asm.regions[0];
                eprintln!(
                    "  [walk] {face_ref:?}: chart roles — outer {} hes, {} inner {:?}, {} separate",
                    loops[*oi].len(),
                    holes.len(),
                    holes.iter().map(|&h| loops[h].len()).collect::<Vec<_>>(),
                    asm.regions.len() - 1
                );
            }
            for (ri, (oi, holes)) in asm.regions.iter().enumerate() {
                if dump {
                    for (tag, lp) in
                        std::iter::once((if ri == 0 { "outer" } else { "separate" }, &loops[*oi]))
                            .chain(holes.iter().map(|&h| ("inner", &loops[h])))
                    {
                        let hes: Vec<String> = lp
                            .iter()
                            .map(|&h| {
                                let n = graph.half_edges[h].origin;
                                let e = graph.half_edges[graph.half_edges[h].twin].origin;
                                let p = graph.nodes[n].position;
                                format!("{h}:n{n}->n{e}({:.2},{:.2},{:.2})", p.x, p.y, p.z)
                            })
                            .collect();
                        eprintln!("  [loopdump-all] {face_ref:?} {tag} [{}]", hes.join(" "));
                    }
                }
                fragments.push(Fragment {
                    face: *face_ref,
                    outer_loop: loops[*oi].clone(),
                    inner_loops: holes.iter().map(|&h| loops[h].clone()).collect(),
                });
            }
            continue;
        }

        // Classify loops by signed area.  The loop with the largest
        // absolute area is the outer boundary.  Smaller loops (less than
        // half the outer area) are inner loops (holes) that belong to
        // the same fragment.  Large loops form separate fragments (they
        // arise when an intersection curve splits a face).
        let areas: Vec<f64> = loops.iter().map(|lp| signed_area(graph, lp)).collect();
        let max_area = areas
            .iter()
            .map(|a| a.abs())
            .fold(f64::NEG_INFINITY, f64::max);
        let threshold = max_area * 0.5;

        // Outer = largest |area|. signed_area samples curve geometry
        // for closed/self-loop edges, so a 1-half-edge closed rim
        // (cylinder cap) reports its true disk area. The historical
        // fallback — "longest loop wins when self-loops are present"
        // — chose a 4-edge interior pocket loop over the 1-edge rim
        // and every pocket cut into a closed-rim cap lost its hole
        // (measured: cap kept with inner 0, pocket walls orphaned,
        // 307 open mesh edges). Length only breaks area ties now.
        let outer_idx = {
            let mut best = 0usize;
            for i in 1..loops.len() {
                let ai = areas[i].abs();
                let ab = areas[best].abs();
                if ai > ab * (1.0 + 1e-9)
                    || ((ai - ab).abs() <= ab * 1e-9 && loops[i].len() > loops[best].len())
                {
                    best = i;
                }
            }
            best
        };

        // Extract the outer loop (swap-remove to avoid shifting).
        let outer_loop = loops.swap_remove(outer_idx);
        // Adjust areas to match new indices after swap_remove.
        let mut remaining_areas = areas;
        remaining_areas.swap_remove(outer_idx);

        let mut inner_loops_out: Vec<Vec<usize>> = Vec::new();
        let mut separate_outers = Vec::new();

        // Collect the half-edge indices in the outer loop for
        // twin-overlap detection.
        let outer_he_set: HashSet<usize> = outer_loop.iter().copied().collect();

        let outer_newell = walk_geometry::loop_newell(graph, &outer_loop);
        // The winding override below applies to WORKPIECE (operand A)
        // faces only. On the TOOL side a large interior ring's both
        // directions serve as candidate caps for the result pieces (a
        // plane-split's cutting box face contributes the cap to each
        // half), and pre-assigning one as the tool-frame's hole
        // starved a result piece of its cap (measured on the
        // solid-by-plane splitter).
        let face_on_workpiece = matches!(face_ref.operand, super::graph::Operand::A);
        for (i, lp) in loops.into_iter().enumerate() {
            // A loop whose half-edges share an underlying graph edge
            // with the outer loop is a face-split boundary, not a hole.
            // Detect this by checking if any half-edge in the loop has
            // its twin in the outer loop.
            let shares_edge_with_outer = lp.iter().any(|&he_idx| {
                let twin = graph.half_edges[he_idx].twin;
                outer_he_set.contains(&twin)
            });
            // Winding RESCUE for large PAIRED holes: the 0.5·area
            // threshold alone called a 76 mm seat ring on a small
            // block face a "separate fragment" purely because the
            // hole was big relative to the face. When BOTH directions
            // of an interior ring were walked, the one winding
            // OPPOSITE the outer (Newell dot < 0) is the hole and its
            // twin is the interior piece — regardless of size. An
            // UNPAIRED big loop keeps the size rule (a sheet-split's
            // interior piece walks once and must stay a piece), and
            // promoting winding to the primary test for small loops
            // was measured to break six tuned suites.
            // Winding override for LARGE holes: the 0.5·area
            // threshold alone called a 76 mm seat ring on a small
            // block face a "separate fragment" purely because the
            // hole was big relative to the face. A big loop winding
            // clearly OPPOSITE the outer (Newell dot < 0) is its
            // hole. Small loops keep the size rule — demoting small
            // same-winding loops was measured to break six suites.
            let separate_by_shape = if remaining_areas[i].abs() >= threshold {
                if face_on_workpiece {
                    let winding_dot = walk_geometry::loop_newell(graph, &lp).dot(&outer_newell);
                    winding_dot >= -1e-9
                } else {
                    true
                }
            } else {
                // CONTAINMENT override: a hole must lie INSIDE the
                // outer loop's region. The size proxy alone called a
                // disjoint root-collar loop a "hole" of a swoop
                // fragment on one mirror side and not the other
                // (area 263 vs threshold 352 — pure size accident).
                // A loop whose nodes are strictly outside the outer
                // polygon can never be its hole; loops inside keep
                // the existing rules. (Disjoint mis-filed holes are
                // ALSO a load-bearing waypoint: hole nesting relocates
                // them to their true regions, and the ones nesting
                // cannot place are PROMOTED after it — PK-343. Do not
                // evict them here.)
                coincident_interface
                    && walk_geometry::loop_outside_outer(graph, &lp, &outer_loop, &outer_newell)
            };
            if popkern_core::trace::trace_enabled("fragment-walk") {
                eprintln!(
                    "    [elect] loop {:?}{} area {:.1} thr {:.1} shares={} sep={}",
                    &lp[..lp.len().min(6)],
                    if lp.len() > 6 { "…" } else { "" },
                    remaining_areas[i],
                    threshold,
                    shares_edge_with_outer,
                    separate_by_shape
                );
            }
            if shares_edge_with_outer || separate_by_shape {
                separate_outers.push(lp);
            } else {
                // Both directions of an interior closed ring reach here.
                // They are NOT duplicates: one (winding OPPOSITE the
                // outer loop) is the hole boundary of this fragment;
                // the other (winding WITH the outer) is the OUTER
                // boundary of the interior piece — a genuine separate
                // fragment. Discarding it as a "duplicate" deleted the
                // interior piece wholesale (measured: a rotated slab's
                // cutting face lost the quad inside the prism, so the
                // subtracted wedge had no rest plane — shell gap = the
                // full prism width).
                let twin_of = |target: &[usize]| {
                    inner_loops_out.iter().position(|existing| {
                        existing.len() == target.len()
                            && target.iter().all(|&he_idx| {
                                let twin = graph.half_edges[he_idx].twin;
                                existing.contains(&twin)
                            })
                    })
                };
                match twin_of(&lp) {
                    None => inner_loops_out.push(lp),
                    // When BOTH directions of the ring were walked,
                    // assign roles by 3D winding: a hole winds
                    // OPPOSITE its outer loop (Newell dot < 0); the
                    // same-winding direction is the interior piece's
                    // outer boundary. First-walked order is chart-
                    // iteration dependent — with honest cap charts it
                    // flipped on boolean-authored faces, the hole
                    // claimed the SAME graph half-edges the other
                    // side's fragments consume, and the assembly
                    // minted duplicate ring edges (measured: a
                    // pre-filleted lip union left both ring copies
                    // single-use, 8 degenerate junction stubs, and
                    // fell back to graph-merge). Winding is decided
                    // against the outer as walked, which is exactly
                    // the direction assembly wires. Indecisive
                    // windings keep first-walked order.
                    Some(prior_idx) => {
                        let dot_new = walk_geometry::loop_newell(graph, &lp).dot(&outer_newell);
                        let dot_prior =
                            walk_geometry::loop_newell(graph, &inner_loops_out[prior_idx])
                                .dot(&outer_newell);
                        if dot_new < 0.0 && dot_prior > 0.0 {
                            let piece = std::mem::replace(&mut inner_loops_out[prior_idx], lp);
                            separate_outers.push(piece);
                        } else {
                            separate_outers.push(lp);
                        }
                    }
                }
            }
        }

        // the loop-dump-all trace topic: dump every fragment's loops regardless of
        // size (PK-107 forensics: comparing 11- vs 13-half-edge outer
        // loops between the 3D-Newell and chart walks).
        walk_trace::dump_election_loops(
            graph,
            face_ref,
            &outer_loop,
            &inner_loops_out,
            &separate_outers,
        );
        let separate_inners = hole_nesting::assign_holes_to_regions(
            graph,
            &outer_newell,
            &mut inner_loops_out,
            &separate_outers,
        );
        let promoted = hole_nesting::promote_disjoint_turnaround_holes(
            graph,
            &mut inner_loops_out,
            &outer_loop,
        );
        if popkern_core::trace::trace_enabled("boolean") {
            eprintln!(
                "  [walk] {face_ref:?}: outer {} hes, {} inner {:?}, {} separate",
                outer_loop.len(),
                inner_loops_out.len(),
                inner_loops_out.iter().map(|l| l.len()).collect::<Vec<_>>(),
                separate_outers.len()
            );
        }
        // PK-438: a disjoint below-body tool protrusion glued in as an
        // inner loop is promoted to its own fragment so classify culls
        // it independently (else it strands the shared flush rim).
        let exterior_promoted = hole_nesting::promote_exterior_disjoint_holes(
            graph,
            &mut inner_loops_out,
            &outer_loop,
            face_ref.operand,
            &opposite_aabbs,
            tolerance,
        );
        for lp in promoted.into_iter().chain(exterior_promoted) {
            fragments.push(Fragment {
                face: *face_ref,
                outer_loop: lp,
                inner_loops: Vec::new(),
            });
        }
        fragments.push(Fragment {
            face: *face_ref,
            outer_loop,
            inner_loops: inner_loops_out,
        });

        for (lp, inners) in separate_outers.into_iter().zip(separate_inners) {
            fragments.push(Fragment {
                face: *face_ref,
                outer_loop: lp,
                inner_loops: inners,
            });
        }
    }

    fragments
}

#[cfg(test)]
mod tests;
