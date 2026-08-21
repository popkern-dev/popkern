//! Loop walking for face fragment extraction.
//!
//! Walks the intersection graph to extract closed loops of half-edges that
//! form face fragments. Each fragment consists of an outer loop (CCW) and
//! zero or more inner loops (CW holes).

mod boundary_predicates;
pub(crate) mod chart_boundary_substitute;
pub(crate) mod chart_branches;
pub(crate) mod chart_cell_fold;
pub(crate) mod chart_embedding;
pub(crate) mod chart_poles;
pub(crate) mod chart_subdivision;
pub(crate) mod chart_walk;
#[cfg(test)]
mod chart_walk_tests;
mod hole_nesting;
pub(crate) mod walk_geometry;

use std::collections::HashSet;

use popkern_core::topology::BRep;

use super::graph::{FaceRef, IntersectionGraph};
use boundary_predicates::{is_self_loop_shared, node_on_face_boundary};
use walk_geometry::{pick_by_angle, signed_area};

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
/// passes I1–I5 walks as a chart permutation; the rest take the legacy
/// walk per face. Graph-merge graphs lack the cellular phase's seam
/// splits, so its periodic faces fall back more often — the per-face
/// fallback absorbs exactly those (PK-108 retirement: the legacy walk
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
/// when its substrate passes I1–I5, with the legacy walk as recorded
/// per-face fallback (PK-104 stage 3; the cellular pipeline sets this,
/// graph-merge keeps the legacy walk — its graphs lack seam splits).
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

    for (face_ref, he_indices) in &graph.face_edges {
        if std::env::var_os("PK_WALK_TRACE").is_some() {
            let n_int = he_indices
                .iter()
                .filter(|&&h| {
                    matches!(
                        graph.half_edges[h].source,
                        crate::merge::graph::EdgeSource::Shared { .. }
                    )
                })
                .count();
            eprintln!(
                "  [walk-pre] {face_ref:?}: {} hes ({n_int} intersection)",
                he_indices.len()
            );
            if std::env::var_os("PK_WALK_TRACE").is_some() {
                for &h in he_indices.iter() {
                    let he = &graph.half_edges[h];
                    let o = he.origin;
                    let e = graph.half_edges[he.twin].origin;
                    let po = graph.nodes[o].position;
                    let pe = graph.nodes[e].position;
                    eprintln!(
                        "      [he] {h} n{o}({:.9},{:.9},{:.9}) -> n{e}({:.9},{:.9},{:.9}) twin {} src {:?} curve {}",
                        po.x,
                        po.y,
                        po.z,
                        pe.x,
                        pe.y,
                        pe.z,
                        he.twin,
                        he.source,
                        he.curve
                            .as_ref()
                            .or(graph.half_edges[he.twin].curve.as_ref())
                            .map(|c| format!("[{:.2},{:.2}]", c.t_start, c.t_end))
                            .unwrap_or_else(|| "none".into())
                    );
                }
            }
        }
        let face_hes: HashSet<usize> = he_indices.iter().copied().collect();
        let mut visited: HashSet<usize> = HashSet::new();
        let mut loops: Vec<Vec<usize>> = Vec::new();

        // PK-104 stage 3: chart-subdivision permutation walk, primary
        // for the cellular pipeline. Faces whose substrate fails the
        // I1–I5 checks fall through to the legacy walk below, with
        // the reason observable under PK_BOOL_TRACE and recorded in
        // `walk_fallbacks`.
        let mut walked_by_chart = false;
        // Chart polylines + σ for the role assembly below (present
        // only when the chart walk produced the orbit partition).
        let mut chart_roles: Option<(Vec<chart_walk::ChartPolyline>, f64)> = None;
        // A/B lever for divergence adjudication: PK_WALK_CHART=0
        // forces the legacy walk everywhere. Since PK-107 every face
        // of a coincident boolean chart-walks, overlap faces
        // included: the substrate's duplicate-cover dedup,
        // gap-substitute completion, and interior-graze filter give
        // the overlap faces the same fragment population the tuned
        // arbitration expects (their chart-emitted web skins classify
        // OnBoundary with correct coorientation).
        let chart_disabled = std::env::var("PK_WALK_CHART").as_deref() == Ok("0");
        if chart_primary
            && !chart_disabled
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
                    // same face): there the legacy 3D Newell election
                    // degenerates BY CONSTRUCTION — two anti-parallel
                    // border copies measure |N| ≈ 0 (a full sphere's
                    // seam orbit lost the outer election to a record
                    // triangle). Ordinary faces keep the legacy
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
                    // planar splitter faces (the tuned legacy
                    // semantics the scoping protects) have no seam
                    // and keep the legacy election.
                    let wrapping_wall = sub.embeddings.values().any(|e| e.seam_like);
                    if seam_spine || wrapping_wall {
                        chart_roles = Some((
                            orbits
                                .iter()
                                .map(|o| chart_walk::orbit_polyline(&sub, o))
                                .collect(),
                            sub.sigma,
                        ));
                    }
                    loops = orbits;
                    walked_by_chart = true;
                    if std::env::var_os("PK_BOOL_TRACE").is_some() {
                        eprintln!(
                            "  [walk] {face_ref:?}: chart orbits {} (sigma {}) filtered {:?}",
                            loops.len(),
                            sub.sigma,
                            sub.report.filtered
                        );
                    }
                } else {
                    walk_fallbacks.push((*face_ref, "chart orbit collision".to_string()));
                    if std::env::var_os("PK_BOOL_TRACE").is_some() {
                        eprintln!("  [walk] {face_ref:?}: chart orbit collision — legacy fallback");
                    }
                }
            } else {
                let issues: Vec<String> = sub.report.issues.iter().map(|i| i.to_string()).collect();
                walk_fallbacks.push((*face_ref, issues.join("; ")));
                if std::env::var_os("PK_BOOL_TRACE").is_some() {
                    eprintln!(
                        "  [walk] {face_ref:?}: chart substrate incomplete — legacy fallback ({})",
                        issues.join("; ")
                    );
                }
            }
        }

        // Sort starting edges: seam-like first, then other non-shared, then shared.
        // Seam-like edges (where the twin is also on this face) form the outer
        // boundary "spine" that connects self-loop circles. Starting from them
        // ensures the outer boundary walk includes expanded cap circles rather
        // than consuming them as standalone loops.
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
                        if std::env::var_os("PK_WALK_TRACE").is_some() {
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
                    node_on_face_boundary(graph, end_node, &face_hes, last_he, last_twin)
                }
            });
            if current_loop.len() >= min_len && closes {
                loops.push(current_loop);
            } else if !closes && std::env::var_os("PK_WALK_TRACE").is_some() {
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
        // `assemble_chart_roles`). Falls through to the legacy
        // election only when no orbit winds positive (substrate lie).
        if let Some((polys, sigma)) = chart_roles.take()
            && let Some(asm) = chart_walk::assemble_chart_roles(graph, &loops, &polys, sigma)
        {
            let dump = std::env::var_os("PK_LOOP_DUMP_ALL").is_some();
            let trace = std::env::var_os("PK_BOOL_TRACE").is_some();
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
                // the existing rules.
                coincident_interface
                    && walk_geometry::loop_outside_outer(graph, &lp, &outer_loop, &outer_newell)
            };
            if std::env::var_os("PK_WALK_TRACE").is_some() {
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

        // PK_LOOP_DUMP_ALL: dump every fragment's loops regardless of
        // size (PK-107 forensics: comparing 11- vs 13-half-edge outer
        // loops between the legacy and chart walks).
        if std::env::var_os("PK_LOOP_DUMP_ALL").is_some() {
            for (tag, lp) in std::iter::once(("outer", &outer_loop))
                .chain(inner_loops_out.iter().map(|l| ("inner", l)))
                .chain(separate_outers.iter().map(|l| ("separate", l)))
            {
                let hes: Vec<String> = lp
                    .iter()
                    .map(|&h| {
                        let n = graph.half_edges[h].origin;
                        let p = graph.nodes[n].position;
                        let shared = matches!(
                            graph.half_edges[h].source,
                            super::graph::EdgeSource::Shared { .. }
                        );
                        format!(
                            "{h}@n{n}({:.1},{:.1},{:.1}){}",
                            p.x,
                            p.y,
                            p.z,
                            if shared { "S" } else { "" }
                        )
                    })
                    .collect();
                eprintln!("  [loopdump-all] {face_ref:?} {tag} [{}]", hes.join(" "));
            }
        }
        if std::env::var_os("PK_LOOP_DUMP").is_some() && outer_loop.len() >= 20 {
            eprintln!("  [loopdump] {face_ref:?} outer {} hes:", outer_loop.len());
            for &h in &outer_loop {
                let n = graph.half_edges[h].origin;
                let p = graph.nodes[n].position;
                let shared = matches!(
                    graph.half_edges[h].source,
                    super::graph::EdgeSource::Shared { .. }
                );
                eprintln!(
                    "    he {h} ({:.1},{:.1},{:.1}) shared={shared}",
                    p.x, p.y, p.z
                );
            }
        }
        let separate_inners = hole_nesting::assign_holes_to_regions(
            graph,
            &outer_newell,
            &mut inner_loops_out,
            &separate_outers,
        );
        if std::env::var_os("PK_BOOL_TRACE").is_some() {
            eprintln!(
                "  [walk] {face_ref:?}: outer {} hes, {} inner {:?}, {} separate",
                outer_loop.len(),
                inner_loops_out.len(),
                inner_loops_out.iter().map(|l| l.len()).collect::<Vec<_>>(),
                separate_outers.len()
            );
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

/// Find the next half-edge around a face after `he_idx`.
///
/// At the end node of `he_idx`, finds the outgoing half-edge belonging to
/// the same face that makes the smallest counter-clockwise angle from the
/// incoming direction.  This uses per-face angular comparison rather than
/// the global radial sort, because intersection nodes sit at the junction
/// of multiple faces with different normals.
///
/// **Cut behavior:** when `he_idx` is a closed self-loop on a `Shared`
/// edge — i.e. an intersection circle that the boolean inserted as a
/// cut through the face — and the previous edge in `current_loop`
/// has its twin among the candidates, the walker bounces back along
/// that twin. This treats the intersection circle as a topological
/// cut: walk in via some boundary edge, traverse the circle, walk
/// back out the same way. Without this rule, the walker keeps going
/// past the cut and merges multiple regions into a single mega-loop
/// (the cube-minus-cylinder Euler-violation bug).
#[allow(clippy::too_many_arguments)]
fn next_around_face(
    graph: &IntersectionGraph,
    he_idx: usize,
    face: &FaceRef,
    face_hes: &HashSet<usize>,
    brep_a: &BRep,
    brep_b: &BRep,
    start_he: Option<usize>,
    visited: &HashSet<usize>,
    current_loop: &[usize],
) -> Option<usize> {
    let twin_idx = graph.half_edges[he_idx].twin;
    let end_node = graph.half_edges[twin_idx].origin;
    let outgoing = &graph.nodes[end_node].outgoing;

    // For self-loop edges (origin == twin.origin), exclude both the current
    // edge AND its twin to avoid re-traversing the same self-loop or U-turning.
    // For normal edges, only exclude the twin (standard U-turn prevention).
    let is_self_loop = graph.half_edges[he_idx].origin == end_node;

    let candidates: Vec<usize> = outgoing
        .iter()
        .copied()
        .filter(|&idx| {
            if !face_hes.contains(&idx) {
                return false;
            }
            // Allow the start_he (loop closure) even if it would normally
            // be excluded by U-turn prevention or visited check.
            if start_he == Some(idx) {
                return true;
            }
            // Exclude already-visited half-edges to prevent the walk
            // from crossing into edges consumed by a previous iteration.
            if visited.contains(&idx) {
                return false;
            }
            // Exclude the twin (U-turn prevention)
            if idx == twin_idx {
                return false;
            }
            // If we're on a self-loop, also exclude ourselves
            if is_self_loop && idx == he_idx {
                return false;
            }
            // When walking a NON-self-loop path, exclude self-loop
            // candidates from other faces. Self-loops whose face matches
            // the current face are boundary edges (e.g., intersection
            // circles on the barrel face); cross-registered self-loops
            // from other faces should still be excluded.
            if !is_self_loop {
                let twin_c = graph.half_edges[idx].twin;
                let is_candidate_self_loop =
                    graph.half_edges[idx].origin == graph.half_edges[twin_c].origin;
                // ZERO-EXTENT self-loops are point-contact records (a
                // tangency between two faces). They bound nothing, and
                // entering one hijacks the host loop: the closure
                // deferral held a ceiling loop open at its closure
                // node, the walk swallowed the tangency plus the twin
                // of the flank it had just traversed, and the spike's
                // flipped claim collided with the neighbouring flank
                // fragment's (measured on the uniform-width flush
                // arm). Their seed walks become one-half-edge
                // zero-extent fragments, which the degenerate cull
                // already removes.
                if is_candidate_self_loop {
                    use popkern_core::geometry::curve::CurveEval as _;
                    let extent = graph.half_edges[idx]
                        .curve
                        .as_ref()
                        .or_else(|| graph.half_edges[twin_c].curve.as_ref())
                        .map(|c| {
                            let origin = graph.nodes[graph.half_edges[idx].origin].position;
                            (1..8)
                                .map(|i| {
                                    let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 8.0;
                                    c.geometry.point_at(t).distance_to(&origin)
                                })
                                .fold(0.0, f64::max)
                        })
                        .unwrap_or(0.0);
                    if extent < 1e-9 {
                        return false;
                    }
                }
                // A self-loop candidate belongs to this walk when
                // EITHER half is assigned to this face: intersection
                // circles are cross-registered exactly so both of a
                // face's fragments (each side of the cut) can use one
                // direction each. Excluding by primary face alone
                // starved the second fragment of its circle boundary.
                let this_face = graph.half_edges[idx].face.as_ref() == Some(face)
                    || graph.half_edges[twin_c].face.as_ref() == Some(face);
                if is_candidate_self_loop && !this_face {
                    return false;
                }
                // INTERIOR-GRAZE FILTER: a candidate that dead-ends
                // (its far node has no onward half-edge on this face)
                // comes in two shapes that need opposite handling.
                // A HALF-SPLIT trim's far endpoint lies ON one of the
                // face's boundary edges — the boundary just was not
                // split there in this phase. Entering it is required:
                // the walk records the cut and downstream completes it
                // (measured: a cube bisected by a plane sheet stopped
                // splitting when these were filtered). An INTERIOR
                // graze record's far endpoint dead-ends strictly
                // inside the face — an embedded body's face tracing a
                // segment whose both sides bound the SAME kept
                // material (a flush gusset top's trace on the arm
                // wall). No result edge should exist there and no rule
                // may enter it: the junction rule already declined
                // these spurs, but the angular fallback picked them
                // anyway, weaving a 10-he unclosable outer that minted
                // a once-used edge.
                if !is_candidate_self_loop {
                    let far = graph.half_edges[twin_c].origin;
                    let has_onward = graph.nodes[far]
                        .outgoing
                        .iter()
                        .any(|&o| o != twin_c && o != idx && face_hes.contains(&o));
                    if !has_onward && !node_on_face_boundary(graph, far, face_hes, idx, twin_c) {
                        return false;
                    }
                }
            }
            true
        })
        .collect();

    if std::env::var_os("PK_WALK_TRACE").is_some() {
        let p = graph.nodes[end_node].position;
        eprintln!(
            "    [next] at ({:.3},{:.3},{:.3}) from he {he_idx}: candidates {candidates:?}",
            p.x, p.y, p.z
        );
    }
    if candidates.is_empty() {
        return None;
    }
    if candidates.len() == 1 {
        return Some(candidates[0]);
    }

    // §7.4 cut behavior: when the previous edge is a closed self-loop
    // on a Shared edge (an intersection circle that the boolean
    // inserted as a cut through this face), bounce back along the
    // twin of the edge that brought us TO this circle. The
    // intersection circle acts as a topological barrier — without
    // this rule, the walker continues past the cut and merges
    // multiple disjoint regions into a single mega-loop.
    //
    // The detection: `he_idx` is a closed self-loop (origin == end),
    // its source is `Shared`, and the SECOND-most-recent edge in
    // `current_loop` (the one before he_idx) has its twin in the
    // candidate set.
    if is_self_loop_shared(graph, he_idx)
        && current_loop.len() >= 2
        && let Some(prev_prev) = current_loop.get(current_loop.len().wrapping_sub(2))
    {
        let prev_prev_twin = graph.half_edges[*prev_prev].twin;
        if candidates.contains(&prev_prev_twin) {
            return Some(prev_prev_twin);
        }
    }

    // Prefer loop closure when start_he is a candidate AND all other
    // candidates are already visited. This avoids closing the loop
    // prematurely (skipping unvisited boundary edges like cap circles)
    // while still preventing extension into cross-registered duplicates.
    if let Some(sh) = start_he
        && candidates.contains(&sh)
    {
        // Only defer closure if there are unvisited candidates whose face
        // matches (actual boundary edges, not cross-registered duplicates).
        // Only defer closure for unvisited NON-SHARED face-matching
        // candidates (actual boundary edges like cap circles). Shared
        // edges (even face-matching ones from cross-registration) should
        // not prevent loop closure.
        let unvisited_face_matching = candidates.iter().any(|&c| {
            c != sh
                && !visited.contains(&c)
                && graph.half_edges[c].face.as_ref() == Some(face)
                && !matches!(
                    graph.half_edges[c].source,
                    super::graph::EdgeSource::Shared { .. }
                )
        });
        // An unvisited SELF-LOOP intersection curve at this node is a
        // topological barrier that bounds this fragment — closing the
        // loop early would leave the fragment without its cut
        // boundary (the wall of a boss/plate union lost its seam
        // circle and the circle became a stray one-edge fragment).
        let unvisited_cut_here = candidates.iter().any(|&c| {
            let twin_c = graph.half_edges[c].twin;
            c != sh
                && !visited.contains(&c)
                && graph.half_edges[c].origin == graph.half_edges[twin_c].origin
                && (graph.half_edges[c].face.as_ref() == Some(face)
                    || graph.half_edges[twin_c].face.as_ref() == Some(face))
        });
        if !unvisited_face_matching && !unvisited_cut_here {
            return Some(sh);
        }
    }

    // At seam T-junctions on curved surfaces (where seam continuation
    // competes with a self-loop circle), prefer the unvisited self-loop.
    // Self-loops with matching face are actual boundary edges; the seam
    // is just the internal spine of a periodic surface.
    let self_loop_candidates: Vec<usize> = candidates
        .iter()
        .copied()
        .filter(|&c| {
            !visited.contains(&c) && {
                let twin_c = graph.half_edges[c].twin;
                graph.half_edges[c].origin == graph.half_edges[twin_c].origin
            }
        })
        .collect();
    if !self_loop_candidates.is_empty() && self_loop_candidates.len() < candidates.len() {
        if self_loop_candidates.len() == 1 {
            return Some(self_loop_candidates[0]);
        }
        return pick_by_angle(graph, he_idx, &self_loop_candidates, face, brep_a, brep_b);
    }

    // TOPOLOGICAL JUNCTION RULE (convention-free): a transversal
    // boundary-to-boundary trim separates two fragments, so every
    // loop through one of its junction nodes must ALTERNATE — a walk
    // arriving along the BOUNDARY turns onto the trim (it bounds the
    // fragment on this side), and a walk arriving along the TRIM
    // continues onto the boundary. Continuing straight past a
    // crossing trim encloses both fragments in one loop (measured: a
    // rest-plane trim on a reflex block/plate junction wall was left
    // unused, both wall regions merged into one 10-half-edge
    // fragment carrying removed-region geometry, and the trim's own
    // walk dead-ended). This is exact whenever the junction has
    // exactly one non-self-loop Shared candidate against non-shared
    // boundary continuations — no surface normal or angle involved.
    // Branch points (several trims) and circle junctions keep the
    // angular pick below.
    let incoming_shared = matches!(
        graph.half_edges[he_idx].source,
        super::graph::EdgeSource::Shared { .. }
    ) && !is_self_loop;
    if !incoming_shared {
        let shared_straight: Vec<usize> = candidates
            .iter()
            .copied()
            .filter(|&c| {
                let twin_c = graph.half_edges[c].twin;
                let self_loop_c = graph.half_edges[c].origin == graph.half_edges[twin_c].origin;
                !self_loop_c
                    && matches!(
                        graph.half_edges[c].source,
                        super::graph::EdgeSource::Shared { .. }
                    )
            })
            .collect();
        if shared_straight.len() == 1 && shared_straight.len() < candidates.len() {
            // SPUR GUARD: a closed loop can never traverse a trim
            // whose far node has no onward half-edge on this face —
            // entering it forces a forbidden twin u-turn. A grazing
            // record (a flush gusset's top-edge line lying ON the
            // wall face) dead-ends mid-face; the rule walked the
            // outer boundary straight into it and the walk stuck.
            // Only take the trim when it connects onward.
            let c = shared_straight[0];
            let twin_c = graph.half_edges[c].twin;
            let far = graph.half_edges[twin_c].origin;
            let has_onward = graph.nodes[far]
                .outgoing
                .iter()
                .any(|&o| o != twin_c && o != c && face_hes.contains(&o));
            if has_onward {
                return Some(c);
            }
        }
    }

    // Use per-face angular comparison for all multi-candidate cases.
    pick_by_angle(graph, he_idx, &candidates, face, brep_a, brep_b)
}

#[cfg(test)]
mod tests;
