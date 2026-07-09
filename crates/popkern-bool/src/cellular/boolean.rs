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
    split_graph_edge,
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
    let fragments = walk_all_fragments_ext(
        &graph,
        original_a,
        original_b,
        !intersection_result.coincidence.face_overlaps.is_empty(),
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
    )?;

    // Build provenance map (available via boolean_with_provenance).
    let _provenance =
        super::assemble::build_face_provenance(&result, &kept, &fragments, &selections);

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
                graph.half_edges[a].twin == fragment.outer_loop[1]
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
        .map(|s| s.geometry.clone());
        let sample = frag_surface
            .and_then(|surf| fragment_interior_sample(graph, fragment, &surf))
            .unwrap_or_else(|| sample_fragment_point(graph, &fragment.outer_loop));

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
        let class = if !needs_offset_pair {
            robust_point_in_solid_final(&sample, other_brep, other_aabbs, tolerance)
        } else {
            use popkern_core::geometry::surface::SurfaceEval;
            let frag_face = match fragment.face.operand {
                Operand::A => brep_a.faces.get(fragment.face.face),
                Operand::B => brep_b.faces.get(fragment.face.face),
            };
            let normal = frag_face
                .and_then(|f| f.surface.as_ref())
                .and_then(|sof| {
                    popkern_core::geometry::uv_projection::project_point_to_uv(
                        &sof.geometry,
                        &sample,
                    )
                    .map(|(u, v)| sof.geometry.normal_at(u, v))
                })
                .and_then(|n| n.normalized());
            match normal {
                Some(n) => {
                    let eps = (tolerance.linear() * 1.0e3).max(1.0e-3);
                    let p_plus = sample + n * eps;
                    let p_minus = sample - n * eps;
                    let c_plus =
                        robust_point_in_solid_final(&p_plus, other_brep, other_aabbs, tolerance);
                    let c_minus =
                        robust_point_in_solid_final(&p_minus, other_brep, other_aabbs, tolerance);
                    if std::env::var_os("PK_SEL_TRACE").is_some() {
                        eprintln!(
                            "      [pair] {:?} outer={} n=({:.2},{:.2},{:.2}) eps={eps:.2e} c_plus={c_plus:?} c_minus={c_minus:?}",
                            fragment.face,
                            fragment.outer_loop.len(),
                            n.x,
                            n.y,
                            n.z
                        );
                    }
                    match (c_plus, c_minus) {
                        (PointClassification::Inside, PointClassification::Inside) => {
                            PointClassification::Inside
                        }
                        (PointClassification::Outside, PointClassification::Outside) => {
                            PointClassification::Outside
                        }
                        _ => PointClassification::OnBoundary,
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
                "  [sel] {:?} outer={} sample=({:.2},{:.2},{:.2}) class={class:?}",
                fragment.face,
                fragment.outer_loop.len(),
                sample.x,
                sample.y,
                sample.z
            );
        }
        let (keep, flip) = match op {
            BoolOperation::Union => (
                class == PointClassification::Outside
                    || (own_operand == Operand::A && class == PointClassification::OnBoundary),
                false,
            ),
            BoolOperation::Intersection => (
                class == PointClassification::Inside
                    || (own_operand == Operand::A && class == PointClassification::OnBoundary),
                false,
            ),
            BoolOperation::Subtraction => match own_operand {
                Operand::A => (
                    class == PointClassification::Outside
                        || class == PointClassification::OnBoundary,
                    false,
                ),
                Operand::B => (class == PointClassification::Inside, true),
            },
        };

        selections.push(FragmentSelection { keep, flip });
    }

    Ok(selections)
}

/// Sample a point from a fragment's outer loop for classification.
/// Interior sample point for a fragment, aware of inner loops and
/// closed-rim outers. Works in the UV chart of the fragment's face:
/// builds pcurve-sampled polygons for the outer and inner loops,
/// scans across the outer polygon's centroid row, and picks the
/// midpoint of the widest interval covered by the outer but not by
/// any inner. Falls back to `None` when pcurves are missing.
fn fragment_interior_sample(
    graph: &IntersectionGraph,
    fragment: &crate::merge::walk::Fragment,
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
) -> Option<popkern_core::geometry::Point3> {
    use popkern_core::geometry::surface::SurfaceEval;

    let needs_care = !fragment.inner_loops.is_empty()
        || fragment.outer_loop.iter().any(|&h| {
            let he = &graph.half_edges[h];
            he.origin == graph.half_edges[he.twin].origin
        });
    if !needs_care {
        return None;
    }

    // Stored pcurves may live on ANY period branch of a periodic
    // chart (a cylinder pcurve at u−2π is the same geometry), while
    // projected samples come back on the principal branch. Mixing
    // branches in one scanline polygon opens a phantom interval
    // spanning the branch gap — the widest-interval pick then lands
    // its midpoint on the far side of the cylinder (measured: an
    // embedded sliver of a flush gusset's arc wall sampled at the
    // antipode, classified Outside, and survived to fight the step
    // face for its seam edge). Shift each half-edge's sample run by
    // whole periods so it lands adjacent to the previous run.
    let u_period = surface.parameter_periodicity().0;
    let branch_aligned = |poly: &mut Vec<(f64, f64)>, run: Vec<(f64, f64)>| {
        if let (Some(p), Some(&(lu, _)), Some(&(fu, _))) = (u_period, poly.last(), run.first()) {
            let k = ((lu - fu) / p).round();
            if k != 0.0 {
                poly.extend(run.into_iter().map(|(u, v)| (u + k * p, v)));
                return;
            }
        }
        poly.extend(run);
    };

    let sample_loop = |hes: &[usize]| -> Option<Vec<(f64, f64)>> {
        use popkern_core::geometry::curve::CurveEval;
        use popkern_core::geometry::uv_projection::{
            project_point_to_uv, project_point_to_uv_near,
        };
        let trace = std::env::var_os("PK_SEL_TRACE").is_some();
        let mut poly: Vec<(f64, f64)> = Vec::new();
        for &raw in hes {
            // Pcurves live in their OWNER face's chart. A fragment
            // loop may carry either direction of an edge (the walker
            // keeps the first-walked direction), so read the side
            // that belongs to THIS face — the twin's pcurve is in
            // the other body's chart and produced garbage scanline
            // polygons (measured: a counterbore floor's star-hole
            // fragment sampled its own hole centre and was dropped).
            let h = if graph.half_edges[raw].face.as_ref() == Some(&fragment.face) {
                raw
            } else {
                let twin = graph.half_edges[raw].twin;
                if graph.half_edges[twin].face.as_ref() == Some(&fragment.face) {
                    twin
                } else {
                    raw
                }
            };
            let mut run: Vec<(f64, f64)> = Vec::new();
            let kind;
            if let Some(pc) = graph.half_edges[h].pcurve.as_ref() {
                kind = "PCURVE";
                for i in 0..8 {
                    let t = pc.t_start + (pc.t_end - pc.t_start) * i as f64 / 8.0;
                    let p = pc.point_at(t);
                    run.push((p.x, p.y));
                }
            } else if let Some(c) = graph.half_edges[h].curve.as_ref() {
                // No pcurve (split edges): project curve samples,
                // with continuation to stay on one branch.
                kind = "CURVE";
                let ceval: &dyn CurveEval = &c.geometry;
                for i in 0..8 {
                    let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 8.0;
                    let p = ceval.point_at(t);
                    let uv = match run.last().or_else(|| poly.last()) {
                        Some(&(pu, pv)) => project_point_to_uv_near(surface, &p, pu, pv)
                            .or_else(|| project_point_to_uv(surface, &p)),
                        None => project_point_to_uv(surface, &p),
                    }?;
                    run.push(uv);
                }
            } else {
                // Straight graph edge without geometry: project the
                // endpoints.
                kind = "ENDPOINTS";
                let a = graph.nodes[graph.half_edges[h].origin].position;
                let b = graph.nodes[graph.half_edges[graph.half_edges[h].twin].origin].position;
                for p in [a, b] {
                    run.push(project_point_to_uv(surface, &p)?);
                }
            }
            if trace && !run.is_empty() {
                eprintln!(
                    "      [selpoly] he {raw} (used {h}) src {:?} {kind} first=({:.3},{:.3}) last=({:.3},{:.3})",
                    graph.half_edges[h].source,
                    run[0].0,
                    run[0].1,
                    run[run.len() - 1].0,
                    run[run.len() - 1].1
                );
            }
            branch_aligned(&mut poly, run);
        }
        (poly.len() >= 3).then_some(poly)
    };
    let outer = sample_loop(&fragment.outer_loop)?;
    let inners: Vec<Vec<(f64, f64)>> = fragment
        .inner_loops
        .iter()
        .filter_map(|l| sample_loop(l))
        .collect();
    if inners.len() != fragment.inner_loops.len() {
        return None;
    }

    // Scan row: outer centroid v.
    let cv = outer.iter().map(|p| p.1).sum::<f64>() / outer.len() as f64;
    let crossings_at = |poly: &[(f64, f64)], v: f64| -> Vec<f64> {
        let mut xs = Vec::new();
        let n = poly.len();
        for i in 0..n {
            let (u0, v0) = poly[i];
            let (u1, v1) = poly[(i + 1) % n];
            if (v0 > v) != (v1 > v) {
                xs.push(u0 + (v - v0) / (v1 - v0) * (u1 - u0));
            }
        }
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        xs
    };
    let outer_x = crossings_at(&outer, cv);
    if outer_x.len() < 2 {
        return None;
    }
    // Candidate intervals inside the outer polygon (even-odd pairs).
    let mut best: Option<(f64, f64)> = None; // (width, mid_u)
    for pair in outer_x.chunks(2) {
        let [a, b] = pair else { continue };
        // Subtract inner coverage: split [a,b] at inner crossings.
        let mut cuts = vec![*a, *b];
        for inner in &inners {
            for x in crossings_at(inner, cv) {
                if x > *a && x < *b {
                    cuts.push(x);
                }
            }
        }
        cuts.sort_by(|p, q| p.partial_cmp(q).unwrap_or(std::cmp::Ordering::Equal));
        for w in cuts.windows(2) {
            let mid = 0.5 * (w[0] + w[1]);
            // Point-in-polygon (even-odd) for the mid against inners.
            let inside_inner = inners.iter().any(|inner| {
                let xs = crossings_at(inner, cv);
                xs.iter().filter(|&&x| x < mid).count() % 2 == 1
            });
            if !inside_inner {
                let width = w[1] - w[0];
                if best.map(|(bw, _)| width > bw).unwrap_or(true) {
                    best = Some((width, mid));
                }
            }
        }
    }
    let (_, mid_u) = best?;
    if std::env::var_os("PK_SEL_TRACE").is_some() {
        eprintln!("      [selpoly] scan cv={cv:.3} outer_x={outer_x:?} -> mid_u={mid_u:.3}");
    }
    Some(surface.point_at(mid_u, cv))
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
                let start = graph.nodes[he.origin].position;
                let twin_idx = he.twin;
                let end = graph.nodes[graph.half_edges[twin_idx].origin].position;
                return popkern_core::geometry::Point3::new(
                    (start.x + end.x) * 0.5,
                    (start.y + end.y) * 0.5,
                    (start.z + end.z) * 0.5,
                );
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
    let closed = c.point_at(t0).distance_to(&c.point_at(t1)) < tolerance.linear() * 100.0;
    if !closed {
        return record.clone();
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
        if c.point_at(t_hit).distance_to(&hit.point) > tolerance.linear() * 100.0 {
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
            control_points.extend(a.control_points.iter().skip(1).copied());
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
