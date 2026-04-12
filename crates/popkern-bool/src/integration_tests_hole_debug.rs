//! Debug tests for cylinder-through-cube hole subtraction.
//!
//! The stress test revealed that cube-cylinder subtraction at tight tolerance
//! (1e-6) produces invalid topology, while at coarse tolerance (1e-2) it works.
//! These tests trace the issue.

use popkern_core::geometry::{Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn default_tol() -> Tolerance {
    Tolerance::default() // 1e-6
}

fn coarse_tol() -> Tolerance {
    Tolerance::default().with_linear(1e-2)
}

fn medium_tol() -> Tolerance {
    Tolerance::default().with_linear(1e-4)
}

fn make_cube_with_hole(tol: Tolerance) -> Result<BRep, crate::error::BoolError> {
    let cube = BRep::unit_cube(tol).unwrap();
    let cyl = BRep::cylinder(0.2, 2.0, tol)
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    crate::subtraction(cube, cyl, tol, None)
}

#[test]
fn hole_at_coarse_tol_works() {
    let result = make_cube_with_hole(coarse_tol()).expect("coarse tol should work");
    let v = result.vertices.len();
    let e = result.edges.len();
    let f = result.faces.len();
    eprintln!("[coarse] V={v} E={e} F={f}");
    assert!(f >= 6, "should have at least 6 faces");
    let inner_loops: usize = result.faces.values().map(|f| f.inner_loops.len()).sum();
    eprintln!("[coarse] inner_loops={inner_loops}");
}

#[test]
fn hole_at_coarse_tol_passes_strict_post_check() {
    // After the merge-phase cut-behavior fix, the cube-minus-cylinder
    // case should produce a topologically clean genus-1 result
    // (V-E+F-(L-F)=0). This test uses boolean_strict to enforce
    // post-check validation; if the merge regresses, this fires.
    let tol = coarse_tol();
    let cube = BRep::unit_cube(tol).unwrap();
    let cyl = BRep::cylinder(0.2, 2.0, tol)
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    let result = crate::boolean_strict(crate::BoolOperation::Subtraction, &cube, &cyl, tol, false);
    assert!(
        result.is_ok(),
        "cube-minus-cylinder should pass strict post-check, got {result:?}"
    );
}

#[test]
#[ignore = "diagnostic-only — see docs/reports/coverage/booleans-vs-requirements.md"]
fn diagnose_merge_fragments() {
    // Run the merge phases manually so we can see exactly how many
    // fragments walk_all_fragments produces for the cylinder side
    // face, and what they look like.
    use crate::intersect::find_intersections;
    use crate::merge::graph::{IntersectionGraph, Operand};
    use crate::merge::graph_build::insert_operand_topology;
    use crate::merge::graph_intersect::insert_intersection_curves;
    use crate::merge::graph_split::{
        EndpointClass, classify_endpoint, remove_collinear_duplicates, split_edges_at_shared_nodes,
        split_graph_edge,
    };
    use crate::merge::radial::sort_all_radial;
    use crate::merge::walk::walk_all_fragments;

    let tol = coarse_tol();
    let cube = BRep::unit_cube(tol).unwrap();
    let cyl = BRep::cylinder(0.2, 2.0, tol)
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();

    let intersection_result = find_intersections(&cube, &cyl, &tol)
        .map_err(|d| format!("intersection: {d:?}"))
        .unwrap();
    eprintln!(
        "[merge-diag] intersection records: {}",
        intersection_result.records.len()
    );

    let snap_tol = tol.linear() * 2.0;
    let mut graph = IntersectionGraph::new();
    insert_operand_topology(&mut graph, &cube, Operand::A, snap_tol);
    insert_operand_topology(&mut graph, &cyl, Operand::B, snap_tol);
    eprintln!(
        "[merge-diag] after operand insert: nodes={} half_edges={}",
        graph.nodes.len(),
        graph.half_edges.len()
    );

    for hit in &intersection_result.edge_face_hits {
        let endpoint_class = classify_endpoint(&graph, &hit.point, snap_tol);
        if let EndpointClass::OnEdge(he_idx, t_param) = endpoint_class {
            split_graph_edge(&mut graph, he_idx, hit.point, t_param);
        }
    }
    eprintln!(
        "[merge-diag] after edge-face hit splits: nodes={} half_edges={}",
        graph.nodes.len(),
        graph.half_edges.len()
    );

    insert_intersection_curves(&mut graph, &intersection_result.records, snap_tol);
    eprintln!(
        "[merge-diag] after intersection curve insert: nodes={} half_edges={}",
        graph.nodes.len(),
        graph.half_edges.len()
    );

    if intersection_result
        .coincidence
        .face_overlaps
        .iter()
        .any(|o| !o.full_overlap)
    {
        split_edges_at_shared_nodes(&mut graph, snap_tol);
        remove_collinear_duplicates(&mut graph);
    }

    sort_all_radial(&mut graph);

    let fragments = walk_all_fragments(&graph, &cube, &cyl);
    eprintln!("[merge-diag] fragments produced: {}", fragments.len());
    for (i, frag) in fragments.iter().enumerate() {
        eprintln!(
            "[merge-diag]   fragment[{}] face={:?} outer_loop_len={} inner_loops={}",
            i,
            frag.face,
            frag.outer_loop.len(),
            frag.inner_loops.len()
        );
        // Print the half-edge sequence for the cylinder side face fragment.
        if frag.face.operand == Operand::B {
            for (j, &he_idx) in frag.outer_loop.iter().enumerate() {
                let he = &graph.half_edges[he_idx];
                let twin = &graph.half_edges[he.twin];
                let start = graph.nodes[he.origin].position;
                let end = graph.nodes[twin.origin].position;
                eprintln!(
                    "[merge-diag]     [{}] he={} ({:.2},{:.2},{:.2}) -> ({:.2},{:.2},{:.2}) src={:?}",
                    j, he_idx, start.x, start.y, start.z, end.x, end.y, end.z, he.source,
                );
            }
        }
    }
}

#[test]
#[ignore = "diagnostic-only — see docs/reports/coverage/booleans-vs-requirements.md"]
fn diagnose_input_topology_counts() {
    // Reference counts for the bodies that go INTO the failing
    // boolean. If the boolean drops or merges entities, comparing
    // input/output deltas pinpoints which phase is responsible.
    let tol = coarse_tol();
    let cube = BRep::unit_cube(tol).unwrap();
    let cyl = BRep::cylinder(0.2, 2.0, tol)
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    eprintln!(
        "[input] cube V={} E={} F={} L={}",
        cube.vertices.len(),
        cube.edges.len(),
        cube.faces.len(),
        cube.loops.len()
    );
    eprintln!(
        "[input] cyl  V={} E={} F={} L={}",
        cyl.vertices.len(),
        cyl.edges.len(),
        cyl.faces.len(),
        cyl.loops.len()
    );
    let cube_lhs = cube.vertices.len() as i64 - cube.edges.len() as i64 + cube.faces.len() as i64
        - (cube.loops.len() as i64 - cube.faces.len() as i64);
    let cyl_lhs = cyl.vertices.len() as i64 - cyl.edges.len() as i64 + cyl.faces.len() as i64
        - (cyl.loops.len() as i64 - cyl.faces.len() as i64);
    eprintln!("[input] cube Euler lhs = {cube_lhs} (expected 2)");
    eprintln!("[input] cyl  Euler lhs = {cyl_lhs} (expected 2)");
}

#[test]
#[ignore = "diagnostic-only — see docs/reports/coverage/booleans-vs-requirements.md"]
fn diagnose_hole_euler_violation() {
    // Run the cube-minus-cylinder case at coarse tol and report the
    // full topology counts. The cube-with-hole result should be a
    // genus-1 solid: V - E + F - (L - F) = 0 = 2(S - G) with S=1, G=1.
    // Anything else is corruption that this diagnostic surfaces.
    let result = make_cube_with_hole(coarse_tol()).expect("coarse tol should work");
    let v = result.vertices.len() as i64;
    let e = result.edges.len() as i64;
    let f = result.faces.len() as i64;
    let l = result.loops.len() as i64;
    let s = result.shells.len() as i64;
    let chi = v - e + f;
    let lhs = v - e + f - (l - f);
    let inner_loops: usize = result.faces.values().map(|fi| fi.inner_loops.len()).sum();

    eprintln!("[diag] V={v} E={e} F={f} L={l} S={s}");
    eprintln!("[diag] chi=V-E+F={chi}");
    eprintln!("[diag] lhs (full Euler) = V-E+F-(L-F) = {lhs}");
    eprintln!("[diag] expected lhs for genus-0 single-shell: 2");
    eprintln!("[diag] expected lhs for genus-1 single-shell: 0");
    eprintln!("[diag] inner_loops total = {inner_loops}");
    eprintln!("[diag] outer_loops = F = {f}");
    eprintln!(
        "[diag] L should equal F + inner_loops = {}",
        f + inner_loops as i64
    );
    eprintln!("[diag] L actual = {l}");

    // Walk and report each face's loop counts.
    for (i, (_, face)) in result.faces.iter().enumerate() {
        eprintln!(
            "[diag] face[{}] inner_loops={} surface={:?}",
            i,
            face.inner_loops.len(),
            face.surface
                .as_ref()
                .map(|s| std::any::type_name_of_val(&s.geometry)),
        );
    }

    // Walk every edge and report which faces use it.
    eprintln!("[diag] edge -> face usage:");
    for (eh, edge) in &result.edges {
        let he = &result.half_edges[edge.half_edge];
        let twin = &result.half_edges[he.twin];
        let face_a = result.loops.get(he.loop_).map(|l| l.face);
        let face_b = result.loops.get(twin.loop_).map(|l| l.face);
        eprintln!(
            "[diag]   edge {:?}: face_a={:?} face_b={:?}",
            eh, face_a, face_b
        );
    }

    // Walk every vertex and report its position.
    eprintln!("[diag] vertices:");
    for (vh, v) in &result.vertices {
        eprintln!(
            "[diag]   {:?}: ({:.3}, {:.3}, {:.3})",
            vh, v.position.x, v.position.y, v.position.z
        );
    }

    // Per-face loop counts: outer loop length + each inner loop length.
    eprintln!("[diag] per-face loop layout:");
    for (fh, face) in &result.faces {
        let outer_len = walk_loop_length(&result, face.outer_loop);
        eprintln!(
            "[diag]   face {:?}: outer_loop_len={} inner_loops={}",
            fh,
            outer_len,
            face.inner_loops.len()
        );
        for &inner in &face.inner_loops {
            let inner_len = walk_loop_length(&result, inner);
            eprintln!("[diag]     inner_loop len={inner_len}");
        }
    }
}

#[test]
#[ignore = "diagnostic-only — see docs/reports/coverage/booleans-vs-requirements.md"]
fn diagnose_orphan_loops() {
    use std::collections::HashSet;
    let result = make_cube_with_hole(coarse_tol()).expect("coarse tol should work");

    // Collect all loop handles referenced by faces.
    let mut owned: HashSet<popkern_core::topology::LoopHandle> = HashSet::new();
    for (_, face) in &result.faces {
        owned.insert(face.outer_loop);
        for &inner in &face.inner_loops {
            owned.insert(inner);
        }
    }

    eprintln!("[diag] total loops in slotmap: {}", result.loops.len());
    eprintln!("[diag] loops referenced by faces: {}", owned.len());

    // Walk every loop in the slotmap and report whether it's owned.
    for (lh, ld) in &result.loops {
        let face_of = ld.face;
        let referenced = owned.contains(&lh);
        eprintln!(
            "[diag]   loop {:?}: face={:?} referenced_by_face={}",
            lh, face_of, referenced
        );
        if !referenced {
            // Orphan loop — walk it and report length and starting half-edge.
            let len = walk_loop_length(&result, lh);
            eprintln!("[diag]     ORPHAN loop length={len}");
        }
    }
}

#[test]
#[ignore = "diagnostic-only — see docs/reports/coverage/booleans-vs-requirements.md"]
fn diagnose_loop_repair_outcome() {
    // Run the boolean and then explicitly run loop_repair on the
    // result. Print the report so we can see whether the pass
    // identified stranded half-edges and whether the removal stuck.
    let mut result = make_cube_with_hole(coarse_tol()).expect("coarse tol should work");
    let report = crate::postprocess::loop_repair::repair_stranded_loops(&mut result);
    eprintln!("[diag] loop_repair report: {report:?}");
    let v = result.vertices.len();
    let e = result.edges.len();
    let f = result.faces.len();
    let l = result.loops.len();
    eprintln!("[diag] post-repair V={v} E={e} F={f} L={l}");
}

#[test]
#[ignore = "diagnostic-only — see docs/reports/coverage/booleans-vs-requirements.md"]
fn diagnose_unreachable_half_edges() {
    use std::collections::HashSet;
    let result = make_cube_with_hole(coarse_tol()).expect("coarse tol should work");

    let mut reachable: HashSet<popkern_core::topology::HalfEdgeHandle> = HashSet::new();
    for (_, face) in &result.faces {
        collect_loop(&result, face.outer_loop, &mut reachable);
        for &inner in &face.inner_loops {
            collect_loop(&result, inner, &mut reachable);
        }
    }

    eprintln!("[diag] half-edges in slotmap = {}", result.half_edges.len());
    eprintln!(
        "[diag] half-edges reachable from face loops = {}",
        reachable.len()
    );

    let mut leftover_count = 0;
    for (heh, he) in &result.half_edges {
        if !reachable.contains(&heh) {
            leftover_count += 1;
            let origin_pos = result.vertices.get(he.origin).map(|v| v.position);
            eprintln!(
                "[diag]   leftover half-edge {:?}: origin={:?}, loop_={:?}",
                heh, origin_pos, he.loop_
            );
        }
    }
    eprintln!("[diag] total leftover half-edges = {leftover_count}");

    // Per-face: count expected (sum of self-loop edges * 2 + manifold edges)
    // versus actual outer-loop reach to find loop-graph inconsistencies.
    eprintln!("[diag] per-face half-edge usage:");
    for (fh, face) in &result.faces {
        let mut declared = 0;
        for (_, edge) in &result.edges {
            let he = &result.half_edges[edge.half_edge];
            let twin = &result.half_edges[he.twin];
            let face_a = result.loops.get(he.loop_).map(|l| l.face);
            let face_b = result.loops.get(twin.loop_).map(|l| l.face);
            if face_a == Some(fh) {
                declared += 1;
            }
            if face_b == Some(fh) {
                declared += 1;
            }
        }
        let outer_walked = walk_loop_length(&result, face.outer_loop);
        let inner_walked: usize = face
            .inner_loops
            .iter()
            .map(|&l| walk_loop_length(&result, l))
            .sum();
        eprintln!(
            "[diag]   face {fh:?}: declared_he={declared} outer_walked={outer_walked} inner_walked={inner_walked} delta={}",
            declared as i64 - outer_walked as i64 - inner_walked as i64
        );
    }
}

fn collect_loop(
    brep: &BRep,
    loop_h: popkern_core::topology::LoopHandle,
    out: &mut std::collections::HashSet<popkern_core::topology::HalfEdgeHandle>,
) {
    let Some(loop_data) = brep.loops.get(loop_h) else {
        return;
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        out.insert(cursor);
        if let Some(he) = brep.half_edges.get(cursor) {
            cursor = he.next;
        } else {
            break;
        }
        if cursor == start {
            break;
        }
    }
}

fn walk_loop_length(brep: &BRep, loop_h: popkern_core::topology::LoopHandle) -> usize {
    let Some(loop_data) = brep.loops.get(loop_h) else {
        return 0;
    };
    let start = loop_data.half_edge;
    let mut count = 0;
    let mut cursor = start;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        count += 1;
        if let Some(he) = brep.half_edges.get(cursor) {
            cursor = he.next;
        } else {
            break;
        }
        if cursor == start {
            break;
        }
    }
    count
}

#[test]
fn hole_at_medium_tol() {
    let result = make_cube_with_hole(medium_tol());
    match result {
        Ok(brep) => {
            let v = brep.vertices.len();
            let e = brep.edges.len();
            let f = brep.faces.len();
            let errors = brep.validate_all();
            eprintln!("[medium] V={v} E={e} F={f} errors={}", errors.len());
            for err in &errors {
                eprintln!("[medium]   {err}");
            }
        }
        Err(e) => {
            eprintln!("[medium] FAILED: {e}");
        }
    }
}

#[test]
fn hole_at_default_tol() {
    let result = make_cube_with_hole(default_tol());
    match result {
        Ok(brep) => {
            let v = brep.vertices.len();
            let e = brep.edges.len();
            let f = brep.faces.len();
            let errors = brep.validate_all();
            eprintln!("[default] V={v} E={e} F={f} errors={}", errors.len());
            for err in &errors {
                eprintln!("[default]   {err}");
            }
            // Trace each loop to find the gap.
            use popkern_core::geometry::CurveEval;
            for (lh, ld) in &brep.loops {
                let start_he = ld.half_edge;
                let mut cursor = start_he;
                let mut first_pos = None;
                let mut end_pos = None;
                let mut count = 0;
                for _ in 0..100 {
                    let he = &brep.half_edges[cursor];
                    let pos = brep.vertices[he.origin].position;
                    if first_pos.is_none() {
                        first_pos = Some(pos);
                    }
                    let edge = &brep.edges[he.edge];
                    end_pos = if let Some(c) = &edge.curve {
                        let is_canonical = edge.half_edge == cursor;
                        Some(if is_canonical {
                            c.geometry.point_at(c.t_end)
                        } else {
                            c.geometry.point_at(c.t_start)
                        })
                    } else {
                        let next_he = &brep.half_edges[he.next];
                        Some(brep.vertices[next_he.origin].position)
                    };
                    count += 1;
                    cursor = he.next;
                    if cursor == start_he {
                        break;
                    }
                }
                if let (Some(fp), Some(ep)) = (first_pos, end_pos) {
                    let gap = fp.distance_to(&ep);
                    if gap > 1e-6 {
                        eprintln!("[gap] loop={lh:?} count={count} gap={gap:.6}");
                        eprintln!(
                            "  first=({:.4},{:.4},{:.4}) last_end=({:.4},{:.4},{:.4})",
                            fp.x, fp.y, fp.z, ep.x, ep.y, ep.z
                        );
                        // Also dump the last edge's curve info
                        cursor = start_he;
                        for _ in 0..count {
                            let he = &brep.half_edges[cursor];
                            let edge = &brep.edges[he.edge];
                            let has_curve = edge.curve.is_some();
                            let origin_pos = brep.vertices[he.origin].position;
                            eprintln!(
                                "  he origin=({:.4},{:.4},{:.4}) curve={has_curve}",
                                origin_pos.x, origin_pos.y, origin_pos.z
                            );
                            cursor = he.next;
                            if cursor == start_he {
                                break;
                            }
                        }
                    }
                }
            }

            // Trace vertex valence and face topology.
            for (vh, vertex) in &brep.vertices {
                let valence = brep
                    .half_edges
                    .iter()
                    .filter(|(_, he)| he.origin == vh)
                    .count();
                if valence < 2 {
                    let p = vertex.position;
                    eprintln!(
                        "[valence] vertex {vh:?} valence={valence} pos=({:.4},{:.4},{:.4})",
                        p.x, p.y, p.z
                    );
                    // Show which half-edges reference this vertex.
                    for (heh, he) in &brep.half_edges {
                        if he.origin == vh {
                            let edge = &brep.edges[he.edge];
                            let has_curve = edge.curve.is_some();
                            let twin_origin = brep.half_edges[he.twin].origin;
                            let is_self_loop = he.origin == twin_origin;
                            eprintln!(
                                "  he={heh:?} twin_origin={twin_origin:?} self_loop={is_self_loop} loop={:?} curve={has_curve}",
                                he.loop_
                            );
                        }
                    }
                }
            }
            // Trace face loops.
            for (fh, face) in &brep.faces {
                let outer_count = {
                    let mut c = 0;
                    let start = brep.loops[face.outer_loop].half_edge;
                    let mut cur = start;
                    for _ in 0..200 {
                        c += 1;
                        cur = brep.half_edges[cur].next;
                        if cur == start {
                            break;
                        }
                    }
                    c
                };
                let inner_counts: Vec<usize> = face
                    .inner_loops
                    .iter()
                    .map(|&lh| {
                        let mut c = 0;
                        let start = brep.loops[lh].half_edge;
                        let mut cur = start;
                        for _ in 0..200 {
                            c += 1;
                            cur = brep.half_edges[cur].next;
                            if cur == start {
                                break;
                            }
                        }
                        c
                    })
                    .collect();
                eprintln!(
                    "[face] {fh:?} outer={outer_count} inner_loops={} inner_counts={inner_counts:?}",
                    face.inner_loops.len()
                );
            }

            // Check intersection curve count
            let tol = default_tol();
            let cube = BRep::unit_cube(tol).unwrap();
            let cyl = BRep::cylinder(0.2, 2.0, tol)
                .unwrap()
                .transform(&Transform3::translation(0.5, 0.5, -0.5))
                .unwrap();
            let isect = crate::intersect::find_intersections(&cube, &cyl, &tol);
            match isect {
                Ok(result) => {
                    eprintln!(
                        "[default] intersection: {} curves, {} coincident overlaps",
                        result.records.len(),
                        result.coincidence.face_overlaps.len()
                    );
                    for (i, rec) in result.records.iter().enumerate() {
                        let (t0, t1) = rec.curve.curve_3d.parameter_range();
                        use popkern_core::geometry::CurveEval;
                        let p0 = rec.curve.curve_3d.point_at(t0);
                        let p1 = rec.curve.curve_3d.point_at(t1);
                        let closed = rec.curve.is_closed;
                        eprintln!(
                            "[default]   curve {i}: face_a={:?} face_b={:?} closed={closed}",
                            rec.face_a, rec.face_b
                        );
                        eprintln!(
                            "[default]     start=({:.3},{:.3},{:.3}) end=({:.3},{:.3},{:.3})",
                            p0.x, p0.y, p0.z, p1.x, p1.y, p1.z
                        );
                    }
                }
                Err(diag) => {
                    eprintln!("[default] intersection FAILED: {:?}", diag.error);
                }
            }
        }
        Err(e) => {
            eprintln!("[default] FAILED: {e}");
        }
    }
}

// ── Overlapping cubes subtraction (stress test 3) ────────────────────

#[test]
fn overlap_cubes_sub_debug() {
    let tol = default_tol();
    let a = BRep::unit_cube(tol).unwrap();
    let b = BRep::unit_cube(tol)
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, 0.0))
        .unwrap();

    // Check intersection details first.
    let isect = crate::intersect::find_intersections(&a, &b, &tol);
    match &isect {
        Ok(result) => {
            eprintln!(
                "[overlap] intersection: {} curves, {} hits, {} coincident",
                result.records.len(),
                result.edge_face_hits.len(),
                result.coincidence.face_overlaps.len()
            );
            for (i, rec) in result.records.iter().enumerate() {
                use popkern_core::geometry::CurveEval;
                let (t0, t1) = rec.curve.curve_3d.parameter_range();
                let p0 = rec.curve.curve_3d.point_at(t0);
                let p1 = rec.curve.curve_3d.point_at(t1);
                eprintln!(
                    "  curve {i}: a={:?} b={:?} tangent={} ({:.2},{:.2},{:.2})->({:.2},{:.2},{:.2})",
                    rec.face_a, rec.face_b, rec.is_tangent, p0.x, p0.y, p0.z, p1.x, p1.y, p1.z
                );
            }
            for ov in &result.coincidence.face_overlaps {
                eprintln!(
                    "  coincident: a={:?} b={:?} orient={:?} full={} center=({:.2},{:.2},{:.2})",
                    ov.face_a,
                    ov.face_b,
                    ov.orientation,
                    ov.full_overlap,
                    ov.overlap_center.x,
                    ov.overlap_center.y,
                    ov.overlap_center.z,
                );
            }
        }
        Err(e) => eprintln!("[overlap] intersection FAILED: {e:?}"),
    }

    // Run subtraction.
    match crate::subtraction(a, b, tol, None) {
        Ok(brep) => {
            let v = brep.vertices.len();
            let e = brep.edges.len();
            let f = brep.faces.len();
            let errors = brep.validate_all();
            eprintln!("[overlap] V={v} E={e} F={f} errors={}", errors.len());
            for err in &errors {
                eprintln!("[overlap]   {err}");
            }
            for (fh, face) in &brep.faces {
                let outer_count = {
                    let mut c = 0;
                    let start = brep.loops[face.outer_loop].half_edge;
                    let mut cur = start;
                    for _ in 0..200 {
                        c += 1;
                        cur = brep.half_edges[cur].next;
                        if cur == start {
                            break;
                        }
                    }
                    c
                };
                if outer_count == 5 || outer_count <= 3 {
                    // Trace vertices of short loops.
                    let start = brep.loops[face.outer_loop].half_edge;
                    let mut cur = start;
                    for _ in 0..10 {
                        let he = &brep.half_edges[cur];
                        let p = brep.vertices[he.origin].position;
                        let twin_he = &brep.half_edges[he.twin];
                        let pe = brep.vertices[twin_he.origin].position;
                        let has_curve = brep.edges[he.edge].curve.is_some();
                        eprintln!(
                            "  vtx ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) curve={has_curve}",
                            p.x, p.y, p.z, pe.x, pe.y, pe.z
                        );
                        cur = he.next;
                        if cur == start {
                            break;
                        }
                    }
                }
                eprintln!("[overlap] face {fh:?} outer={outer_count}");
            }
        }
        Err(e) => eprintln!("[overlap] subtraction FAILED: {e}"),
    }
}
