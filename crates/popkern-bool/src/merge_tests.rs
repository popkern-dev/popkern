//! Tests for Stage 5 merge.
#[cfg(test)]
use std::collections::HashSet;

use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

use crate::merge::merge_faces;
use crate::select::SelectionResult;
use crate::split::{FragmentMap, SplitResult};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Build a SplitResult with no actual splitting — both operands are
/// passed through unchanged with identity fragment maps.
fn identity_split(a: BRep, b: BRep) -> SplitResult {
    let fragments_a = FragmentMap {
        origin: a.faces.keys().map(|fh| (fh, fh)).collect(),
        intersection_edges: HashSet::new(),
    };
    let fragments_b = FragmentMap {
        origin: b.faces.keys().map(|fh| (fh, fh)).collect(),
        intersection_edges: HashSet::new(),
    };
    SplitResult {
        a,
        b,
        fragments_a,
        fragments_b,
    }
}

/// Build a SelectionResult that keeps all faces from both operands.
fn select_all(split: &SplitResult) -> SelectionResult {
    let mut keep = HashSet::new();
    for fh in split.a.faces.keys() {
        keep.insert(fh);
    }
    for fh in split.b.faces.keys() {
        keep.insert(fh);
    }
    SelectionResult {
        keep,
        flip: HashSet::new(),
    }
}

/// Build a SelectionResult that keeps only faces from operand A.
fn select_a_only(split: &SplitResult) -> SelectionResult {
    let keep: HashSet<_> = split.a.faces.keys().collect();
    SelectionResult {
        keep,
        flip: HashSet::new(),
    }
}

#[test]
fn merge_single_cube_copies_all_topology() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    let empty = BRep::new(tol());

    let split = identity_split(cube, empty);
    let selection = select_a_only(&split);
    let result = merge_faces(&split, &selection, &tol()).expect("merge");

    // A unit cube has 8 vertices, 12 edges, 6 faces.
    assert_eq!(result.vertices.len(), 8, "vertices");
    assert_eq!(result.faces.len(), 6, "faces");
    assert_eq!(result.edges.len(), 12, "edges");
    assert_eq!(result.half_edges.len(), 24, "half-edges");
}

#[test]
fn merge_single_cube_stitches_all_edges() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    let empty = BRep::new(tol());

    let split = identity_split(cube, empty);
    let selection = select_a_only(&split);
    let result = merge_faces(&split, &selection, &tol()).expect("merge");

    // After stitching, no half-edge should have self-twin.
    let unpaired: Vec<_> = result
        .half_edges
        .iter()
        .filter(|(heh, he)| he.twin == *heh)
        .collect();
    assert!(
        unpaired.is_empty(),
        "found {} unpaired half-edges after merge",
        unpaired.len()
    );
}

#[test]
fn merge_two_disjoint_cubes() {
    let cube_a = BRep::unit_cube(tol()).expect("cube a");
    let mut cube_b = BRep::unit_cube(tol()).expect("cube b");

    // Move cube_b far away so no vertex dedup occurs.
    for (_, v) in &mut cube_b.vertices {
        v.position = v.position + Vec3::new(10.0, 0.0, 0.0);
    }

    let split = identity_split(cube_a, cube_b);
    let selection = select_all(&split);
    let result = merge_faces(&split, &selection, &tol()).expect("merge");

    // Two cubes: 16 vertices, 24 edges, 12 faces.
    assert_eq!(result.vertices.len(), 16, "vertices");
    assert_eq!(result.faces.len(), 12, "faces");
    assert_eq!(result.edges.len(), 24, "edges");
}

#[test]
fn merge_empty_selection() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    let empty = BRep::new(tol());

    let split = identity_split(cube, empty);
    let selection = SelectionResult {
        keep: HashSet::new(),
        flip: HashSet::new(),
    };
    let result = merge_faces(&split, &selection, &tol()).expect("merge");

    assert_eq!(result.faces.len(), 0);
    assert_eq!(result.vertices.len(), 0);
}

#[test]
fn merge_with_flip_reverses_surface_orientation() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    let empty = BRep::new(tol());

    let split = identity_split(cube, empty);

    // Select all faces from A and flip them all.
    let keep: HashSet<_> = split.a.faces.keys().collect();
    let flip = keep.clone();
    let selection = SelectionResult { keep, flip };

    let result = merge_faces(&split, &selection, &tol()).expect("merge");

    // All faces should have flipped orientation.
    for (_, face) in &result.faces {
        if let Some(ref surf) = face.surface {
            // Original cube surfaces have orientation = true.
            // After flip, they should be false.
            assert!(!surf.orientation, "expected flipped orientation");
        }
    }
}

#[test]
fn merge_deduplicates_coincident_vertices() {
    // Two cubes sharing a face — vertices on the shared face should
    // be deduplicated.
    let cube_a = BRep::unit_cube(tol()).expect("cube a");
    let mut cube_b = BRep::unit_cube(tol()).expect("cube b");

    // Move cube_b by 1 in x so they share the face at x=1.
    for (_, v) in &mut cube_b.vertices {
        v.position = v.position + Vec3::new(1.0, 0.0, 0.0);
    }

    let split = identity_split(cube_a, cube_b);
    let selection = select_all(&split);
    let result = merge_faces(&split, &selection, &tol()).expect("merge");

    // Two cubes sharing 4 vertices: 8 + 8 - 4 = 12 unique vertices.
    assert_eq!(result.vertices.len(), 12, "vertices should be deduplicated");
}

#[test]
fn merge_result_has_solid_and_shell() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    let empty = BRep::new(tol());

    let split = identity_split(cube, empty);
    let selection = select_a_only(&split);
    let result = merge_faces(&split, &selection, &tol()).expect("merge");

    assert_eq!(result.solids.len(), 1);
    assert_eq!(result.shells.len(), 1);

    let solid_h = result.solids.keys().next().expect("solid handle");
    let shell_h = result.solids[solid_h].shell;
    assert!(result.shells.contains_key(shell_h));
    assert_eq!(result.shells[shell_h].solid, solid_h);
}

/// Build a SelectionResult that keeps only faces from operand B.
fn select_b_only(split: &SplitResult) -> SelectionResult {
    let keep: HashSet<_> = split.b.faces.keys().collect();
    SelectionResult {
        keep,
        flip: HashSet::new(),
    }
}

#[test]
fn merge_b_faces_sets_first_face_from_b() {
    // Keep only B's faces — this exercises the first_face from B branch
    // in merge_faces (lines that set first_face when it was None after
    // processing A's empty selection).
    let empty = BRep::new(tol());
    let cube = BRep::unit_cube(tol()).expect("cube");

    let split = identity_split(empty, cube);
    let selection = select_b_only(&split);
    let result = merge_faces(&split, &selection, &tol()).expect("merge");

    // B is the unit cube — 6 faces.
    assert_eq!(result.faces.len(), 6);
    // Shell should point at a real face.
    let shell_h = result.shells.keys().next().expect("shell");
    let face_h = result.shells[shell_h].face;
    assert!(
        result.faces.contains_key(face_h),
        "shell face should be valid"
    );
}

/// Build a BRep with one face that has an inner loop, using kemr.
///
/// Topology (spur tree): v0 -bridge- v1 -spur- v2; v0 -spur- v3
/// After kemr on the bridge edge, the face has outer loop [v0->v3->v0]
/// and inner loop [v1->v2->v1].
fn brep_with_inner_loop() -> BRep {
    use popkern_core::geometry::Point3;
    let mut brep = BRep::new(tol());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let _mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(2.0, 0.0, 0.0))
        .expect("mev2");
    let _mev3 = brep
        .mev(mev1.half_edge, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    // Split into outer loop + inner loop.
    brep.kemr(mev1.half_edge).expect("kemr");
    brep
}

#[test]
fn merge_face_with_inner_loop_copies_inner_loops() {
    // A face with an inner loop should have its inner loop copied into the result.
    let brep_a = brep_with_inner_loop();
    let brep_b = BRep::new(tol());

    assert_eq!(brep_a.loops.len(), 2, "should have outer + inner loop");
    // Find the face that has an inner loop.
    let has_inner = brep_a.faces.values().any(|f| !f.inner_loops.is_empty());
    assert!(has_inner, "brep should have a face with an inner loop");

    let split = identity_split(brep_a, brep_b);
    let selection = select_a_only(&split);
    let result = merge_faces(&split, &selection, &tol()).expect("merge");

    // The result should also have a face with an inner loop.
    let has_inner_in_result = result.faces.values().any(|f| !f.inner_loops.is_empty());
    assert!(has_inner_in_result, "result should preserve the inner loop");
    // Should have 2 loops: outer + inner.
    assert_eq!(result.loops.len(), 2, "outer + inner loop should be copied");
}

#[test]
fn merge_loop_half_edges_form_cycle() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    let empty = BRep::new(tol());

    let split = identity_split(cube, empty);
    let selection = select_a_only(&split);
    let result = merge_faces(&split, &selection, &tol()).expect("merge");

    // Verify every loop's half-edges form a proper cycle.
    for (loop_h, loop_) in &result.loops {
        let start = loop_.half_edge;
        let mut current = start;
        let mut count = 0;
        loop {
            let he = &result.half_edges[current];
            // Check prev/next consistency.
            assert_eq!(
                result.half_edges[he.next].prev, current,
                "next->prev should be self"
            );
            assert_eq!(he.loop_, loop_h, "half-edge should belong to this loop");
            current = he.next;
            count += 1;
            if current == start {
                break;
            }
            assert!(count < 1000, "infinite loop detected");
        }
        assert!(count >= 3, "loop should have at least 3 half-edges");
    }
}
