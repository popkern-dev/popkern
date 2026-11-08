use super::*;

#[test]
fn surface_complex_accepts_t012_components() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T012.obj"
    ))
    .unwrap();
    let complex =
        build_control_net_surface_complex(&net, PolyhedralCoverageOptions::default()).unwrap();
    assert_eq!(complex.candidate_splines, 15);
    assert_eq!(complex.components.len(), 15);
    assert!(complex.rejections.is_empty());
    let sewing = analyze_surface_complex_sewing(&complex, 1.0e-6);
    assert_eq!(sewing.components, 15);
    assert!(sewing.boundary_edges > 0);
}

/// PK-1119: T012's components are refined differently, so a coarse
/// boundary edge faces a finer chain whose middle vertices lie on its
/// curve. Split there, every such seam welds: the sheet is joined
/// through edges, never only through a vertex (unsplit, 32 pairs welded,
/// 68 boundary edges left and the 8-face part hung on vertices, which
/// re-read from STEP as separate parts).
#[test]
fn weld_splits_t012_boundary_t_junctions_so_no_part_hangs_on_a_vertex() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T012.obj"
    ))
    .unwrap();
    let options = PolyhedralCoverageOptions::default();
    let complex = build_control_net_surface_complex(&net, options).unwrap();
    let welded = weld_clean_surface_complex_boundaries(&complex, options).unwrap();
    // One piece, joined through edges: no part hangs on a vertex.
    assert_eq!(edge_joined_face_sets(&welded.brep), 1);
    assert_eq!(welded.report.t_junction_splits, 12);
    assert_eq!(welded.report.welded_pairs, 54);
    assert_eq!(welded.report.declined_pairs, 0);
    assert_eq!(welded.report.remaining_boundary_edges, 36);
    assert_eq!(welded.brep.piece_face_sets().len(), 1);
    assert!(welded.brep.validate().is_ok());
}

#[test]
fn surface_complex_preserves_partial_hand_rejections() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/hand.obj"
    ))
    .unwrap();
    let complex =
        build_control_net_surface_complex(&net, PolyhedralCoverageOptions::default()).unwrap();
    assert_eq!(complex.candidate_splines, 380);
    // 371/9 since PK-317: one Polar candidate's historical rejection was a
    // pole-degenerate normal artifact of the reference package's table
    // noise (see coverage_report_accepts_full_cylinder); the regenerated
    // exactly-symmetric polar masks accept it. Nine genuine rejections
    // remain and keep this path exercised.
    assert_eq!(complex.components.len(), 371);
    assert_eq!(complex.rejections.len(), 9);
    assert!(
        complex
            .rejections
            .iter()
            .all(|rejection| !rejection.reason.is_empty())
    );
}

#[test]
fn weld_clean_boundaries_reduces_tee_open_edges_transactionally() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/tee.obj"
    ))
    .unwrap();
    let options = PolyhedralCoverageOptions::default();
    let complex = build_control_net_surface_complex(&net, options).unwrap();
    let before = analyze_surface_complex_sewing(&complex, options.c0_tolerance);
    let welded = weld_clean_surface_complex_boundaries(&complex, options).unwrap();
    assert_eq!(welded.report.welded_pairs, before.matched_boundary_pairs);
    assert_eq!(
        welded.report.remaining_boundary_edges,
        before.unmatched_boundary_edges
    );
    assert!(welded.report.welded_pairs > 0);
    // Its interfaces meet vertex to vertex: no edge is split.
    assert_eq!(welded.report.t_junction_splits, 0);
    assert_eq!(welded.report.trim_repairs, 0);
    assert_eq!(welded.report.topology_errors, 0);
    // Since the regular-quad quarter split (PK-198) the tee's mixed
    // regular/star interfaces weld 1:1 and the cage closes fully —
    // the historical open edges and topology warnings were the
    // dual-cell/quarter-cell tiling mismatch, not real geometry.
    assert_eq!(welded.report.remaining_boundary_edges, 0);
    assert_eq!(welded.report.declined_pairs, 0);
    assert_eq!(welded.report.topology_warnings, 0);
    assert!(welded.report.export_ready);
    assert!(!welded.report.expected_open_sheet_boundaries);
    assert_eq!(welded.report.export_blockers, 0);
    assert_eq!(welded.report.missing_edge_curves, 0);
    assert_eq!(welded.report.trim_blockers, 0);
    assert!(welded.brep.validate().is_ok());
    let original_boundary_edges: usize = complex
        .components
        .iter()
        .map(|component| component.brep.free_edges().len())
        .sum();
    assert_eq!(original_boundary_edges, before.boundary_edges);
}

/// The face sets of `brep` joined through shared edges only (the
/// pieces a STEP import reads a shell as).
fn edge_joined_face_sets(brep: &crate::topology::BRep) -> usize {
    use std::collections::HashMap;
    let faces: Vec<_> = brep.faces.keys().collect();
    let index: HashMap<_, usize> = faces.iter().enumerate().map(|(i, &f)| (f, i)).collect();
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    let mut by_edge = HashMap::new();
    for he in brep.half_edges.values() {
        let Some(&i) = brep.loops.get(he.loop_).and_then(|l| index.get(&l.face)) else {
            continue;
        };
        let j = *by_edge.entry(he.edge).or_insert(i);
        let (a, b) = (root(&mut parent, i), root(&mut parent, j));
        parent[a] = b;
    }
    (0..faces.len())
        .map(|i| root(&mut parent, i))
        .collect::<std::collections::HashSet<_>>()
        .len()
}
