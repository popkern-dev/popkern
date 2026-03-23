use crate::geometry::Point3;
#[cfg(test)]
#[allow(unused_imports)]
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

// ==========================================================================
// split_edge tests
// ==========================================================================

#[test]
fn split_edge_simple() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    let result = brep.split_edge(mev1.edge, 0.5).unwrap();

    assert_eq!(brep.vertices.len(), 3);
    assert_eq!(brep.edges.len(), 2);
    assert_eq!(brep.half_edges.len(), 4);
    let v = &brep.vertices[result.vertex];
    assert!((v.position.x - 0.5).abs() < 1e-10);
}

#[test]
fn split_edge_on_cube() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let v_before = brep.vertices.len();
    let e_before = brep.edges.len();
    let edge = brep.edges.keys().next().unwrap();

    brep.split_edge(edge, 0.5).unwrap();

    assert_eq!(brep.vertices.len(), v_before + 1);
    assert_eq!(brep.edges.len(), e_before + 1);
}

#[test]
fn split_edge_invalid_handle() {
    let mut brep = BRep::new(Tolerance::default());
    assert!(brep.split_edge(Default::default(), 0.5).is_err());
}

#[test]
fn split_edge_half_edge_chain_valid() {
    // After split, walking the loop should still form a cycle
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let _mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();

    brep.split_edge(mev1.edge, 0.5).unwrap();

    // Walk the loop -- should return to start
    let start = brep.loops[mvfs.loop_].half_edge;
    let mut he = start;
    let mut count = 0;
    loop {
        he = brep.half_edges[he].next;
        count += 1;
        if he == start || count > 100 {
            break;
        }
    }
    assert_eq!(he, start, "Loop should be closed after split_edge");
}

#[test]
fn split_edge_t_out_of_range() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    assert!(brep.split_edge(mev1.edge, 0.0).is_err());
    assert!(brep.split_edge(mev1.edge, 1.0).is_err());
    assert!(brep.split_edge(mev1.edge, -0.1).is_err());
    assert!(brep.split_edge(mev1.edge, 1.1).is_err());
}

// ==========================================================================
// split_face tests
// ==========================================================================

#[test]
fn split_face_with_split_edges() {
    let mut brep = BRep::new(Tolerance::default());
    // Build a square face: v0(0,0) -> v1(1,0) -> v2(1,1) -> v3(0,1)
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

    let f_before = brep.faces.len();

    // Split edge v0->v1 at midpoint -> vertex at (0.5, 0, 0)
    let edge_01 = mev1.edge;
    let split1 = brep.split_edge(edge_01, 0.5).unwrap();

    // Split edge v2->v3 at midpoint -> vertex at (0.5, 1, 0)
    let edge_23 = mev3.edge;
    let split2 = brep.split_edge(edge_23, 0.5).unwrap();

    // Trim curve connecting the two new vertices
    let p1 = brep.vertices[split1.vertex].position;
    let p2 = brep.vertices[split2.vertex].position;
    let trim_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p1, p2],
        weights: vec![1.0, 1.0],
    };
    let trim_uv = trim_3d.clone();

    let result = brep
        .split_face(mef_r.face, &trim_uv, &trim_3d, &Tolerance::default())
        .unwrap();

    assert_eq!(brep.faces.len(), f_before + 1);
    assert!(brep.faces.contains_key(result.face_a));
    assert!(brep.faces.contains_key(result.face_b));
    assert!(!result.new_edges.is_empty());
}

#[test]
fn split_face_vertex_not_found_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

    // Trim curve with endpoints NOT on the face boundary
    let trim = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(5.0, 5.0, 0.0), Point3::new(6.0, 6.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    assert!(
        brep.split_face(mef_r.face, &trim, &trim, &Tolerance::default())
            .is_err()
    );
}

#[test]
fn split_face_preserves_surface() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let had_surface = brep.faces[face].surface.is_some();

    // Split an edge on this face, then split_face
    let loop_handle = brep.faces[face].outer_loop;
    let he1 = brep.loops[loop_handle].half_edge;
    let he2 = brep.half_edges[he1].next;
    let edge1 = brep.half_edges[he1].edge;
    let edge2 = brep.half_edges[he2].edge;

    let s1 = brep.split_edge(edge1, 0.5).unwrap();
    let s2 = brep.split_edge(edge2, 0.5).unwrap();

    let p1 = brep.vertices[s1.vertex].position;
    let p2 = brep.vertices[s2.vertex].position;
    let trim = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p1, p2],
        weights: vec![1.0, 1.0],
    };

    let result = brep
        .split_face(face, &trim, &trim, &Tolerance::default())
        .unwrap();

    // Both faces should have surface geometry if original did
    if had_surface {
        assert!(brep.faces[result.face_a].surface.is_some());
        assert!(brep.faces[result.face_b].surface.is_some());
    }
}

#[test]
fn split_face_invalid_face_handle() {
    let mut brep = BRep::new(Tolerance::default());
    let trim = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    assert!(
        brep.split_face(Default::default(), &trim, &trim, &Tolerance::default())
            .is_err()
    );
}

/// Test the second-pass loop walk where entry is found but exit is not on first pass.
/// This exercises lines 274-292 (the he_entry.is_some() && he_exit.is_none() branch).
#[test]
fn split_face_reversed_endpoint_order() {
    let mut brep = BRep::new(Tolerance::default());
    // Build a square face: v0(0,0) -> v1(2,0) -> v2(2,2) -> v3(0,2)
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(2.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(2.0, 2.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 2.0, 0.0))
        .unwrap();
    let mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

    // Split two opposite edges to get midpoints
    let split1 = brep.split_edge(mev1.edge, 0.5).unwrap(); // midpoint on v0->v1 at (1,0,0)
    let split2 = brep.split_edge(mev3.edge, 0.5).unwrap(); // midpoint on v2->v3 at (1,2,0)

    let p1 = brep.vertices[split1.vertex].position;
    let p2 = brep.vertices[split2.vertex].position;

    // Create trim curve in REVERSE order (end->start) relative to the
    // loop walk. The first pass finds entry but not exit (or vice-versa),
    // triggering the second-pass code.
    let trim_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p2, p1], // reversed from the loop walk order
        weights: vec![1.0, 1.0],
    };
    let trim_uv = trim_3d.clone();

    let result = brep
        .split_face(mef_r.face, &trim_uv, &trim_3d, &Tolerance::default())
        .unwrap();
    assert!(brep.faces.contains_key(result.face_a));
    assert!(brep.faces.contains_key(result.face_b));
}

/// Test where neither entry nor exit found on first pass.
#[test]
fn split_face_no_entry_no_exit_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

    // Trim curve endpoints far from any vertex
    let trim = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(50.0, 50.0, 0.0), Point3::new(60.0, 60.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let err = brep.split_face(mef_r.face, &trim, &trim, &Tolerance::default());
    assert!(err.is_err());
}

/// Test where exit is found on first pass but entry is not (lines 293-311).
#[test]
fn split_face_exit_found_entry_second_pass() {
    let mut brep = BRep::new(Tolerance::default());
    // Build a square: v0(0,0) -> v1(2,0) -> v2(2,2) -> v3(0,2)
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(2.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(2.0, 2.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 2.0, 0.0))
        .unwrap();
    let mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

    // Split edges to get midpoints
    let split1 = brep.split_edge(mev1.edge, 0.5).unwrap(); // (1,0,0)
    let split2 = brep.split_edge(mev3.edge, 0.5).unwrap(); // (1,2,0)

    let p1 = brep.vertices[split1.vertex].position;
    let p2 = brep.vertices[split2.vertex].position;

    // Construct trim curve where the first point encountered in the loop is the
    // "end" point, so first pass finds exit but not entry.
    // The loop visits: v0 -> split1_vertex -> v1 -> v2 -> split2_vertex -> v3
    // If curve_start = p2 (split2) and curve_end = p1 (split1),
    // and the first vertex matched in the loop walk is split1 (matching curve_end),
    // then first pass finds exit but not entry, triggering second pass for entry.
    let trim_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p2, p1],
        weights: vec![1.0, 1.0],
    };
    let trim_uv = trim_3d.clone();

    // This should work (second pass finds entry)
    let result = brep.split_face(mef_r.face, &trim_uv, &trim_3d, &Tolerance::default());
    assert!(result.is_ok());
}

/// Test split_edge with curve geometry (exercises the NURBS split path).
#[test]
fn split_edge_with_curve_geometry() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();

    // Find an edge with curve geometry
    let edge_h = brep
        .edges
        .keys()
        .find(|&eh| brep.edges[eh].curve.is_some())
        .unwrap();

    let result = brep.split_edge(edge_h, 0.5).unwrap();

    // Both edges should have curve geometry after split
    assert!(brep.edges[result.edge_a].curve.is_some());
    assert!(brep.edges[result.edge_b].curve.is_some());
}

/// Test that split_edge preserves pcurves on all half-edges.
#[test]
fn split_edge_preserves_pcurves() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let he = cube.edges[edge].half_edge;
    assert!(
        cube.half_edges[he].pcurve.is_some(),
        "cube half-edge should have pcurve"
    );

    let result = cube.split_edge(edge, 0.5).unwrap();

    // Both resulting edges' primary half-edges should have pcurves.
    let he_a = cube.edges[result.edge_a].half_edge;
    let he_b = cube.edges[result.edge_b].half_edge;
    assert!(
        cube.half_edges[he_a].pcurve.is_some(),
        "edge_a half-edge missing pcurve"
    );
    assert!(
        cube.half_edges[he_b].pcurve.is_some(),
        "edge_b half-edge missing pcurve"
    );

    // Twins should also have pcurves.
    let twin_a = cube.half_edges[he_a].twin;
    let twin_b = cube.half_edges[he_b].twin;
    assert!(
        cube.half_edges[twin_a].pcurve.is_some(),
        "edge_a twin missing pcurve"
    );
    assert!(
        cube.half_edges[twin_b].pcurve.is_some(),
        "edge_b twin missing pcurve"
    );
}

/// Test that split_edge pcurve endpoints are geometrically consistent.
#[test]
fn split_edge_pcurve_endpoints_consistent() {
    use crate::geometry::curve::CurveEval;

    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let he = cube.edges[edge].half_edge;
    let original_pcurve = cube.half_edges[he].pcurve.clone().unwrap();
    let orig_range = original_pcurve.parameter_range();
    let orig_start = original_pcurve.point_at(orig_range.0);
    let orig_end = original_pcurve.point_at(orig_range.1);

    let result = cube.split_edge(edge, 0.5).unwrap();

    // he_a pcurve should start at original start
    let he_a = cube.edges[result.edge_a].half_edge;
    let pcurve_a = cube.half_edges[he_a].pcurve.as_ref().unwrap();
    let range_a = pcurve_a.parameter_range();
    let start_a = pcurve_a.point_at(range_a.0);
    assert!(
        (start_a.x - orig_start.x).abs() < 1e-10 && (start_a.y - orig_start.y).abs() < 1e-10,
        "pcurve_a start should match original start"
    );

    // he_b pcurve should end at original end
    let he_b = cube.edges[result.edge_b].half_edge;
    let pcurve_b = cube.half_edges[he_b].pcurve.as_ref().unwrap();
    let range_b = pcurve_b.parameter_range();
    let end_b = pcurve_b.point_at(range_b.1);
    assert!(
        (end_b.x - orig_end.x).abs() < 1e-10 && (end_b.y - orig_end.y).abs() < 1e-10,
        "pcurve_b end should match original end"
    );

    // The endpoint of pcurve_a should match the start of pcurve_b (continuity at split)
    let end_a = pcurve_a.point_at(range_a.1);
    let start_b = pcurve_b.point_at(range_b.0);
    assert!(
        (end_a.x - start_b.x).abs() < 1e-10 && (end_a.y - start_b.y).abs() < 1e-10,
        "pcurve_a end should match pcurve_b start at split point"
    );
}

/// Test that split_face assigns pcurves to the new splitting edge.
#[test]
fn split_face_assigns_pcurves_to_new_edge() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = brep.faces.keys().next().unwrap();

    let loop_handle = brep.faces[face].outer_loop;
    let he1 = brep.loops[loop_handle].half_edge;
    let he2 = brep.half_edges[he1].next;
    let edge1 = brep.half_edges[he1].edge;
    let edge2 = brep.half_edges[he2].edge;

    let s1 = brep.split_edge(edge1, 0.5).unwrap();
    let s2 = brep.split_edge(edge2, 0.5).unwrap();

    let p1 = brep.vertices[s1.vertex].position;
    let p2 = brep.vertices[s2.vertex].position;
    let trim = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p1, p2],
        weights: vec![1.0, 1.0],
    };

    let result = brep
        .split_face(face, &trim, &trim, &Tolerance::default())
        .unwrap();

    // The new edge's half-edges should both have pcurves.
    let new_edge = result.new_edges[0];
    let he_new = brep.edges[new_edge].half_edge;
    let he_twin = brep.half_edges[he_new].twin;
    assert!(
        brep.half_edges[he_new].pcurve.is_some(),
        "new edge half-edge should have pcurve"
    );
    assert!(
        brep.half_edges[he_twin].pcurve.is_some(),
        "new edge twin should have pcurve"
    );
}

/// Test that split_edge without pcurves leaves None (no crash).
#[test]
fn split_edge_no_pcurve_stays_none() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    // These half-edges have no pcurves
    let he = brep.edges[mev1.edge].half_edge;
    assert!(brep.half_edges[he].pcurve.is_none());

    let result = brep.split_edge(mev1.edge, 0.5).unwrap();

    let he_a = brep.edges[result.edge_a].half_edge;
    let he_b = brep.edges[result.edge_b].half_edge;
    assert!(brep.half_edges[he_a].pcurve.is_none());
    assert!(brep.half_edges[he_b].pcurve.is_none());
}
