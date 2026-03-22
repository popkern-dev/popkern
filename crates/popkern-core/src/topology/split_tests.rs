#[cfg(test)]
#[allow(unused_imports)]
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::Point3;
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
    assert!(brep
        .split_face(mef_r.face, &trim, &trim, &Tolerance::default())
        .is_err());
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
