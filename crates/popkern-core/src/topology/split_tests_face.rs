#![cfg(test)]
use super::*;

#[test]
fn split_face_with_split_edges() {
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

    let f_before = brep.faces.len();
    let split1 = brep.split_edge(mev1.edge, 0.5).unwrap();
    let split2 = brep.split_edge(mev3.edge, 0.5).unwrap();

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

#[test]
fn split_face_reversed_endpoint_order() {
    let mut brep = BRep::new(Tolerance::default());
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

    let split1 = brep.split_edge(mev1.edge, 0.5).unwrap();
    let split2 = brep.split_edge(mev3.edge, 0.5).unwrap();
    let p1 = brep.vertices[split1.vertex].position;
    let p2 = brep.vertices[split2.vertex].position;

    let trim_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p2, p1],
        weights: vec![1.0, 1.0],
    };
    let trim_uv = trim_3d.clone();

    let result = brep
        .split_face(mef_r.face, &trim_uv, &trim_3d, &Tolerance::default())
        .unwrap();
    assert!(brep.faces.contains_key(result.face_a));
    assert!(brep.faces.contains_key(result.face_b));
}

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

    let trim = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(50.0, 50.0, 0.0), Point3::new(60.0, 60.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let err = brep.split_face(mef_r.face, &trim, &trim, &Tolerance::default());
    assert!(err.is_err());
}

#[test]
fn split_face_exit_found_entry_second_pass() {
    let mut brep = BRep::new(Tolerance::default());
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

    let split1 = brep.split_edge(mev1.edge, 0.5).unwrap();
    let split2 = brep.split_edge(mev3.edge, 0.5).unwrap();
    let p1 = brep.vertices[split1.vertex].position;
    let p2 = brep.vertices[split2.vertex].position;

    let trim_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p2, p1],
        weights: vec![1.0, 1.0],
    };
    let trim_uv = trim_3d.clone();

    let result = brep.split_face(mef_r.face, &trim_uv, &trim_3d, &Tolerance::default());
    assert!(result.is_ok());
}

#[test]
fn split_edge_with_curve_geometry() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge_h = brep
        .edges
        .keys()
        .find(|&eh| brep.edges[eh].curve.is_some())
        .unwrap();

    let result = brep.split_edge(edge_h, 0.5).unwrap();
    assert!(brep.edges[result.edge_a].curve.is_some());
    assert!(brep.edges[result.edge_b].curve.is_some());
}

#[test]
fn split_edge_preserves_pcurves() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let he = cube.edges[edge].half_edge;
    assert!(cube.half_edges[he].pcurve.is_some());

    let result = cube.split_edge(edge, 0.5).unwrap();
    let he_a = cube.edges[result.edge_a].half_edge;
    let he_b = cube.edges[result.edge_b].half_edge;
    assert!(cube.half_edges[he_a].pcurve.is_some());
    assert!(cube.half_edges[he_b].pcurve.is_some());

    let twin_a = cube.half_edges[he_a].twin;
    let twin_b = cube.half_edges[he_b].twin;
    assert!(cube.half_edges[twin_a].pcurve.is_some());
    assert!(cube.half_edges[twin_b].pcurve.is_some());
}

#[test]
fn split_edge_pcurve_endpoints_consistent() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let he = cube.edges[edge].half_edge;
    let original_pcurve = cube.half_edges[he].pcurve.clone().unwrap();
    let orig_range = original_pcurve.parameter_range();
    let orig_start = original_pcurve.point_at(orig_range.0);
    let orig_end = original_pcurve.point_at(orig_range.1);

    let result = cube.split_edge(edge, 0.5).unwrap();

    let he_a = cube.edges[result.edge_a].half_edge;
    let pcurve_a = cube.half_edges[he_a].pcurve.as_ref().unwrap();
    let range_a = pcurve_a.parameter_range();
    let start_a = pcurve_a.point_at(range_a.0);
    assert!((start_a.x - orig_start.x).abs() < 1e-10);
    assert!((start_a.y - orig_start.y).abs() < 1e-10);

    let he_b = cube.edges[result.edge_b].half_edge;
    let pcurve_b = cube.half_edges[he_b].pcurve.as_ref().unwrap();
    let range_b = pcurve_b.parameter_range();
    let end_b = pcurve_b.point_at(range_b.1);
    assert!((end_b.x - orig_end.x).abs() < 1e-10);
    assert!((end_b.y - orig_end.y).abs() < 1e-10);

    let end_a = pcurve_a.point_at(range_a.1);
    let start_b = pcurve_b.point_at(range_b.0);
    assert!((end_a.x - start_b.x).abs() < 1e-10);
    assert!((end_a.y - start_b.y).abs() < 1e-10);
}

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

    let new_edge = result.new_edges[0];
    let he_new = brep.edges[new_edge].half_edge;
    let he_twin = brep.half_edges[he_new].twin;
    assert!(brep.half_edges[he_new].pcurve.is_some());
    assert!(brep.half_edges[he_twin].pcurve.is_some());
}

#[test]
fn split_edge_no_pcurve_stays_none() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    let he = brep.edges[mev1.edge].half_edge;
    assert!(brep.half_edges[he].pcurve.is_none());

    let result = brep.split_edge(mev1.edge, 0.5).unwrap();
    let he_a = brep.edges[result.edge_a].half_edge;
    let he_b = brep.edges[result.edge_b].half_edge;
    assert!(brep.half_edges[he_a].pcurve.is_none());
    assert!(brep.half_edges[he_b].pcurve.is_none());
}

/// PK-428: the inner-loop matcher (`match_inner_loop_cut`) may return
/// entry/exit SWAPPED relative to the trim's parameterization — the
/// winding rule owns the order. The trim pcurves must then follow each
/// half-edge's OWN traversal direction; the historical
/// forward-on-`half_edge` assignment put a forward pcurve on a
/// backwards half-edge for the swapped order, and every consumer that
/// samples a loop through its pcurves read geometry running against
/// the traversal. Both trim parameterizations are exercised so
/// whichever order the winding rule picks, the swapped branch is hit
/// by one of them.
#[test]
fn split_face_inner_loop_cut_pcurves_follow_half_edge_direction() {
    use crate::geometry::SurfaceOnFace;
    use crate::geometry::surface::SurfaceGeometry;
    use crate::topology::FaceHandle;
    use crate::topology::trim::InnerLoopPiece;

    let tolerance = Tolerance::default();
    let line = |a: Point3, b: Point3| NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    };
    let polyline = |pts: &[Point3]| {
        let n = pts.len();
        let mut knots = vec![0.0, 0.0];
        for i in 1..n - 1 {
            knots.push(i as f64 / (n - 1) as f64);
        }
        knots.push(1.0);
        knots.push(1.0);
        NurbsCurve {
            degree: 1,
            knots,
            control_points: pts.to_vec(),
            weights: vec![1.0; n],
        }
    };
    // Square face in the z=0 plane with an identity chart (uv = xy),
    // outer loop CCW, plus a square hole (CW) built from four pieces
    // so its corners are real loop vertices.
    let build = || -> (BRep, FaceHandle) {
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
        let face = mef_r.face;
        brep.faces[face].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: crate::geometry::Vec3::new(0.0, 0.0, 1.0),
                u_axis: crate::geometry::Vec3::new(1.0, 0.0, 0.0),
                v_axis: crate::geometry::Vec3::new(0.0, 1.0, 0.0),
            },
            orientation: true,
        });
        // Authored so the STORED face-side loop winds like a real
        // boolean hole (opposite the outer loop) — the winding rule
        // only swaps entry/exit against that storage direction.
        let hole = [
            Point3::new(0.3, 0.3, 0.0),
            Point3::new(0.7, 0.3, 0.0),
            Point3::new(0.7, 0.7, 0.0),
            Point3::new(0.3, 0.7, 0.0),
        ];
        let pieces: Vec<InnerLoopPiece> = (0..4)
            .map(|i| {
                let c = line(hole[i], hole[(i + 1) % 4]);
                InnerLoopPiece {
                    curve_3d: c.clone(),
                    curve_uv: c,
                }
            })
            .collect();
        brep.add_inner_loop_from_pieces(face, &pieces, &tolerance)
            .unwrap();
        (brep, face)
    };
    // A V-shaped trim in the material below the hole, both ends on
    // hole corners.
    let v_path = [
        Point3::new(0.3, 0.3, 0.0),
        Point3::new(0.5, 0.1, 0.0),
        Point3::new(0.7, 0.3, 0.0),
    ];
    for reversed in [false, true] {
        let mut pts = v_path.to_vec();
        if reversed {
            pts.reverse();
        }
        let trim = polyline(&pts);
        let (mut brep, face) = build();
        let result = brep
            .split_face(face, &trim, &trim, &tolerance)
            .unwrap_or_else(|e| panic!("inner-loop cut (reversed={reversed}): {e}"));
        let edge = result.new_edges[0];
        let he = brep.edges[edge].half_edge;
        let twin = brep.half_edges[he].twin;
        for h in [he, twin] {
            let origin = brep.vertices[brep.half_edges[h].origin].position;
            let dest = brep.vertices[brep.half_edges[brep.half_edges[h].twin].origin].position;
            let pcurve = brep.half_edges[h].pcurve.as_ref().expect("trim pcurve");
            let (t0, t1) = pcurve.parameter_range();
            let s = pcurve.point_at(t0);
            let e = pcurve.point_at(t1);
            let d_fwd = ((s.x - origin.x).powi(2) + (s.y - origin.y).powi(2)).sqrt()
                + ((e.x - dest.x).powi(2) + (e.y - dest.y).powi(2)).sqrt();
            assert!(
                d_fwd < 1.0e-9,
                "pcurve runs against its half-edge (reversed={reversed}): \
                 pcurve ({:.3},{:.3})->({:.3},{:.3}) vs he ({:.3},{:.3})->({:.3},{:.3})",
                s.x,
                s.y,
                e.x,
                e.y,
                origin.x,
                origin.y,
                dest.x,
                dest.y
            );
        }
    }
}
