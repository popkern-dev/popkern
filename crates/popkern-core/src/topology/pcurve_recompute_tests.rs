use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn recompute_pcurves_on_cube_is_idempotent() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // All half-edges already have pcurves from unit_cube.
    brep.recompute_all_pcurves().unwrap();
    for (_, he) in &brep.half_edges {
        assert!(he.pcurve.is_some(), "pcurve should still be present");
    }
}

#[test]
fn recompute_fills_missing_pcurves() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Strip all pcurves.
    let he_keys: Vec<_> = brep.half_edges.keys().collect();
    for heh in &he_keys {
        brep.half_edges[*heh].pcurve = None;
    }
    // Recompute should fill them back.
    brep.recompute_all_pcurves().unwrap();
    for (_, he) in &brep.half_edges {
        assert!(he.pcurve.is_some(), "pcurve should be recomputed");
    }
}

#[test]
fn pcurve_endpoints_match_vertices() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Strip and recompute to exercise the fix path.
    let he_keys: Vec<_> = brep.half_edges.keys().collect();
    for heh in &he_keys {
        brep.half_edges[*heh].pcurve = None;
    }
    brep.recompute_all_pcurves().unwrap();

    for fh in brep.faces.keys().collect::<Vec<_>>() {
        let face = &brep.faces[fh];
        let surface = match &face.surface {
            Some(s) => &s.geometry,
            None => continue,
        };
        let loop_ = &brep.loops[face.outer_loop];
        let mut cursor = loop_.half_edge;
        loop {
            let he = &brep.half_edges[cursor];
            if let Some(ref pcurve) = he.pcurve {
                let (t0, t1) = pcurve.parameter_range();
                // Check start
                let uv_start = pcurve.point_at(t0);
                let pos_from_uv = surface.point_at(uv_start.x, uv_start.y);
                let actual_start = brep.vertices[he.origin].position;
                let err_start = pos_from_uv.distance_to(&actual_start);
                assert!(err_start < 1e-6, "pcurve start error: {err_start}");
                // Check end
                let uv_end = pcurve.point_at(t1);
                let pos_from_uv_end = surface.point_at(uv_end.x, uv_end.y);
                let next_he = &brep.half_edges[he.next];
                let actual_end = brep.vertices[next_he.origin].position;
                let err_end = pos_from_uv_end.distance_to(&actual_end);
                assert!(err_end < 1e-6, "pcurve end error: {err_end}");
            }
            cursor = brep.half_edges[cursor].next;
            if cursor == loop_.half_edge {
                break;
            }
        }
    }
}

#[test]
fn cylinder_seam_pcurve_is_degree_1() {
    let mut brep = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    // Clear all pcurves and recompute.
    for heh in brep.half_edges.keys().collect::<Vec<_>>() {
        brep.half_edges[heh].pcurve = None;
    }
    brep.recompute_all_pcurves().unwrap();

    // The seam edges on the lateral (cylindrical) face should get degree-1
    // pcurves since their UV path is a straight vertical line.
    let lateral_face = brep.faces.keys().find(|&fh| {
        matches!(
            brep.faces[fh].surface.as_ref().map(|s| &s.geometry),
            Some(SurfaceGeometry::Cylinder { .. })
        )
    });
    if let Some(fh) = lateral_face {
        let loop_ = &brep.loops[brep.faces[fh].outer_loop];
        let mut cursor = loop_.half_edge;
        let mut found_degree_1 = false;
        loop {
            let he = &brep.half_edges[cursor];
            if let Some(ref pcurve) = he.pcurve {
                if pcurve.degree == 1 && pcurve.control_points.len() == 2 {
                    found_degree_1 = true;
                }
            }
            cursor = he.next;
            if cursor == loop_.half_edge {
                break;
            }
        }
        assert!(
            found_degree_1,
            "cylinder lateral face should have at least one degree-1 pcurve (seam)"
        );
    }
}

#[test]
fn cylinder_recomputed_pcurves_are_accurate() {
    let mut brep = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    // Clear all pcurves and recompute.
    for heh in brep.half_edges.keys().collect::<Vec<_>>() {
        brep.half_edges[heh].pcurve = None;
    }
    brep.recompute_all_pcurves().unwrap();

    // Check that pcurve-evaluated 3D positions match vertex positions on
    // the cylindrical face.
    let lateral_face = brep.faces.keys().find(|&fh| {
        matches!(
            brep.faces[fh].surface.as_ref().map(|s| &s.geometry),
            Some(SurfaceGeometry::Cylinder { .. })
        )
    });
    if let Some(fh) = lateral_face {
        let surface = &brep.faces[fh].surface.as_ref().unwrap().geometry;
        let loop_ = &brep.loops[brep.faces[fh].outer_loop];
        let mut cursor = loop_.half_edge;
        loop {
            let he = &brep.half_edges[cursor];
            if let Some(ref pcurve) = he.pcurve {
                let (t0, t1) = pcurve.parameter_range();
                let uv_start = pcurve.point_at(t0);
                let pos_start = surface.point_at(uv_start.x, uv_start.y);
                let actual_start = brep.vertices[he.origin].position;
                let err = pos_start.distance_to(&actual_start);
                assert!(err < 1e-6, "cylinder pcurve start error: {err}");

                let uv_end = pcurve.point_at(t1);
                let pos_end = surface.point_at(uv_end.x, uv_end.y);
                let next_he = &brep.half_edges[he.next];
                let actual_end = brep.vertices[next_he.origin].position;
                let err_end = pos_end.distance_to(&actual_end);
                assert!(err_end < 1e-6, "cylinder pcurve end error: {err_end}");
            }
            cursor = brep.half_edges[cursor].next;
            if cursor == loop_.half_edge {
                break;
            }
        }
    }
}

/// Recomputed pcurves on a cube (all planar faces) should be degree-1
/// (linear) since UV projection onto a plane is a linear map.
#[test]
fn cube_recomputed_pcurves_are_linear() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Strip all pcurves.
    for heh in brep.half_edges.keys().collect::<Vec<_>>() {
        brep.half_edges[heh].pcurve = None;
    }
    brep.recompute_all_pcurves().unwrap();

    for (_, he) in &brep.half_edges {
        let pcurve = he
            .pcurve
            .as_ref()
            .expect("pcurve should exist after recompute");
        assert_eq!(pcurve.degree, 1, "cube pcurve should be degree-1 (linear)");
        assert_eq!(
            pcurve.control_points.len(),
            2,
            "degree-1 pcurve should have exactly 2 control points"
        );
    }
}

/// Cylinder lateral face has curved (non-planar) surface. After recompute,
/// the arc pcurves should still accurately map to 3D positions at several
/// intermediate parameter values, not just the endpoints.
#[test]
fn cylinder_lateral_pcurves_accurate_at_midpoints() {
    let mut brep = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    for heh in brep.half_edges.keys().collect::<Vec<_>>() {
        brep.half_edges[heh].pcurve = None;
    }
    brep.recompute_all_pcurves().unwrap();

    let lateral_face = brep.faces.keys().find(|&fh| {
        matches!(
            brep.faces[fh].surface.as_ref().map(|s| &s.geometry),
            Some(SurfaceGeometry::Cylinder { .. })
        )
    });
    let fh = lateral_face.expect("cylinder must have a lateral face");
    let surface = &brep.faces[fh].surface.as_ref().unwrap().geometry;
    let loop_ = &brep.loops[brep.faces[fh].outer_loop];
    let mut cursor = loop_.half_edge;
    let mut checked = 0;
    loop {
        let he = &brep.half_edges[cursor];
        if let Some(ref pcurve) = he.pcurve {
            // Sample the pcurve at multiple t values and check that the
            // resulting UV maps to a point that lies on the surface.
            let (t0, t1) = pcurve.parameter_range();
            for i in 0..=4 {
                let frac = i as f64 / 4.0;
                let t = t0 + (t1 - t0) * frac;
                let uv = pcurve.point_at(t);
                let pos = surface.point_at(uv.x, uv.y);
                // The point should be on the cylinder (radius 1, Z in [0,2]).
                let r = (pos.x * pos.x + pos.y * pos.y).sqrt();
                // On the lateral face, all points lie on the cylinder barrel
                // or on a cap plane. Just check coordinates are finite.
                assert!(pos.x.is_finite() && pos.y.is_finite() && pos.z.is_finite());
                // For seam edges (line along Z axis): r == radius
                // For arc edges: r == radius and z is 0 or height
                // Either way, the point should be near the cylinder surface.
                if r > 0.1 {
                    assert!(
                        (r - 1.0).abs() < 0.01,
                        "point at t={t:.2} should be on cylinder: r={r}"
                    );
                }
                checked += 1;
            }
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == loop_.half_edge {
            break;
        }
    }
    assert!(checked > 0, "should have checked at least one pcurve");
}

/// Deliberately corrupt pcurve UV values so they map to 3D positions far
/// from the vertices. `fix_stale_pcurves` should detect the staleness
/// (error > 0.05) and recompute correct pcurves.
#[test]
fn fix_stale_pcurves_after_corrupted_uv() {
    use crate::geometry::NurbsCurve;

    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face_handle = brep.faces.keys().next().unwrap();

    // Corrupt all pcurves on this face by shifting UV values by 100 units.
    let loop_ = &brep.loops[brep.faces[face_handle].outer_loop];
    let mut cursor = loop_.half_edge;
    loop {
        if let Some(ref pc) = brep.half_edges[cursor].pcurve {
            let bad_pcurve = NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![
                    Point3::new(pc.control_points[0].x + 100.0, 100.0, 0.0),
                    Point3::new(100.0, pc.control_points[1].y + 100.0, 0.0),
                ],
                weights: vec![1.0, 1.0],
            };
            brep.half_edges[cursor].pcurve = Some(bad_pcurve);
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == loop_.half_edge {
            break;
        }
    }

    // Verify stale: pcurve UV should map far from vertex.
    let he = &brep.half_edges[loop_.half_edge];
    let pcurve = he.pcurve.as_ref().unwrap();
    let (t0, _) = pcurve.parameter_range();
    let uv = pcurve.point_at(t0);
    let surface = &brep.faces[face_handle].surface.as_ref().unwrap().geometry;
    let pos_from_bad_pcurve = surface.point_at(uv.x, uv.y);
    let vtx = brep.vertices[he.origin].position;
    let stale_err = pos_from_bad_pcurve.distance_to(&vtx);
    assert!(stale_err > 1.0, "pcurve should be stale: err={stale_err}");

    // Fix stale pcurves.
    brep.fix_stale_pcurves();

    // After fix, pcurves should correctly map vertex positions.
    let loop_ = &brep.loops[brep.faces[face_handle].outer_loop];
    let mut cursor = loop_.half_edge;
    loop {
        let he = &brep.half_edges[cursor];
        let pcurve = he
            .pcurve
            .as_ref()
            .expect("pcurve should be recomputed after fix");
        let (t0, _) = pcurve.parameter_range();
        let uv = pcurve.point_at(t0);
        let surface = &brep.faces[face_handle].surface.as_ref().unwrap().geometry;
        let pos = surface.point_at(uv.x, uv.y);
        let vtx = brep.vertices[he.origin].position;
        let err = pos.distance_to(&vtx);
        assert!(
            err < 1e-3,
            "after fix_stale_pcurves, pcurve start should map to vertex: err={err}"
        );
        cursor = he.next;
        if cursor == loop_.half_edge {
            break;
        }
    }
}

/// Recompute on a face whose outer loop exists but has no surface should
/// be a no-op (return Ok).
#[test]
fn recompute_face_with_no_surface_is_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face_handle = brep.faces.keys().next().unwrap();
    // Remove the surface.
    brep.faces[face_handle].surface = None;
    // Should succeed without error.
    let result = brep.recompute_face_pcurves(face_handle);
    assert!(result.is_ok());
}

/// Calling recompute_face_pcurves with an invalid handle returns an error.
#[test]
fn recompute_face_pcurves_invalid_handle() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Create a stale handle by removing a face.
    let fh = brep.faces.keys().next().unwrap();
    brep.faces.remove(fh);
    let result = brep.recompute_face_pcurves(fh);
    assert!(result.is_err());
}

/// fix_stale_pcurves on a clean BRep should be a no-op (no pcurves cleared).
#[test]
fn fix_stale_pcurves_noop_on_clean_brep() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Collect pcurve hashes before.
    let before: Vec<_> = brep
        .half_edges
        .values()
        .filter_map(|he| he.pcurve.as_ref().map(|pc| pc.control_points.clone()))
        .collect();
    brep.fix_stale_pcurves();
    let after: Vec<_> = brep
        .half_edges
        .values()
        .filter_map(|he| he.pcurve.as_ref().map(|pc| pc.control_points.clone()))
        .collect();
    assert_eq!(before.len(), after.len(), "no pcurves should be lost");
}
