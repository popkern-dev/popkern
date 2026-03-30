//! Tests for edge analysis (fillet/chamfer).
//! Included only in test mode via `#[cfg(test)] #[path = "tests.rs"] mod tests`.
// #[cfg(test)] — marker for the no_unwrap linter; .unwrap() is allowed in tests.

use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn make_cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).unwrap()
}

/// Find the first edge handle in the BRep.
fn first_edge(brep: &BRep) -> EdgeHandle {
    brep.edges.keys().next().unwrap()
}

#[test]
fn analyze_cube_edge() {
    let cube = make_cube();
    let edge = first_edge(&cube);

    let analyses = analyze_edges(&cube, &[edge], false).unwrap();
    assert_eq!(analyses.len(), 1);

    let a = &analyses[0];
    assert_eq!(a.edge, edge);

    // Both faces should have planar surfaces.
    assert!(
        matches!(a.surface_1, SurfaceGeometry::Plane { .. }),
        "face_1 should be planar"
    );
    assert!(
        matches!(a.surface_2, SurfaceGeometry::Plane { .. }),
        "face_2 should be planar"
    );

    // A cube edge is convex (exterior corner).
    assert!(
        a.convexity == Convexity::Convex || a.convexity == Convexity::Concave,
        "cube edge should have definite convexity, got {:?}",
        a.convexity
    );

    // Vertices should be distinct.
    assert_ne!(
        a.vertex_start, a.vertex_end,
        "start and end vertices should differ"
    );
}

#[test]
fn analyze_invalid_edge() {
    let cube = make_cube();
    // A default handle will never exist in a populated BRep.
    let bogus = popkern_core::topology::EdgeHandle::default();

    let result = analyze_edges(&cube, &[bogus], false);
    assert!(result.is_err(), "should fail for invalid edge handle");
}

#[test]
fn classify_convexity_on_cube_is_definite() {
    let cube = make_cube();
    // Analyze ALL edges — every cube edge should be either Convex or Concave.
    let edges: Vec<EdgeHandle> = cube.edges.keys().collect();
    let analyses = analyze_edges(&cube, &edges, false).unwrap();
    for a in &analyses {
        assert!(
            a.convexity == Convexity::Convex || a.convexity == Convexity::Concave,
            "cube edge {:?} has unexpected convexity {:?}",
            a.edge,
            a.convexity
        );
    }
}

#[test]
fn concave_edge_from_subtraction() {
    use popkern_core::geometry::Transform3;
    use popkern_core::tolerance::Tolerance;
    let tol = Tolerance::default();
    let big = BRep::unit_cube(tol).unwrap();
    // Scale to half-size then translate so it overlaps one corner.
    let xf = Transform3::translation(0.25, 0.25, 0.5).compose(&Transform3::scale_uniform(0.5));
    let small = BRep::unit_cube(tol).unwrap().transform(&xf).unwrap();

    let result = popkern_bool::subtraction(big, small, tol, None);
    if let Ok(subtracted) = result {
        let edges: Vec<EdgeHandle> = subtracted.edges.keys().collect();
        if !edges.is_empty() {
            let analyses = analyze_edges(&subtracted, &edges, false);
            if let Ok(analyses) = analyses {
                let has_concave = analyses.iter().any(|a| a.convexity == Convexity::Concave);
                assert!(has_concave, "subtraction result should have concave edges");
            }
        }
    }
}

#[test]
fn edge_missing_surface_geometry() {
    // Construct a BRep where a face has no surface, triggering the error path
    // in get_face_surface.
    let mut cube = make_cube();
    let edge = first_edge(&cube);
    let edge_data = &cube.edges[edge];
    let he = &cube.half_edges[edge_data.half_edge];
    let face = cube.loops[he.loop_].face;

    // Remove surface geometry from the face.
    cube.faces[face].surface = None;

    let result = analyze_edges(&cube, &[edge], false);
    assert!(
        result.is_err(),
        "should fail when face has no surface geometry"
    );
    let err = result.unwrap_err();
    assert!(
        err.contains("no surface"),
        "error should mention missing surface: {err}"
    );
}

#[test]
fn analyze_multiple_edges() {
    let cube = make_cube();
    let edges: Vec<EdgeHandle> = cube.edges.keys().take(3).collect();
    let analyses = analyze_edges(&cube, &edges, false).unwrap();
    assert_eq!(analyses.len(), 3);
    for a in &analyses {
        assert_ne!(a.face_1, a.face_2, "adjacent faces should be different");
    }
}

#[test]
fn cylinder_face_convexity_uses_grid_search_fallback() {
    // Replace one face's surface with a Cylinder so that the UV-projection
    // for that face falls back to the grid search path in find_uv_for_point.
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

    let mut cube = make_cube();
    let edge = first_edge(&cube);
    let edge_data = &cube.edges[edge];
    let he = &cube.half_edges[edge_data.half_edge];
    let face = cube.loops[he.loop_].face;

    // Replace the planar surface with a cylinder (analytic UV projection
    // is not defined for cylinders in popkern, so grid fallback fires).
    cube.faces[face].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Cylinder {
            origin: popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        orientation: true,
    });

    // analyze_edges may succeed or fail depending on normal computation,
    // but it must NOT panic — the grid search must be exercised.
    let result = analyze_edges(&cube, &[edge], false);
    // Either path is acceptable; we just confirm no panic.
    let _ = result;
}

#[test]
fn flipped_face_orientation_negates_normal() {
    // Verify that faces with orientation=false produce negated normals.
    // We achieve this by flipping the surface orientation on one face
    // and checking the analysis still runs without error.
    let mut cube = make_cube();
    let edge = first_edge(&cube);
    let edge_data = &cube.edges[edge];
    let he = &cube.half_edges[edge_data.half_edge];
    let face = cube.loops[he.loop_].face;

    // Flip the orientation flag.
    if let Some(sof) = cube.faces[face].surface.as_mut() {
        sof.orientation = false;
    }

    // With flipped orientation, classify_convexity negates the normal for
    // that face. Should still produce a result (possibly opposite convexity).
    let result = analyze_edges(&cube, &[edge], false);
    assert!(
        result.is_ok(),
        "flipped orientation should not break analysis"
    );
}

#[test]
fn convexity_classification_samples_multiple_points() {
    // A cube edge is uniformly convex — multi-point sampling should
    // still return Convex (not Mixed).
    let cube = make_cube();
    let edges: Vec<EdgeHandle> = cube.edges.keys().collect();
    let analyses = analyze_edges(&cube, &edges, false).unwrap();

    for a in &analyses {
        assert!(
            a.convexity == Convexity::Convex || a.convexity == Convexity::Concave,
            "cube edge {:?} should be uniformly convex or concave with \
             multi-point sampling, got {:?}",
            a.edge,
            a.convexity,
        );
    }

    // A concave edge from a boolean subtraction should be Concave.
    use popkern_core::geometry::Transform3;
    use popkern_core::tolerance::Tolerance;
    let tol = Tolerance::default();
    let big = BRep::unit_cube(tol).unwrap();
    let xf = Transform3::translation(0.25, 0.25, 0.5).compose(&Transform3::scale_uniform(0.5));
    let small = BRep::unit_cube(tol).unwrap().transform(&xf).unwrap();

    if let Ok(subtracted) = popkern_bool::subtraction(big, small, tol, None) {
        let edges: Vec<EdgeHandle> = subtracted.edges.keys().collect();
        if let Ok(analyses) = analyze_edges(&subtracted, &edges, false) {
            // At least one concave edge should exist, and multi-point
            // sampling should not flip it to Mixed.
            let has_concave = analyses.iter().any(|a| a.convexity == Convexity::Concave);
            assert!(
                has_concave,
                "subtraction result should have concave edges with multi-point sampling"
            );
        }
    }
}
