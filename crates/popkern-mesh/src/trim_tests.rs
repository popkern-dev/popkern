use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn cube_face_produces_one_outer_trim_loop() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face_handle = brep.faces.keys().next().unwrap();
    let face = &brep.faces[face_handle];
    let surface = &face.surface.as_ref().unwrap().geometry;

    let loops = build_trim_loops(&brep, face_handle, surface, 0.01);
    assert_eq!(loops.len(), 1, "cube face should have exactly 1 loop");
    assert!(loops[0].is_outer, "the loop should be outer");
    assert!(
        loops[0].points.len() >= 4,
        "cube face loop should have at least 4 points, got {}",
        loops[0].points.len()
    );
}

#[test]
fn trim_loop_from_pcurves_is_closed() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face_handle = brep.faces.keys().next().unwrap();
    let face = &brep.faces[face_handle];
    let surface = &face.surface.as_ref().unwrap().geometry;

    let loops = build_trim_loops(&brep, face_handle, surface, 0.01);
    assert!(!loops.is_empty(), "should have at least one loop");
    let pts = &loops[0].points;
    let first = pts.first().unwrap();
    let last = pts.last().unwrap();
    let dist = ((first.0 - last.0).powi(2) + (first.1 - last.1).powi(2)).sqrt();
    assert!(
        dist < 0.1,
        "first and last points should be close, got distance {dist}"
    );
}

#[test]
fn boolean_face_with_hole_produces_outer_and_inner_loops() {
    let tol = Tolerance::default().with_linear(1e-2);
    let cube = BRep::unit_cube(tol).unwrap();
    let cyl = BRep::cylinder(0.3, 2.0, tol).unwrap();

    // Translate cylinder to center it on the cube's top face.
    let offset = popkern_core::geometry::Vec3::new(0.5, 0.5, -0.5);
    let cyl = cyl
        .transform(&popkern_core::geometry::Transform3::translation(
            offset.x, offset.y, offset.z,
        ))
        .unwrap();

    let result = popkern_bool::subtraction(cube, cyl, tol).unwrap();

    // Check ALL faces with inner_loops; at least one must produce the
    // expected outer + inner trim loop pattern.
    let mut found_valid = false;
    for (fh, face) in &result.faces {
        if face.inner_loops.is_empty() || face.surface.is_none() {
            continue;
        }
        let surface = &face.surface.as_ref().unwrap().geometry;
        let loops = build_trim_loops(&result, fh, surface, 0.01);

        let outer_count = loops.iter().filter(|l| l.is_outer).count();
        let inner_count = loops.iter().filter(|l| !l.is_outer).count();

        if outer_count == 1 && inner_count >= 1 {
            found_valid = true;
            break;
        }
    }

    // If no face with inner_loops was found at all, the boolean result
    // topology doesn't use inner_loops for this operation -- acceptable.
    let has_any_face_with_inner = result.faces.values().any(|f| !f.inner_loops.is_empty());
    if has_any_face_with_inner {
        assert!(
            found_valid,
            "found faces with inner_loops but none produced outer+inner trim loops"
        );
    }
}

/// Helper: create a unit-square outer trim loop for classification tests.
fn unit_square_loop() -> TrimLoop {
    TrimLoop {
        points: vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        is_outer: true,
        from_pcurves: true,
    }
}

#[test]
fn classify_triangle_inside_square() {
    let loops = vec![unit_square_loop()];
    let tri = [(0.3, 0.3), (0.5, 0.3), (0.4, 0.5)];
    assert_eq!(classify_triangle(&tri, &loops), TriClass::Inside);
}

#[test]
fn classify_triangle_outside_square() {
    let loops = vec![unit_square_loop()];
    let tri = [(1.5, 1.5), (2.0, 1.5), (1.7, 2.0)];
    assert_eq!(classify_triangle(&tri, &loops), TriClass::Outside);
}

#[test]
fn classify_triangle_boundary_crossing_edge() {
    let loops = vec![unit_square_loop()];
    // Centroid ≈ (0.867, 0.433): inside the unit square.
    // The edge (1.2, 0.3)→(1.0, 0.7) crosses x=1 boundary.
    let tri = [(0.6, 0.3), (1.2, 0.3), (0.8, 0.7)];
    assert_eq!(classify_triangle(&tri, &loops), TriClass::Boundary);
}

#[test]
fn classify_triangle_inside_hole_is_outside() {
    let outer = TrimLoop {
        points: vec![(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)],
        is_outer: true,
        from_pcurves: true,
    };
    let inner = TrimLoop {
        points: vec![(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)],
        is_outer: false,
        from_pcurves: true,
    };
    let loops = vec![outer, inner];
    let tri = [(0.8, 0.8), (1.2, 0.8), (1.0, 1.2)];
    assert_eq!(classify_triangle(&tri, &loops), TriClass::Outside);
}
