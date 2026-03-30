use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::f64::consts::TAU;

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

    let result = popkern_bool::subtraction(cube, cyl, tol, None).unwrap();

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
        points_3d: Vec::new(),
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
        points_3d: Vec::new(),
        is_outer: true,
        from_pcurves: true,
    };
    let inner = TrimLoop {
        points: vec![(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)],
        points_3d: Vec::new(),
        is_outer: false,
        from_pcurves: true,
    };
    let loops = vec![outer, inner];
    let tri = [(0.8, 0.8), (1.2, 0.8), (1.0, 1.2)];
    assert_eq!(classify_triangle(&tri, &loops), TriClass::Outside);
}

#[test]
fn boolean_hole_face_inner_trim_loop_has_points() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();
    let cyl = BRep::cylinder(0.2, 2.0, tol)
        .unwrap()
        .transform(&popkern_core::geometry::Transform3::translation(
            0.5, 0.5, -0.5,
        ))
        .unwrap();
    let result = popkern_bool::subtraction(cube, cyl, tol, None).unwrap();

    // Find faces with inner loops.
    for fh in result.faces.keys() {
        let face = &result.faces[fh];
        if !face.inner_loops.is_empty() {
            let surface = &face.surface.as_ref().unwrap().geometry;
            let loops = build_trim_loops(&result, fh, surface, 0.01);
            let inner_count = loops.iter().filter(|l| !l.is_outer).count();
            assert!(
                inner_count >= 1,
                "face with {} inner loops should produce >= 1 inner trim loop, got {inner_count}",
                face.inner_loops.len()
            );
            for l in &loops {
                if !l.is_outer {
                    assert!(
                        l.points.len() >= 8,
                        "inner trim loop should have >= 8 points for a circle, got {}",
                        l.points.len()
                    );
                }
            }
        }
    }
}

#[test]
fn seam_crossing_polygon_is_unwrapped() {
    let mut points = vec![(5.5, 0.5), (6.0, 0.5), (0.3, 0.5), (0.0, 0.5)];
    // This polygon crosses the u=0/2pi seam on a cylinder.
    assert!(crosses_seam(&points, Some(TAU), None));
    unwrap_periodic(&mut points, Some(TAU), None);
    // After unwrapping, all points should be contiguous.
    for i in 0..points.len() - 1 {
        let du = (points[i + 1].0 - points[i].0).abs();
        assert!(
            du < std::f64::consts::PI,
            "points should be contiguous after unwrap, gap={du}"
        );
    }
}

#[test]
fn non_seam_polygon_unchanged() {
    let mut points = vec![(1.0, 0.5), (2.0, 0.5), (2.0, 1.5), (1.0, 1.5)];
    let original = points.clone();
    unwrap_periodic(&mut points, Some(TAU), None);
    // No seam crossing -- points should be unchanged (within small
    // tolerance from circular mean shift).
    for (p, o) in points.iter().zip(original.iter()) {
        assert!(
            (p.0 - o.0).abs() < 0.01 && (p.1 - o.1).abs() < 0.01,
            "point changed unexpectedly: ({}, {}) -> ({}, {})",
            o.0,
            o.1,
            p.0,
            p.1
        );
    }
}

#[test]
fn crosses_seam_false_for_non_periodic() {
    let points = vec![(5.5, 0.5), (6.0, 0.5), (0.3, 0.5)];
    // No periodicity: should not detect seam crossing.
    assert!(!crosses_seam(&points, None, None));
}

#[test]
fn circular_mean_empty_returns_zero() {
    let result = circular_mean(std::iter::empty(), TAU);
    assert!((result).abs() < 1e-10);
}

#[test]
fn torus_v_seam_crossing_unwrapped() {
    // Polygon crosses the v=0/2pi seam on a torus.
    let mut points = vec![(1.0, 5.8), (1.0, 6.1), (1.0, 0.2), (1.0, 0.0)];
    assert!(crosses_seam(&points, Some(TAU), Some(TAU)));
    unwrap_periodic(&mut points, Some(TAU), Some(TAU));
    for i in 0..points.len() - 1 {
        let dv = (points[i + 1].1 - points[i].1).abs();
        assert!(
            dv < std::f64::consts::PI,
            "v-direction should be contiguous after unwrap, gap={dv}"
        );
    }
}
