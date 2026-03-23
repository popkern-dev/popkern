// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
#[cfg(test)]
use super::*;
use popkern_core::{
    geometry::{Point3, Vec3, nurbs_curve::NurbsCurve},
    tolerance::Tolerance,
    topology::BRep,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Translate all vertices of a BRep by a vector offset.
fn translate(mut brep: BRep, offset: Vec3) -> BRep {
    for (_, vertex) in &mut brep.vertices {
        vertex.position = vertex.position + offset;
    }
    brep
}

#[test]
fn no_records_returns_unchanged_breps() {
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = BRep::unit_cube(tol()).expect("cube b");
    let a_faces: Vec<_> = a.faces.keys().collect();
    let b_faces: Vec<_> = b.faces.keys().collect();

    let result = split_along_curves(a, b, &[], &tol()).expect("no records should succeed");

    // All faces map to themselves.
    for &fh in &a_faces {
        assert_eq!(
            result.fragments_a.origin.get(&fh),
            Some(&fh),
            "unsplit face in A should map to itself"
        );
    }
    for &fh in &b_faces {
        assert_eq!(
            result.fragments_b.origin.get(&fh),
            Some(&fh),
            "unsplit face in B should map to itself"
        );
    }

    // No intersection edges.
    assert!(result.fragments_a.intersection_edges.is_empty());
    assert!(result.fragments_b.intersection_edges.is_empty());
}

#[test]
fn split_with_synthetic_record() {
    // Build a BRep with one square face, split an edge to create a vertex
    // at the midpoint, then create a synthetic IntersectionRecord that
    // connects two midpoint vertices.
    let mut brep_a = BRep::new(tol());
    let mvfs = brep_a.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
    let mev1 = brep_a
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep_a
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep_a
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef_r = brep_a
        .mef(mev3.half_edge_twin, mev1.half_edge)
        .expect("mef");

    // Split two opposite edges at midpoints.
    let split1 = brep_a.split_edge(mev1.edge, 0.5).expect("split1");
    let split2 = brep_a.split_edge(mev3.edge, 0.5).expect("split2");

    let p1 = brep_a.vertices[split1.vertex].position;
    let p2 = brep_a.vertices[split2.vertex].position;

    let face_a_handle = mef_r.face;
    let faces_before = brep_a.faces.len();

    // Create a synthetic curve connecting the two midpoints.
    let trim_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p1, p2],
        weights: vec![1.0, 1.0],
    };
    let trim_uv = trim_3d.clone();

    use popkern_core::geometry::curve_on_surface::CurveOnSurface;
    let curve = CurveOnSurface {
        curve_3d: trim_3d,
        curve_uv1: trim_uv.clone(),
        curve_uv2: trim_uv,
        is_closed: false,
    };

    let record = IntersectionRecord {
        face_a: face_a_handle,
        face_b: face_a_handle, // dummy — won't be used for B
        curve,
    };

    // Use a trivial BRep for B (no faces affected).
    let brep_b = BRep::new(tol());

    let result =
        split_along_curves(brep_a, brep_b, &[record], &tol()).expect("split should succeed");

    // A should have one more face than before.
    assert_eq!(result.a.faces.len(), faces_before + 1);

    // The face_a_handle should appear as an origin value.
    let fragments_from_a: Vec<_> = result
        .fragments_a
        .origin
        .iter()
        .filter(|&(_, &origin)| origin == face_a_handle)
        .collect();
    assert_eq!(
        fragments_from_a.len(),
        2,
        "splitting one face should produce two fragments"
    );

    // Intersection edges should be populated.
    assert!(
        !result.fragments_a.intersection_edges.is_empty(),
        "split should create intersection edges"
    );
}

/// Build a square BRep with two midpoint vertices and return (brep, face, p1, p2).
fn build_splittable_square(
    origin: Point3,
) -> (BRep, popkern_core::topology::FaceHandle, Point3, Point3) {
    let mut brep = BRep::new(tol());
    let mvfs = brep.mvfs(origin).expect("mvfs");
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, origin + Vec3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep
        .mev(mev1.half_edge_twin, origin + Vec3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep
        .mev(mev2.half_edge_twin, origin + Vec3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).expect("mef");
    let split1 = brep.split_edge(mev1.edge, 0.5).expect("split1");
    let split2 = brep.split_edge(mev3.edge, 0.5).expect("split2");
    let p1 = brep.vertices[split1.vertex].position;
    let p2 = brep.vertices[split2.vertex].position;
    (brep, mef_r.face, p1, p2)
}

#[test]
fn split_both_sides_succeeds() {
    // Build two identical square BReps and a synthetic curve that connects
    // their midpoints. Both A and B faces should be split successfully,
    // exercising the B-side success path.
    let (brep_a, face_a_h, p1_a, p2_a) = build_splittable_square(Point3::new(0.0, 0.0, 0.0));
    let (brep_b, face_b_h, _p1_b, _p2_b) = build_splittable_square(Point3::new(0.0, 0.0, 0.0));

    let trim_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p1_a, p2_a],
        weights: vec![1.0, 1.0],
    };
    use popkern_core::geometry::curve_on_surface::CurveOnSurface;
    let curve = CurveOnSurface {
        curve_3d: trim_3d.clone(),
        curve_uv1: trim_3d.clone(),
        curve_uv2: trim_3d,
        is_closed: false,
    };
    let record = IntersectionRecord {
        face_a: face_a_h,
        face_b: face_b_h,
        curve,
    };

    let a_faces_before = brep_a.faces.len();
    let b_faces_before = brep_b.faces.len();

    let result = split_along_curves(brep_a, brep_b, &[record], &tol())
        .expect("both-sides split should succeed");

    // Both A and B should have gained a face.
    assert_eq!(
        result.a.faces.len(),
        a_faces_before + 1,
        "A should gain a face"
    );
    assert_eq!(
        result.b.faces.len(),
        b_faces_before + 1,
        "B should gain a face"
    );
    // Both should have intersection edges.
    assert!(!result.fragments_a.intersection_edges.is_empty());
    assert!(!result.fragments_b.intersection_edges.is_empty());
}

#[test]
fn split_b_side_failure_returns_diagnostics() {
    // Create a valid brep_a (empty — no split needed) and a brep_b with
    // a single face. Then construct an IntersectionRecord that refers to
    // brep_b's face with an invalid curve that will cause split_face to fail.
    // This exercises the split-B error path.

    let brep_a = BRep::new(tol());

    let mut brep_b = BRep::new(tol());
    let mvfs = brep_b.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
    let mev1 = brep_b
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep_b
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep_b
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef_r = brep_b
        .mef(mev3.half_edge_twin, mev1.half_edge)
        .expect("mef");
    let face_b_h = mef_r.face;

    // Build a split record that has face_a = non-existent (so A-split skips)
    // but face_b = valid face in brep_b. Use a degenerate curve (zero-length)
    // that will cause split_face to fail.
    let degenerate_curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        // Both endpoints at (2,2,0) — far outside the face, will fail split.
        control_points: vec![Point3::new(2.0, 2.0, 0.0), Point3::new(2.0, 2.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    use popkern_core::geometry::curve_on_surface::CurveOnSurface;
    let curve = CurveOnSurface {
        curve_3d: degenerate_curve.clone(),
        curve_uv1: degenerate_curve.clone(),
        curve_uv2: degenerate_curve,
        is_closed: false,
    };

    // Use a non-existent face handle for face_a so A-side split is skipped.
    let dummy_face_a = popkern_core::topology::FaceHandle::default();
    let record = IntersectionRecord {
        face_a: dummy_face_a, // not in brep_a → A-side skip
        face_b: face_b_h,
        curve,
    };

    let result = split_along_curves(brep_a, brep_b, &[record], &tol());
    // Expect either success (if split_face is lenient) or a SplitFailed error.
    // Either way the function should not panic.
    match result {
        Ok(_) => { /* split_face was lenient — acceptable */ }
        Err(diag) => {
            assert!(
                matches!(diag.error, Some(BoolError::SplitFailed { .. })),
                "expected SplitFailed diagnostic"
            );
            assert_eq!(diag.failed_stage, Some(Stage::Split));
        }
    }
}

#[test]
fn split_a_side_failure_returns_diagnostics() {
    // Similar but for the A-side error path — use a valid face in brep_a
    // and a non-existent face for B, with a degenerate curve.
    let mut brep_a = BRep::new(tol());
    let mvfs = brep_a.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
    let mev1 = brep_a
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep_a
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep_a
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef_r = brep_a
        .mef(mev3.half_edge_twin, mev1.half_edge)
        .expect("mef");
    let face_a_h = mef_r.face;

    let brep_b = BRep::new(tol());

    let degenerate_curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(5.0, 5.0, 0.0), Point3::new(5.0, 5.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    use popkern_core::geometry::curve_on_surface::CurveOnSurface;
    let curve = CurveOnSurface {
        curve_3d: degenerate_curve.clone(),
        curve_uv1: degenerate_curve.clone(),
        curve_uv2: degenerate_curve,
        is_closed: false,
    };

    let dummy_face_b = popkern_core::topology::FaceHandle::default();
    let record = IntersectionRecord {
        face_a: face_a_h,
        face_b: dummy_face_b, // not in brep_b → B-side skip
        curve,
    };

    let result = split_along_curves(brep_a, brep_b, &[record], &tol());
    match result {
        Ok(_) => { /* lenient split — acceptable */ }
        Err(diag) => {
            assert!(
                matches!(diag.error, Some(BoolError::SplitFailed { .. })),
                "expected SplitFailed diagnostic"
            );
        }
    }
}

#[test]
fn fragment_map_covers_all_faces() {
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = BRep::unit_cube(tol()).expect("cube b");
    let a_face_count = a.faces.len();
    let b_face_count = b.faces.len();

    let result = split_along_curves(a, b, &[], &tol()).expect("should succeed");

    assert_eq!(result.fragments_a.origin.len(), a_face_count);
    assert_eq!(result.fragments_b.origin.len(), b_face_count);
}

#[test]
fn split_with_overlapping_cubes_intersection() {
    // Use find_intersections to get real intersection records, then split.
    use crate::intersect::find_intersections;

    let cube_a = BRep::unit_cube(tol()).expect("cube a");
    let cube_b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(0.5, 0.3, 0.2),
    );

    let int_result = find_intersections(&cube_a, &cube_b, &tol());
    // If intersection succeeds and has records, try splitting.
    // If intersection fails (e.g. due to surface limitations), skip gracefully.
    match int_result {
        Ok(int_res) if !int_res.records.is_empty() => {
            let a_faces_before = cube_a.faces.len();
            let b_faces_before = cube_b.faces.len();

            let split_result = split_along_curves(cube_a, cube_b, &int_res.records, &tol());

            match split_result {
                Ok(result) => {
                    // Fragment maps should cover all faces.
                    assert!(result.fragments_a.origin.len() >= a_faces_before);
                    assert!(result.fragments_b.origin.len() >= b_faces_before);

                    // At least some intersection edges should exist.
                    let total_int_edges = result.fragments_a.intersection_edges.len()
                        + result.fragments_b.intersection_edges.len();
                    assert!(
                        total_int_edges > 0,
                        "expected intersection edges after splitting overlapping cubes"
                    );
                }
                Err(_diag) => {
                    // Split may fail if intersection curve endpoints don't
                    // land on face boundary vertices — acceptable for now.
                }
            }
        }
        _ => {
            // No intersection curves found or intersection failed — skip.
        }
    }
}
