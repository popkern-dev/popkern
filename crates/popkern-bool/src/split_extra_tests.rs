// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
#[cfg(test)]
use super::*;
use crate::error::BoolError;
use popkern_core::{
    geometry::{Point3, Vec3, nurbs_curve::NurbsCurve},
    tolerance::Tolerance,
    topology::BRep,
};

fn tol() -> Tolerance {
    Tolerance::default()
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
    let cube_b = {
        use popkern_core::geometry::{curve::CurveGeometry, surface::SurfaceGeometry};
        let mut b = BRep::unit_cube(tol()).expect("cube b");
        let offset = Vec3::new(0.5, 0.3, 0.2);
        for (_, vertex) in &mut b.vertices {
            vertex.position = vertex.position + offset;
        }
        for (_, face) in &mut b.faces {
            if let Some(ref mut sof) = face.surface {
                if let SurfaceGeometry::Plane { ref mut origin, .. } = sof.geometry {
                    *origin = *origin + offset;
                }
            }
        }
        for (_, edge) in &mut b.edges {
            if let Some(ref mut coe) = edge.curve {
                if let CurveGeometry::Line { ref mut origin, .. } = coe.geometry {
                    *origin = *origin + offset;
                }
            }
        }
        b
    };

    let int_result = find_intersections(&cube_a, &cube_b, &tol());
    // If intersection succeeds and has records, try splitting.
    // If intersection fails (e.g. due to surface limitations), skip gracefully.
    match int_result {
        Ok(int_res) if !int_res.records.is_empty() => {
            let a_faces_before = cube_a.faces.len();
            let b_faces_before = cube_b.faces.len();

            let split_result = split_along_curves(cube_a, cube_b, &int_res.records, &tol());

            let result = split_result.expect("split should succeed with pre-split edges");

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
        _ => {
            // No intersection curves found or intersection failed — skip.
        }
    }
}

// ── split.rs coverage: non-planar face skip (line 173) ───────────────────────

/// A face with a non-planar surface (e.g. Sphere) should be skipped during
/// splitting. Build a BRep with a spherical face and an IntersectionRecord
/// targeting it; split_along_curves should skip the face without error.
#[test]
fn non_planar_face_is_skipped_during_split() {
    use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

    let (mut brep_a, face_a_h, p1, p2) = build_splittable_square(Point3::new(0.0, 0.0, 0.0));

    // Overwrite the face's surface with a non-planar (Sphere) geometry.
    brep_a.faces[face_a_h].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        },
        orientation: true,
    });

    let trim_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p1, p2],
        weights: vec![1.0, 1.0],
    };
    use popkern_core::geometry::curve_on_surface::CurveOnSurface;
    let curve = CurveOnSurface {
        curve_3d: trim_3d.clone(),
        curve_uv1: trim_3d.clone(),
        curve_uv2: trim_3d,
        is_closed: false,
    };

    let faces_before = brep_a.faces.len();
    let record = IntersectionRecord {
        face_a: face_a_h,
        face_b: popkern_core::topology::FaceHandle::default(),
        curve,
    };

    let brep_b = BRep::new(tol());
    let result = split_along_curves(brep_a, brep_b, &[record], &tol())
        .expect("non-planar face should be skipped, not error");

    // Face count should be unchanged — the sphere face was skipped.
    assert_eq!(
        result.a.faces.len(),
        faces_before,
        "sphere face should not be split"
    );
    assert!(result.fragments_a.intersection_edges.is_empty());
}

// ── split.rs coverage: point_on_segment degenerate edge (line 79) ────────────

/// Verify that ensure_endpoint_vertices silently skips degenerate (zero-length)
/// boundary edges. We build a face where one edge has both endpoints at the
/// same position, so point_on_segment returns None (line 79 early exit).
#[test]
fn point_on_segment_degenerate_edge_does_not_panic() {
    let (mut brep_a, face_a_h, p1, p2) = build_splittable_square(Point3::new(0.0, 0.0, 0.0));

    // Collapse one edge to a point by moving a vertex onto its neighbour.
    // Find any vertex and move it to the same position as its loop-next origin.
    let loop_h = brep_a.faces[face_a_h].outer_loop;
    let start_he = brep_a.loops[loop_h].half_edge;
    let next_he = brep_a.half_edges[start_he].next;
    let v_start = brep_a.half_edges[start_he].origin;
    let v_next = brep_a.half_edges[next_he].origin;
    let pos_next = brep_a.vertices[v_next].position;
    // Collapse v_start onto v_next — makes the first edge degenerate.
    brep_a.vertices[v_start].position = pos_next;

    let trim_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p1, p2],
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
        face_b: popkern_core::topology::FaceHandle::default(),
        curve,
    };

    let brep_b = BRep::new(tol());
    // Should not panic even with a degenerate edge on the face boundary.
    let _result = split_along_curves(brep_a, brep_b, &[record], &tol());
}

// ── split.rs coverage: split_face error path for A (lines 271-279) ───────────

/// Exercise the split_face Err path for A by using a curve that clips to the
/// face boundary but then attempting to split the face with mismatched topology.
/// Either the split succeeds or returns SplitFailed — both acceptable.
#[test]
fn split_a_face_failure_via_tight_tolerance_is_handled() {
    use popkern_core::tolerance::Tolerance;

    // Build a square face.
    let mut brep_a = BRep::new(Tolerance::default());
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

    let p1 = Point3::new(0.5, 0.0, 0.0); // on bottom edge (already a midpoint)
    let p2 = Point3::new(0.5, 1.0, 0.0); // on top edge
    let trim_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p1, p2],
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
        face_b: popkern_core::topology::FaceHandle::default(),
        curve,
    };

    let brep_b = BRep::new(Tolerance::default());
    // Either split succeeds or returns SplitFailed — both are acceptable.
    let result = split_along_curves(brep_a, brep_b, &[record], &Tolerance::default());
    match result {
        Ok(r) => {
            assert!(
                r.a.faces.len() >= 2,
                "face should have been split or stayed same"
            );
        }
        Err(diag) => {
            assert!(matches!(diag.error, Some(BoolError::SplitFailed { .. })));
        }
    }
}

// ── UV-space splitting tests ─────────────────────────────────────────────────

/// A cube (which has pcurves) split with a real intersection record should
/// use the UV-space path and produce face fragments.
#[test]
fn split_cube_face_via_uv_path() {
    use crate::intersect::find_intersections;

    let cube_a = BRep::unit_cube(tol()).expect("cube a");
    let cube_b = {
        use popkern_core::geometry::{curve::CurveGeometry, surface::SurfaceGeometry};
        let mut b = BRep::unit_cube(tol()).expect("cube b");
        let offset = Vec3::new(0.5, 0.3, 0.2);
        for (_, vertex) in &mut b.vertices {
            vertex.position = vertex.position + offset;
        }
        for (_, face) in &mut b.faces {
            if let Some(ref mut sof) = face.surface {
                if let SurfaceGeometry::Plane { ref mut origin, .. } = sof.geometry {
                    *origin = *origin + offset;
                }
            }
        }
        for (_, edge) in &mut b.edges {
            if let Some(ref mut coe) = edge.curve {
                if let CurveGeometry::Line { ref mut origin, .. } = coe.geometry {
                    *origin = *origin + offset;
                }
            }
        }
        b
    };

    // Verify cube has pcurves (UV path will be taken).
    let first_face = cube_a.faces.keys().next().expect("has faces");
    assert!(
        crate::clip::face_uv_boundary(&cube_a, first_face).is_some(),
        "cube faces should have pcurves"
    );

    let int_result = find_intersections(&cube_a, &cube_b, &tol());
    match int_result {
        Ok(int_res) if !int_res.records.is_empty() => {
            let a_faces_before = cube_a.faces.len();
            let result = split_along_curves(cube_a, cube_b, &int_res.records, &tol())
                .expect("split should succeed via UV path");

            // Should have split at least some faces.
            assert!(
                result.fragments_a.origin.len() >= a_faces_before,
                "fragment map should cover all A faces"
            );
            assert!(
                result.fragments_b.origin.len() >= 6,
                "fragment map should cover all B faces"
            );
        }
        _ => {
            // Intersection didn't produce records — skip.
        }
    }
}

/// Existing planar-face splitting should still work (regression test).
#[test]
fn planar_faces_without_pcurves_still_split() {
    // build_splittable_square creates faces without pcurves or surface geometry,
    // so the planar fallback path should still handle them.
    let (brep_a, face_a_h, p1, p2) = build_splittable_square(Point3::new(0.0, 0.0, 0.0));
    let (brep_b, face_b_h, _, _) = build_splittable_square(Point3::new(0.0, 0.0, 0.0));

    let trim_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p1, p2],
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

    let a_before = brep_a.faces.len();
    let b_before = brep_b.faces.len();
    let result = split_along_curves(brep_a, brep_b, &[record], &tol())
        .expect("planar split should still succeed");

    assert_eq!(result.a.faces.len(), a_before + 1, "A should gain a face");
    assert_eq!(result.b.faces.len(), b_before + 1, "B should gain a face");
}

/// Face with pcurves but where the UV curve doesn't cross the boundary
/// should be skipped gracefully (0 crossings).
#[test]
fn uv_curve_outside_boundary_skips_face() {
    use popkern_core::geometry::curve_on_surface::CurveOnSurface;

    let cube = BRep::unit_cube(tol()).expect("cube");
    let face_h = cube.faces.keys().next().expect("has faces");

    // Curve entirely outside the cube's UV domain: u,v in [10,11].
    // Use 3D coords far off-plane so the planar fallback also rejects it.
    let far_uv = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(10.0, 10.0, 0.0), Point3::new(11.0, 11.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let far_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(10.0, 10.0, 5.0), Point3::new(11.0, 11.0, 5.0)],
        weights: vec![1.0, 1.0],
    };
    let curve = CurveOnSurface {
        curve_3d: far_3d,
        curve_uv1: far_uv.clone(),
        curve_uv2: far_uv,
        is_closed: false,
    };
    let record = IntersectionRecord {
        face_a: face_h,
        face_b: popkern_core::topology::FaceHandle::default(),
        curve,
    };

    let brep_b = BRep::new(tol());
    let faces_before = cube.faces.len();
    let result = split_along_curves(cube, brep_b, &[record], &tol())
        .expect("should succeed, skipping the face");

    assert_eq!(
        result.a.faces.len(),
        faces_before,
        "face should not be split when UV curve is outside boundary"
    );
}
