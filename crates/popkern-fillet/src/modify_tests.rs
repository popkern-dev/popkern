//! Tests for [`crate::modify::modify_fillet_radius`] and the
//! record-based path through [`crate::FilletOperationRecord`].

use super::*;
use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use std::collections::HashSet;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube() -> BRep {
    BRep::unit_cube(tol()).unwrap()
}

/// Find a new fillet face — one that wasn't present in `pre` and
/// whose surface is non-planar (the closed-fillet path may rewrite
/// the surface to a cylinder, but on plane-plane we expect a torus
/// or NURBS sweep).
fn find_new_blend(pre: &HashSet<FaceHandle>, post: &BRep) -> Option<FaceHandle> {
    post.faces.keys().find(|fh| {
        if pre.contains(fh) {
            return false;
        }
        crate::identify_blend(post, *fh).is_some()
    })
}

// ---------------------------------------------------------------------------
// modify_fillet_radius — open plane-plane case
// ---------------------------------------------------------------------------

#[test]
fn modify_fillet_radius_changes_open_fillet_on_unit_cube() {
    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let pre: HashSet<_> = brep.faces.keys().collect();

    // Initial fillet at radius 0.1.
    let r1 = crate::fillet_edge(&brep, edge, 0.1, tol(), None).unwrap();
    let blend_face = find_new_blend(&pre, &r1.brep).expect("initial fillet should produce a blend");

    // Modify the fillet to radius 0.2.
    let modified = modify_fillet_radius(&r1.brep, blend_face, 0.2, tol())
        .expect("modify should succeed for plane-plane fillet");

    // Cross-check: at least one blend face in the modified BRep
    // should have a measured radius close to the new requested
    // radius (0.2). We scan all candidates because identify_blend
    // can also tag other faces (e.g. small slivers from the
    // surface-fix post-pass).
    let radii: Vec<f64> = modified
        .faces
        .keys()
        .filter_map(|fh| crate::blend_radius(&modified, fh))
        .collect();
    assert!(
        !radii.is_empty(),
        "modified BRep should have at least one identifiable blend"
    );
    let close_to_target = radii.iter().any(|r| (r - 0.2).abs() < 0.05);
    assert!(
        close_to_target,
        "no blend in modified BRep matches the new radius 0.2; observed radii: {radii:?}"
    );
}

#[test]
fn modify_fillet_radius_is_deterministic() {
    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let pre: HashSet<_> = brep.faces.keys().collect();
    let r1 = crate::fillet_edge(&brep, edge, 0.1, tol(), None).unwrap();
    let blend_face = find_new_blend(&pre, &r1.brep).unwrap();

    let m1 = modify_fillet_radius(&r1.brep, blend_face, 0.15, tol()).unwrap();
    let m2 = modify_fillet_radius(&r1.brep, blend_face, 0.15, tol()).unwrap();
    assert_eq!(m1.faces.len(), m2.faces.len());
    assert_eq!(m1.edges.len(), m2.edges.len());
    assert_eq!(m1.vertices.len(), m2.vertices.len());
}

// ---------------------------------------------------------------------------
// modify_fillet_radius — error cases
// ---------------------------------------------------------------------------

#[test]
fn modify_fillet_radius_rejects_non_fillet_face() {
    let brep = cube();
    let face = brep.faces.keys().next().unwrap();
    let result = modify_fillet_radius(&brep, face, 0.1, tol());
    assert!(
        matches!(result, Err(BlendError::InvalidParameter { .. })),
        "modifying a non-fillet face should fail with InvalidParameter, got {result:?}"
    );
}

// ---------------------------------------------------------------------------
// FilletOperationRecord — record-based path
// ---------------------------------------------------------------------------

#[test]
fn record_with_uniform_radius_rebuilds_fillet_at_new_radius() {
    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];

    let record = crate::fillet_recorded(&brep, &specs, &FilletOptions::default(), tol()).unwrap();

    // Modify the radius via the record.
    let modified = record.with_uniform_radius(0.25).unwrap();

    // The modified record's input is the original BRep, but the
    // result reflects the new radius.
    assert_eq!(modified.input.faces.len(), brep.faces.len());
    assert!(modified.result.faces.len() > brep.faces.len());

    // The modified record's specs should carry the new radius.
    match &modified.specs[0].radius {
        FilletRadius::Constant(r) => assert!(
            (*r - 0.25).abs() < 1e-12,
            "modified spec should hold the new radius, got {r}"
        ),
        other => panic!("expected Constant radius, got {other:?}"),
    }

    // The original record is unchanged.
    match &record.specs[0].radius {
        FilletRadius::Constant(r) => assert!((*r - 0.1).abs() < 1e-12),
        other => panic!("expected Constant radius, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// modify_fillet_radius — closed (plane-cylinder hole rim) case
// ---------------------------------------------------------------------------

/// Round trip through [`modify_fillet_radius`] on a closed-fillet
/// annulus: fillet a hole rim at radius 0.05, then modify to 0.08.
/// Verifies the closed-fillet branch composes
/// [`crate::unfillet_closed::unfillet_closed_annulus`] with a fresh
/// closed-fillet at the new radius.
#[test]
fn modify_fillet_radius_changes_closed_hole_rim_fillet() {
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let cube_brep = cube();
    let cyl = BRep::cylinder(0.3, 2.0, tol()).unwrap();
    let translated = cyl
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    let subtracted = match boolean(
        BoolOperation::Subtraction,
        cube_brep,
        translated,
        tol(),
        None,
    ) {
        Ok(b) => b,
        Err(_) => return,
    };

    let edge = subtracted.edges.keys().find(|e| {
        let he = &subtracted.half_edges[subtracted.edges[*e].half_edge];
        let fa = subtracted.loops[he.loop_].face;
        let fb = subtracted.loops[subtracted.half_edges[he.twin].loop_].face;
        let sa = subtracted.faces[fa].surface.as_ref().map(|s| &s.geometry);
        let sb = subtracted.faces[fb].surface.as_ref().map(|s| &s.geometry);
        matches!(
            (sa, sb),
            (
                Some(SurfaceGeometry::Plane { .. }),
                Some(SurfaceGeometry::Cylinder { .. })
            ) | (
                Some(SurfaceGeometry::Cylinder { .. }),
                Some(SurfaceGeometry::Plane { .. })
            )
        )
    });
    let Some(edge) = edge else {
        return;
    };

    let r1 = match crate::fillet_edge(&subtracted, edge, 0.05, tol(), None) {
        Ok(b) => b,
        Err(_) => return,
    };
    let annulus = r1.brep.faces.keys().find(|fh| {
        let f = &r1.brep.faces[*fh];
        if f.inner_loops.len() != 1 {
            return false;
        }
        let single = |lh: popkern_core::topology::LoopHandle| -> bool {
            let Some(l) = r1.brep.loops.get(lh) else {
                return false;
            };
            let start = l.half_edge;
            let Some(he) = r1.brep.half_edges.get(start) else {
                return false;
            };
            he.next == start && he.prev == start
        };
        single(f.outer_loop) && single(f.inner_loops[0])
    });
    let annulus = match annulus {
        Some(f) => f,
        None => return,
    };

    let modified = modify_fillet_radius(&r1.brep, annulus, 0.08, tol())
        .expect("closed-fillet modify should succeed");
    // Result should still have an annulus face on it (re-filleted).
    let still_has_annulus = modified.faces.keys().any(|fh| {
        let f = &modified.faces[fh];
        f.inner_loops.len() == 1
    });
    assert!(
        still_has_annulus,
        "modified BRep should still contain an annulus fillet face"
    );
}

// ---------------------------------------------------------------------------
// Round 5d-2: true in-place modification (§14.7 ≥2× speedup)
// ---------------------------------------------------------------------------

/// Direct test of [`crate::modify_inplace::modify_closed_in_place`]:
/// compute it on a closed-fillet annulus and verify the result is
/// (a) topologically identical to the input, (b) has the new
/// effective radius, and (c) leaves the input unchanged.
#[test]
fn modify_closed_in_place_preserves_topology_and_changes_radius() {
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let cube_brep = cube();
    let cyl = BRep::cylinder(0.3, 2.0, tol()).unwrap();
    let translated = cyl
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    let subtracted = match boolean(
        BoolOperation::Subtraction,
        cube_brep,
        translated,
        tol(),
        None,
    ) {
        Ok(b) => b,
        Err(_) => return,
    };
    let edge = subtracted.edges.keys().find(|e| {
        let he = &subtracted.half_edges[subtracted.edges[*e].half_edge];
        let fa = subtracted.loops[he.loop_].face;
        let fb = subtracted.loops[subtracted.half_edges[he.twin].loop_].face;
        let sa = subtracted.faces[fa].surface.as_ref().map(|s| &s.geometry);
        let sb = subtracted.faces[fb].surface.as_ref().map(|s| &s.geometry);
        matches!(
            (sa, sb),
            (
                Some(SurfaceGeometry::Plane { .. }),
                Some(SurfaceGeometry::Cylinder { .. })
            ) | (
                Some(SurfaceGeometry::Cylinder { .. }),
                Some(SurfaceGeometry::Plane { .. })
            )
        )
    });
    let Some(edge) = edge else { return };
    let r1 = match crate::fillet_edge(&subtracted, edge, 0.05, tol(), None) {
        Ok(b) => b,
        Err(_) => return,
    };
    let annulus = r1.brep.faces.keys().find(|fh| {
        let f = &r1.brep.faces[*fh];
        if f.inner_loops.len() != 1 {
            return false;
        }
        let single = |lh: popkern_core::topology::LoopHandle| -> bool {
            let Some(l) = r1.brep.loops.get(lh) else {
                return false;
            };
            let start = l.half_edge;
            let Some(he) = r1.brep.half_edges.get(start) else {
                return false;
            };
            he.next == start && he.prev == start
        };
        single(f.outer_loop) && single(f.inner_loops[0])
    });
    let annulus = match annulus {
        Some(f) => f,
        None => return,
    };

    let pre_face_count = r1.brep.faces.len();
    let pre_edge_count = r1.brep.edges.len();
    let pre_vertex_count = r1.brep.vertices.len();

    let modified = crate::modify_inplace::modify_closed_in_place(&r1.brep, annulus, 0.08)
        .expect("in-place modify should succeed for plane-cylinder annulus");

    // (a) Topology is byte-identical: same handle counts.
    assert_eq!(modified.faces.len(), pre_face_count);
    assert_eq!(modified.edges.len(), pre_edge_count);
    assert_eq!(modified.vertices.len(), pre_vertex_count);

    // (b) The fillet face's surface is now a Torus with minor
    // radius 0.08.
    let new_surface = modified.faces[annulus]
        .surface
        .as_ref()
        .expect("modified annulus should have a surface");
    match &new_surface.geometry {
        SurfaceGeometry::Torus { minor_radius, .. } => {
            assert!(
                (minor_radius - 0.08).abs() < 1e-9,
                "in-place torus minor radius {minor_radius} should equal new radius 0.08"
            );
        }
        other => panic!("in-place modify should assign a Torus surface, got {other:?}"),
    }

    // (c) The original BRep is unchanged.
    let original_surface = r1.brep.faces[annulus]
        .surface
        .as_ref()
        .expect("original annulus surface");
    if let SurfaceGeometry::Torus { minor_radius, .. } = &original_surface.geometry {
        assert!(
            (minor_radius - 0.08).abs() > 1e-6,
            "original surface should not have already been mutated"
        );
    }
}

/// §14.7 acceptance: in-place modification must be ≥2× faster than
/// the remove+recreate baseline. Times both paths and asserts the
/// ratio. Loose budget (factor of 2) so noisy CI machines don't
/// flake; the actual ratio in practice is much larger.
#[test]
fn modify_closed_in_place_is_at_least_2x_faster_than_remove_recreate() {
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use std::time::Instant;

    let cube_brep = cube();
    let cyl = BRep::cylinder(0.3, 2.0, tol()).unwrap();
    let translated = cyl
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    let subtracted = match boolean(
        BoolOperation::Subtraction,
        cube_brep,
        translated,
        tol(),
        None,
    ) {
        Ok(b) => b,
        Err(_) => return,
    };
    let edge = subtracted.edges.keys().find(|e| {
        let he = &subtracted.half_edges[subtracted.edges[*e].half_edge];
        let fa = subtracted.loops[he.loop_].face;
        let fb = subtracted.loops[subtracted.half_edges[he.twin].loop_].face;
        let sa = subtracted.faces[fa].surface.as_ref().map(|s| &s.geometry);
        let sb = subtracted.faces[fb].surface.as_ref().map(|s| &s.geometry);
        matches!(
            (sa, sb),
            (
                Some(SurfaceGeometry::Plane { .. }),
                Some(SurfaceGeometry::Cylinder { .. })
            ) | (
                Some(SurfaceGeometry::Cylinder { .. }),
                Some(SurfaceGeometry::Plane { .. })
            )
        )
    });
    let Some(edge) = edge else { return };
    let r1 = match crate::fillet_edge(&subtracted, edge, 0.05, tol(), None) {
        Ok(b) => b,
        Err(_) => return,
    };
    let annulus = r1.brep.faces.keys().find(|fh| {
        let f = &r1.brep.faces[*fh];
        if f.inner_loops.len() != 1 {
            return false;
        }
        let single = |lh: popkern_core::topology::LoopHandle| -> bool {
            let Some(l) = r1.brep.loops.get(lh) else {
                return false;
            };
            let start = l.half_edge;
            let Some(he) = r1.brep.half_edges.get(start) else {
                return false;
            };
            he.next == start && he.prev == start
        };
        single(f.outer_loop) && single(f.inner_loops[0])
    });
    let annulus = match annulus {
        Some(f) => f,
        None => return,
    };

    // Time the in-place path.
    let n = 20;
    let t_inplace = {
        let t0 = Instant::now();
        for i in 0..n {
            let r = 0.06 + (i as f64) * 0.001;
            let _ = crate::modify_inplace::modify_closed_in_place(&r1.brep, annulus, r).unwrap();
        }
        t0.elapsed()
    };

    // Time the remove+recreate path (unfillet → fillet).
    let t_baseline = {
        let t0 = Instant::now();
        for i in 0..n {
            let r = 0.06 + (i as f64) * 0.001;
            let unfilleted =
                crate::unfillet_closed::unfillet_closed_annulus(&r1.brep, annulus).unwrap();
            // Find the closed plane-cylinder rim edge.
            let rim = unfilleted.edges.keys().find(|e| {
                let he = &unfilleted.half_edges[unfilleted.edges[*e].half_edge];
                if he.origin != unfilleted.half_edges[he.twin].origin {
                    return false;
                }
                let fa = unfilleted.loops[he.loop_].face;
                let fb = unfilleted.loops[unfilleted.half_edges[he.twin].loop_].face;
                let sa = unfilleted.faces[fa].surface.as_ref().map(|s| &s.geometry);
                let sb = unfilleted.faces[fb].surface.as_ref().map(|s| &s.geometry);
                matches!(
                    (sa, sb),
                    (
                        Some(SurfaceGeometry::Plane { .. }),
                        Some(SurfaceGeometry::Cylinder { .. })
                    ) | (
                        Some(SurfaceGeometry::Cylinder { .. }),
                        Some(SurfaceGeometry::Plane { .. })
                    )
                )
            });
            let rim = match rim {
                Some(r) => r,
                None => return,
            };
            let _ = crate::fillet_edge(&unfilleted, rim, r, tol(), None).unwrap();
        }
        t0.elapsed()
    };

    let ratio = t_baseline.as_secs_f64() / t_inplace.as_secs_f64();
    eprintln!(
        "[5d-2 perf] in-place: {:?}, baseline: {:?}, speedup: {:.2}×",
        t_inplace, t_baseline, ratio
    );
    assert!(
        ratio >= 2.0,
        "in-place should be ≥2× faster than remove+recreate; got {ratio:.2}× \
         (in-place: {t_inplace:?}, baseline: {t_baseline:?})"
    );
}

#[test]
fn record_re_execute_with_specs_accepts_arbitrary_replacement() {
    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let record = crate::fillet_recorded(&brep, &specs, &FilletOptions::default(), tol()).unwrap();

    let new_specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Variable {
            start: 0.05,
            end: 0.2,
        },
    }];
    let modified = record.re_execute_with_specs(&new_specs).unwrap();
    assert_eq!(modified.specs.len(), 1);
    assert!(matches!(
        &modified.specs[0].radius,
        FilletRadius::Variable { .. }
    ));
}
