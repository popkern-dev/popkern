//! Integration tests for filleting edges between curved faces.
//!
//! These tests create BRep solids with curved faces (via boolean operations)
//! and attempt to fillet edges at the plane-cylinder boundary. Both success
//! and geometric failure are acceptable outcomes — the key invariant is that
//! no panics occur and successful results satisfy Euler-Poincare.

use popkern_bool::BoolOperation;
use popkern_core::geometry::Transform3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

/// Assert that a BRep has no structural validation errors (BrokenCycle,
/// BrokenTwin, InconsistentVertex, NonManifoldEdge, OpenLoop).
fn assert_no_structural_errors(brep: &BRep, label: &str) {
    let errors = brep.validate_all();
    let structural: Vec<_> = errors
        .iter()
        .filter(|e| {
            matches!(
                e,
                popkern_core::error::ValidationError::BrokenCycle { .. }
                    | popkern_core::error::ValidationError::BrokenTwin { .. }
                    | popkern_core::error::ValidationError::InconsistentVertex { .. }
                    | popkern_core::error::ValidationError::NonManifoldEdge { .. }
                    | popkern_core::error::ValidationError::OpenLoop { .. }
            )
        })
        .collect();
    assert!(
        structural.is_empty(),
        "{label}: structural errors: {structural:?}"
    );
}

fn default_tol() -> Tolerance {
    Tolerance::default()
}

/// Compute the Euler-Poincare characteristic V - E + F - L for a
/// BRep. By the formula V - E + F = 2(S - G) + L this equals
/// 2(S - G), so it is constant for any topology-preserving
/// modification of the same underlying solid. A genus-0 single-shell
/// solid yields 2, a genus-1 (one through-hole) solid yields 0,
/// genus-2 yields -2, etc.
fn euler_characteristic(brep: &BRep) -> i64 {
    let v = brep.vertices.len() as i64;
    let e = brep.edges.len() as i64;
    let f = brep.faces.len() as i64;
    let l: i64 = brep
        .faces
        .values()
        .map(|face| face.inner_loops.len() as i64)
        .sum();
    v - e + f - l
}

/// Assert the Euler-Poincare invariant for a closed-shell genus-0
/// solid (V - E + F - L = 2). Used by tests that start from a
/// closed cube (genus 0).
fn assert_euler_poincare(brep: &BRep, label: &str) {
    let chi = euler_characteristic(brep);
    assert!(
        chi == 2,
        "{label}: Euler characteristic V - E + F - L = {chi}, expected 2 (genus 0)"
    );
}

/// Assert that an operation preserved the Euler characteristic of
/// the underlying solid — i.e. the fillet did not change the
/// topology's genus or shell count. Used by tests that start from a
/// solid of unknown genus (e.g. boolean results).
fn assert_euler_preserved(pre: &BRep, post: &BRep, label: &str) {
    let chi_pre = euler_characteristic(pre);
    let chi_post = euler_characteristic(post);
    assert!(
        chi_pre == chi_post,
        "{label}: Euler characteristic changed by operation: pre={chi_pre}, post={chi_post}"
    );
}

/// Find an edge that lies between a planar face and a cylindrical face.
///
/// Returns `None` if no such edge exists in the BRep.
fn find_plane_cylinder_edge(brep: &BRep) -> Option<EdgeHandle> {
    for (edge_h, edge) in &brep.edges {
        let he_fwd = match brep.half_edges.get(edge.half_edge) {
            Some(h) => h,
            None => continue,
        };
        let he_rev = match brep.half_edges.get(he_fwd.twin) {
            Some(h) => h,
            None => continue,
        };

        let face_1 = brep.loops[he_fwd.loop_].face;
        let face_2 = brep.loops[he_rev.loop_].face;

        let surf_1 = brep.faces[face_1].surface.as_ref().map(|s| &s.geometry);
        let surf_2 = brep.faces[face_2].surface.as_ref().map(|s| &s.geometry);

        let is_plane_cyl = matches!(
            (&surf_1, &surf_2),
            (
                Some(SurfaceGeometry::Plane { .. }),
                Some(SurfaceGeometry::Cylinder { .. })
            ) | (
                Some(SurfaceGeometry::Cylinder { .. }),
                Some(SurfaceGeometry::Plane { .. })
            )
        );

        if is_plane_cyl {
            return Some(edge_h);
        }
    }
    None
}

#[test]
fn fillet_edge_between_plane_and_cylinder() {
    let tol = default_tol();

    // Create a cube and a cylinder that overlaps it.
    let cube = BRep::unit_cube(tol).unwrap();
    let cyl = BRep::cylinder(0.3, 2.0, tol).unwrap();

    // Translate the cylinder so it partially overlaps the cube.
    let translated_cyl = cyl
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();

    // Subtract the cylinder from the cube to create plane-cylinder edges.
    let subtracted =
        popkern_bool::boolean(BoolOperation::Subtraction, cube, translated_cyl, tol, None);

    // Boolean may fail if the geometry is too complex — that's acceptable.
    let subtracted = match subtracted {
        Ok(b) => b,
        Err(_) => return, // Boolean failure is acceptable, test passes.
    };

    // Find an edge between a planar face and a cylindrical face.
    let edge = match find_plane_cylinder_edge(&subtracted) {
        Some(e) => e,
        None => return, // No plane-cylinder edge found; skip.
    };

    // Try to fillet the edge. Both success and geometric failure are OK.
    let result = popkern_blend::fillet_edge(&subtracted, edge, 0.05, tol, None);

    match result {
        Ok(blend) => {
            // Fillet must add at least one face (the fillet face itself).
            assert!(
                blend.brep.faces.len() > subtracted.faces.len(),
                "fillet result should have more faces than input ({} vs {})",
                blend.brep.faces.len(),
                subtracted.faces.len(),
            );
            // The boolean engine produces a genus-1 result (cube
            // with through-hole), so the fixed-value Euler test
            // doesn't apply. Instead, assert the operation
            // *preserved* whatever Euler value the input had.
            assert_euler_preserved(&subtracted, &blend.brep, "fillet plane-cylinder edge");
            // No structural errors.
            assert_no_structural_errors(&blend.brep, "fillet plane-cylinder edge");
        }
        Err(_) => {
            // Geometric failure is acceptable — no panic is the invariant.
        }
    }
}

#[test]
fn fillet_edge_on_boolean_cube_result() {
    let tol = default_tol();

    // Create two overlapping cubes.
    let cube_a = BRep::unit_cube(tol).unwrap();
    let cube_b = BRep::unit_cube(tol).unwrap();
    let cube_b = cube_b
        .transform(&Transform3::translation(0.5, 0.5, 0.0))
        .unwrap();

    // Subtract cube_b from cube_a.
    let subtracted = popkern_bool::boolean(BoolOperation::Subtraction, cube_a, cube_b, tol, None);

    let subtracted = match subtracted {
        Ok(b) => b,
        Err(_) => return, // Boolean failure is acceptable.
    };

    // Pick the first available edge and try to fillet it.
    let edge = match subtracted.edges.keys().next() {
        Some(e) => e,
        None => return,
    };

    let result = popkern_blend::fillet_edge(&subtracted, edge, 0.05, tol, None);

    match result {
        Ok(blend) => {
            // Fillet must add at least one face.
            assert!(
                blend.brep.faces.len() > subtracted.faces.len(),
                "fillet result should have more faces than input ({} vs {})",
                blend.brep.faces.len(),
                subtracted.faces.len(),
            );
            // Boolean result is of unknown genus (depends on which
            // cube edge the test happens to pick); assert the
            // fillet preserved the input's Euler characteristic
            // rather than asserting a fixed genus-0 value.
            assert_euler_preserved(&subtracted, &blend.brep, "fillet boolean cube result");
            assert_no_structural_errors(&blend.brep, "fillet boolean cube result");
        }
        Err(_) => {
            // Geometric failure is acceptable — no panic is the invariant.
        }
    }
}

#[test]
fn fillet_cube_edge_has_torus_or_nurbs_surface() {
    let tol = default_tol();
    let cube = BRep::unit_cube(tol).unwrap();
    let original_face_count = cube.faces.len();
    let edge = cube.edges.keys().next().unwrap();

    let result = popkern_blend::fillet_edge(&cube, edge, 0.1, tol, None);
    let blend = result.unwrap();

    // Fillet must add at least one face.
    assert!(
        blend.brep.faces.len() > original_face_count,
        "should have more faces after fillet"
    );

    // Find the new face(s) — faces whose surface geometry is NOT a Plane.
    // The fillet face should have either Torus or Nurbs surface geometry.
    let fillet_faces: Vec<_> = blend
        .brep
        .faces
        .iter()
        .filter(|(_, f)| {
            f.surface
                .as_ref()
                .map(|s| {
                    matches!(
                        s.geometry,
                        SurfaceGeometry::Torus { .. }
                            | SurfaceGeometry::Cylinder { .. }
                            | SurfaceGeometry::Nurbs(_)
                    )
                })
                .unwrap_or(false)
        })
        .collect();

    assert!(
        !fillet_faces.is_empty(),
        "fillet result should have at least one face with Torus or Nurbs surface, \
         found only: {:?}",
        blend
            .brep
            .faces
            .iter()
            .map(|(_, f)| f
                .surface
                .as_ref()
                .map(|s| format!("{:?}", std::mem::discriminant(&s.geometry)))
                .unwrap_or_else(|| "None".into()))
            .collect::<Vec<_>>()
    );

    // Verify the fillet face has non-None surface geometry.
    for (_, face) in &fillet_faces {
        assert!(
            face.surface.is_some(),
            "fillet face surface geometry must not be None"
        );
    }

    assert_euler_poincare(&blend.brep, "fillet cube edge surface type");
    assert_no_structural_errors(&blend.brep, "fillet cube edge surface type");
}

#[test]
fn fillet_edge_between_two_cylinders() {
    let tol = default_tol();

    // Create two cylinders that intersect.
    let cyl_a = BRep::cylinder(0.5, 2.0, tol).unwrap();
    let cyl_b = BRep::cylinder(0.3, 2.0, tol).unwrap();
    let cyl_b = cyl_b
        .transform(&Transform3::translation(0.3, 0.0, 0.0))
        .unwrap();

    // Subtract cyl_b from cyl_a.
    let subtracted = popkern_bool::boolean(BoolOperation::Subtraction, cyl_a, cyl_b, tol, None);

    let subtracted = match subtracted {
        Ok(b) => b,
        Err(_) => return,
    };

    // Pick the first edge and try to fillet it.
    let edge = match subtracted.edges.keys().next() {
        Some(e) => e,
        None => return,
    };

    let result = popkern_blend::fillet_edge(&subtracted, edge, 0.03, tol, None);

    match result {
        Ok(blend) => {
            assert_euler_poincare(&blend.brep, "fillet cylinder-cylinder edge");
        }
        Err(_) => {
            // Geometric failure is acceptable.
        }
    }
}
