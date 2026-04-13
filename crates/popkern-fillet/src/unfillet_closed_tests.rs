//! Tests for [`crate::unfillet_closed`] — round-trip fillet+unfillet
//! on closed (cylindrical hole rim) edges.
//!
//! Kept in its own file so neither the production module nor the
//! main `unfillet_tests.rs` exceeds the 500-line structural cap.

use crate::unfillet::unfillet_face;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Round-trip test for Round 5c-4: fillet a cylindrical hole rim,
/// then unfillet it. The result should restore the same face,
/// vertex, and edge count as the pre-fillet input.
#[test]
fn round_5c4_closed_fillet_round_trip() {
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let cube_brep = BRep::unit_cube(tol()).unwrap();
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

    let pre_face_count = subtracted.faces.len();
    let pre_edge_count = subtracted.edges.len();
    let pre_vertex_count = subtracted.vertices.len();

    let result = match crate::fillet_edge(&subtracted, edge, 0.05, tol(), None) {
        Ok(b) => b,
        Err(_) => return,
    };

    // Find the new closed-fillet annulus face — both loops are
    // single-element self-closed, the topological signature the
    // closed-fillet construction produces. Surface-kind-agnostic
    // because the surface-fix post-pass may overwrite the NURBS
    // sweep with a nearby cylinder.
    let annulus_face = result.brep.faces.keys().find(|fh| {
        let f = &result.brep.faces[*fh];
        if f.inner_loops.len() != 1 {
            return false;
        }
        let single = |lh: popkern_core::topology::LoopHandle| -> bool {
            let Some(l) = result.brep.loops.get(lh) else {
                return false;
            };
            let start = l.half_edge;
            let Some(he) = result.brep.half_edges.get(start) else {
                return false;
            };
            he.next == start && he.prev == start
        };
        single(f.outer_loop) && single(f.inner_loops[0])
    });
    let annulus_face = match annulus_face {
        Some(f) => f,
        None => panic!("could not locate closed-fillet annulus face after fillet"),
    };

    let restored = match unfillet_face(&result.brep, annulus_face) {
        Ok(b) => b,
        Err(e) => panic!("unfillet of closed annulus failed: {e:?}"),
    };
    assert_eq!(
        restored.faces.len(),
        pre_face_count,
        "round trip should restore original face count: pre={pre_face_count}, post-fillet={}, post-unfillet={}",
        result.brep.faces.len(),
        restored.faces.len()
    );
    assert_eq!(
        restored.vertices.len(),
        pre_vertex_count,
        "vertex count should match pre-fillet"
    );
    assert_eq!(
        restored.edges.len(),
        pre_edge_count,
        "edge count should match pre-fillet"
    );
}
