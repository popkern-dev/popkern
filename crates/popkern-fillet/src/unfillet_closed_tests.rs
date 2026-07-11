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
    let subtracted = boolean(
        BoolOperation::Subtraction,
        cube_brep,
        translated,
        tol(),
        None,
    )
    .expect("cube-minus-cylinder must subtract (round_5c4 fixture)");

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
    let edge = edge.expect("subtraction must leave a plane-cylinder hole rim");

    let pre_face_count = subtracted.faces.len();
    let pre_edge_count = subtracted.edges.len();
    let pre_vertex_count = subtracted.vertices.len();
    let pre_faces: std::collections::HashSet<_> = subtracted.faces.keys().collect();

    // The Err arm used to `return`, passing VACUOUSLY. The fillet
    // must succeed.
    let result = crate::fillet_edge(&subtracted, edge, 0.05, tol(), None)
        .expect("closed-fillet on the hole rim must succeed");

    // The construction ships the SEAMED band (two contact rings
    // joined by a cross-section seam — apply_closed_fillet step 14),
    // not the pre-seam annulus. Locate it by its signature and
    // round-trip through unfillet_face, which un-seams (kemr) and
    // reverses the annulus surgery.
    // Search NEW faces only: the ring+seam+ring+seam signature also
    // matches the ORIGINAL cylinder wall (any periodic wall between
    // two closed rims has it) — the band is the one the fillet added.
    let band_face = result
        .brep
        .faces
        .keys()
        .find(|fh| {
            !pre_faces.contains(fh)
                && crate::unfillet_closed::seamed_band_bridge(&result.brep, *fh).is_some()
        })
        .expect("could not locate the seamed closed-fillet band face after fillet");

    let restored = match unfillet_face(&result.brep, band_face) {
        Ok(b) => b,
        Err(e) => panic!("unfillet of the seamed closed-fillet band failed: {e:?}"),
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
