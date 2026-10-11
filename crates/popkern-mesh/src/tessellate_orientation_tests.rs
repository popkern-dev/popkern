//! Orientation votes see each face's boundary in one fixed order.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::TessellationOptions;

/// Two extractions of the same face boundary list the same segments in
/// the same order. Each counting map gets its own random hash seed, so
/// without the sort two calls disagree (PK-796).
#[test]
fn boundary_segments_come_out_in_one_order() {
    let brep = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    let mesh = crate::tessellate(&brep, &TessellationOptions::default()).unwrap();
    let face = mesh.face_map[0];
    let faces = crate::face_pass::tessellate_face_pass(
        &brep,
        &[face],
        Default::default(),
        &crate::boundary_first::sample_all_edges(&brep, &TessellationOptions::default()),
        &TessellationOptions::default(),
        false,
        None,
    );
    let fm = match faces.into_iter().next().unwrap() {
        crate::face_pass::FacePassOutcome::BoundaryFirst(fm) => fm,
        crate::face_pass::FacePassOutcome::UniformGrid { mesh, .. } => mesh.unwrap(),
    };
    let first = super::boundary_segments(&fm);
    assert!(first.len() > 16, "a boundary with {} segments", first.len());
    for _ in 0..8 {
        assert_eq!(super::boundary_segments(&fm), first);
    }
}
