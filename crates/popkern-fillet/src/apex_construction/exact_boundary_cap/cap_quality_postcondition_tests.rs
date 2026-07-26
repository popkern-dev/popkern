use super::*;
use popkern_core::tolerance::Tolerance;

#[test]
fn a_cap_with_no_internal_seam_reports_no_measurement() {
    // Not "perfectly continuous" — simply not described by this
    // metric. A caller that read `None` as a pass would reinstate the
    // silence-is-success bug this module exists to remove.
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube");
    let single: Vec<FaceHandle> = brep.faces.iter().map(|(f, _)| f).take(1).collect();
    assert!(worst_mid_seam_misalignment(&brep, &single).is_none());
}

#[test]
fn adjacent_planar_faces_measure_their_real_dihedral() {
    // Two cube faces meet at 90 degrees, so |n_a x n_b| = 1. This
    // pins that the metric reads a genuine crease rather than
    // returning zero for everything.
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube");
    let faces: Vec<FaceHandle> = brep.faces.iter().map(|(f, _)| f).collect();
    let worst = worst_mid_seam_misalignment(&brep, &faces).expect("cube faces share seams");
    assert!(
        worst > 0.5,
        "perpendicular faces must read as a large misalignment, got {worst:.3e}"
    );
}
