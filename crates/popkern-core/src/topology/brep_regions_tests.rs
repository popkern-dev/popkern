use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle};

#[test]
fn a_cube_has_one_solid_region_bounded_by_its_one_face_set() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let regions = cube.solid_regions();
    assert_eq!(regions.len(), 1);
    let shells = cube.back_shells();
    assert_eq!(shells.len(), 1);
    assert_eq!(cube.back_shell_count(), 1);
    assert_eq!(cube.region_shells(regions[0]), shells);
    assert_eq!(cube.region_outer_shell(regions[0]), Some(shells[0]));
    assert_eq!(cube.shell_region(shells[0]), Some(regions[0]));
    assert_eq!(cube.shell_faces(shells[0]).len(), 6);
    for f in cube.faces.keys() {
        assert_eq!(cube.face_back_shell(f), Some(shells[0]));
    }
}

/// A cavity's shell joins the region after its outer shell, which it
/// does not replace.
#[test]
fn a_region_lists_its_outer_shell_first() {
    let mut body = BRep::new(Tolerance::default());
    let (region, outer) = body.add_solid_region(FaceHandle::default());
    let cavity = body.add_region_shell(region, FaceHandle::default());
    assert_eq!(body.region_outer_shell(region), Some(outer));
    assert_eq!(body.region_shells(region), vec![outer, cavity]);
    assert_eq!(body.shell_region(cavity), Some(region));
}
