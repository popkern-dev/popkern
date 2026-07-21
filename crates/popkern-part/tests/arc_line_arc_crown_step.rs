//! PK-155: STEP fidelity of an extruded profile whose crown is an
//! arc-line-arc run (two r=2 corner arcs flanking a short 1 mm flat —
//! the peg-lip blade-top class). The mixed run compiles to one
//! rational degree-2 side-wall surface; the exported STEP entity for
//! that surface was rejected by OCC's reader ("Entity has illegal
//! type / Surface has not been created"), dropping the face.
//!
//! This test generates the fixture and holds the popkern-side gates
//! (build, validate, STEP round-trip). The OCC acceptance leg runs
//! via tools/occ-oracle/occ_check.py on the emitted file.

use popkern_part::{Plane, Profile};

#[test]
fn arc_line_arc_crown_round_trips() {
    let r = 2.0_f64;
    let d = r - r / 2.0_f64.sqrt();
    // Cross-section: 5 wide, crown at y=10, 1 mm flat between the
    // two corner arcs (flat from x=2 to x=3).
    let part = Profile::on(Plane::XY)
        .path(0.0, 0.0)
        .line_to(5.0, 0.0)
        .line_to(5.0, 8.0)
        .arc_to((5.0 - d, 10.0 - d), (3.0, 10.0))
        .line_to(2.0, 10.0)
        .arc_to((d, 10.0 - d), (0.0, 8.0))
        .close()
        .extrude(4.0)
        .expect("arc-line-arc crown profile must extrude");

    part.brep().validate().expect("crown solid must validate");

    let dir = std::path::Path::new("../../examples/output/fixtures");
    std::fs::create_dir_all(dir).expect("fixture output dir");
    let path = dir.join("arc_line_arc_crown.step");
    part.export_step(&path).expect("crown solid must export");
    part.assert_step_roundtrip(&path)
        .expect("crown solid must round-trip through STEP");
}
