//! PK-393 — subtraction of a cylinder flush with the host's faces
//! serves a valid 7-face body in every flush configuration.
//!
//! The tool cap's rim circle lies ON the host face; its UV clip on the
//! tool lateral classified all-outside (3 seam-degenerate crossings)
//! and the imprint curve was dropped: the host face kept no hole and
//! the wall self-intersected it. A closed record whose intervals all
//! classify outside while hugging the trim boundary is now kept as the
//! flush rim it is.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_part::prelude::*;

fn check(label: &str, tool_z0: f64, tool_h: f64) {
    let plate = Profile::on(Plane::XY)
        .rect(2.4, 1.6)
        .extrude(0.4)
        .expect("plate");
    let tool = Profile::on(Plane::xy_at(tool_z0))
        .circle(0.35)
        .extrude(tool_h)
        .expect("tool");
    let cut = plate.subtract(&tool).expect("subtract serves");
    assert_eq!(cut.face_count(), 7, "{label}: plate 6 + one hole wall");
    let issues = cut.brep().validate_all();
    assert!(issues.is_empty(), "{label}: {issues:?}");
    assert_eq!(
        cut.open_edge_count().expect("mesh"),
        0,
        "{label}: watertight"
    );
    let b = cut.brep();
    // Both the top and bottom faces carry the hole (two loops each).
    for z in [0.0, 0.4] {
        let face = b
            .faces
            .iter()
            .find(|(h, f)| {
                matches!(
                    f.surface.as_ref().map(|s| &s.geometry),
                    Some(SurfaceGeometry::Plane { .. })
                ) && {
                    let pts = b.face_boundary_sample_points(*h, 4);
                    !pts.is_empty() && pts.iter().all(|p| (p.z - z).abs() < 1e-9)
                }
            })
            .map(|(h, _)| h)
            .unwrap_or_else(|| panic!("{label}: face at z={z}"));
        let loops = b.loops.values().filter(|l| l.face == face).count();
        assert_eq!(loops, 2, "{label}: face at z={z} carries the hole rim");
    }
}

#[test]
fn flush_cut_matrix_serves_valid_bodies() {
    check("flush both", 0.0, 0.4);
    check("flush bottom, proud top", 0.0, 0.5);
    check("proud bottom, flush top", -0.1, 0.5);
    check("proud both", -0.1, 0.6);
}
