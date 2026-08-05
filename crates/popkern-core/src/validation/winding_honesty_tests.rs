//! Calibration for the winding-flag honesty audit: primitive
//! constructors are winding-flag consistent (the plane-winding
//! invariant), so every determinate verdict must be honest — and
//! flipping one stored flag must flip exactly that face's verdict.

use crate::tolerance::Tolerance;
use crate::topology::BRep;

fn assert_all_honest(brep: &BRep, label: &str) {
    let audit = brep.audit_winding_honesty();
    assert!(!audit.is_empty(), "{label}: audit produced no verdicts");
    for v in &audit {
        assert_ne!(
            v.honest,
            Some(false),
            "{label}: face {:?} dishonest (area {:.3e}, flag {}, reason {:?})",
            v.face,
            v.signed_uv_area,
            v.orientation,
            v.reason
        );
    }
    let determinate = audit.iter().filter(|v| v.honest.is_some()).count();
    assert!(
        determinate * 2 >= audit.len(),
        "{label}: fewer than half the faces produced determinate verdicts \
         ({determinate}/{})",
        audit.len()
    );
}

#[test]
fn unit_cube_is_winding_honest() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    assert_all_honest(&cube, "unit_cube");
    assert!(cube.dishonest_winding_faces().is_empty());
}

#[test]
fn cylinder_is_winding_honest_across_periodic_seam() {
    let cylinder = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    assert_all_honest(&cylinder, "cylinder");
}

#[test]
fn sphere_wall_seam_spine_is_indeterminate_not_dishonest() {
    // A full sphere wall's outer loop is a seam spine: net-zero
    // enclosed UV area. The audit must decline the verdict rather
    // than amplify roundoff into "dishonest".
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    assert!(sphere.dishonest_winding_faces().is_empty());
    let audit = sphere.audit_winding_honesty();
    assert!(
        audit.iter().any(|v| v.honest.is_none()
            && matches!(v.reason, "degenerate-uv-area" | "pole-touching-loop")),
        "sphere wall must decline with a chart-degeneracy reason, got {:?}",
        audit
            .iter()
            .map(|v| (v.honest, v.reason))
            .collect::<Vec<_>>()
    );
}

#[test]
fn flipping_one_flag_is_detected_as_dishonest() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let victim = cube
        .audit_winding_honesty()
        .into_iter()
        .find(|v| v.honest == Some(true))
        .map(|v| v.face)
        .expect("cube has a determinate honest face");
    let sof = cube
        .faces
        .get_mut(victim)
        .and_then(|f| f.surface.as_mut())
        .expect("victim face has surface geometry");
    sof.orientation = !sof.orientation;
    let dishonest = cube.dishonest_winding_faces();
    assert_eq!(
        dishonest,
        vec![victim],
        "exactly the flipped face must report dishonest"
    );
}
