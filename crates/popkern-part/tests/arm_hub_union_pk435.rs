//! PK-435 repro: union a straight-sided arm prism into a round hub.
//! The arm's two long edges cross the hub's cylindrical wall obliquely
//! and its inner end sits inside the hub. Before the fix this declined
//! with "face split failed: ... trim curve end point does not match any
//! vertex on the face boundary" — the arm-side / hub-wall crossing
//! vertices were never minted on the cylinder lateral face boundary.
use popkern_part::prelude::*;

/// Hub dimensions mirror `planetary_gearbox.rhai`'s spider hub.
const HUB_R: f64 = 0.24;
const HUB_Z: f64 = 0.44;
const ARM_HW: f64 = 0.09; // arm half-width
const ARM_OUT: f64 = 0.9; // arm reach

fn hub() -> Part {
    Profile::on(Plane::xy_at(0.0))
        .circle(HUB_R)
        .extrude(HUB_Z)
        .expect("hub cylinder")
}

/// One radial arm along +x: a rectangle from inside the hub (t0) out to
/// ARM_OUT, half-width ARM_HW. Its two long edges cross the hub wall.
fn arm() -> Part {
    let t0 = (HUB_R * HUB_R - ARM_HW * ARM_HW).sqrt();
    Profile::on(Plane::xy_at(0.0))
        .polygon(&[
            (t0 - 0.05, -ARM_HW),
            (ARM_OUT, -ARM_HW),
            (ARM_OUT, ARM_HW),
            (t0 - 0.05, ARM_HW),
        ])
        .extrude(HUB_Z)
        .expect("arm prism")
}

/// Same arm rotated to +y (90°): its wall crossings land near 68°..112°,
/// well clear of the cylinder seam at angle 0°. Isolates the oblique
/// straight-edge-through-round-wall crossing from any seam-straddle.
fn arm_plus_y() -> Part {
    let t0 = (HUB_R * HUB_R - ARM_HW * ARM_HW).sqrt();
    Profile::on(Plane::xy_at(0.0))
        .polygon(&[
            (ARM_HW, t0 - 0.05),
            (ARM_HW, ARM_OUT),
            (-ARM_HW, ARM_OUT),
            (-ARM_HW, t0 - 0.05),
        ])
        .extrude(HUB_Z)
        .expect("arm prism +y")
}

/// Mid-height arm (z inside the hub, NOT flush with the caps): the
/// PK-150 prism-through-wall case, which already serves. Confirms the
/// PK-435 trigger is specifically the flush coincident caps.
fn arm_mid_height() -> Part {
    let t0 = (HUB_R * HUB_R - ARM_HW * ARM_HW).sqrt();
    Profile::on(Plane::xy_at(0.10))
        .polygon(&[
            (t0 - 0.05, -ARM_HW),
            (ARM_OUT, -ARM_HW),
            (ARM_OUT, ARM_HW),
            (t0 - 0.05, ARM_HW),
        ])
        .extrude(0.24)
        .expect("arm prism mid-height")
}

#[test]
fn arm_mid_height_into_hub_union_serves_watertight() {
    let out = hub()
        .union(&arm_mid_height())
        .expect("mid-height arm ∪ hub union serves");
    let issues = out.brep().validate_all();
    assert!(issues.is_empty(), "validate_all: {issues:?}");
    assert_eq!(out.open_edge_count().expect("mesh"), 0, "watertight");
}

// PK-439 closed the flush-cap class at root: the 3D-Newell walk now
// measures true departure tangents (twin-owned curves borrowed, not
// chord-approximated), walks curved faces material-on-left like the
// planar convention, only accepts STRUCTURAL boundary dead-ends as
// half-split chains, and samples all-shared fragments in the
// fragment's own face chart. With the PK-435 coincidence fixes
// (material-side orientation over curve-sampled centroids +
// curve-aware partial-overlap classify), the union serves watertight.
#[test]
fn arm_into_hub_union_serves_watertight() {
    let out = hub().union(&arm()).expect("arm ∪ hub union serves");
    let issues = out.brep().validate_all();
    assert!(issues.is_empty(), "validate_all: {issues:?}");
    assert_eq!(out.open_edge_count().expect("mesh"), 0, "watertight");
}

#[test]
fn arm_off_seam_into_hub_union_serves_watertight() {
    let out = hub()
        .union(&arm_plus_y())
        .expect("arm(+y) ∪ hub union serves");
    let issues = out.brep().validate_all();
    assert!(issues.is_empty(), "validate_all: {issues:?}");
    assert_eq!(out.open_edge_count().expect("mesh"), 0, "watertight");
}
