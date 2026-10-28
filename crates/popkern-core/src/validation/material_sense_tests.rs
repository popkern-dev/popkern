//! Material-sense resolution pins.

use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn unit_cube_material_sense_is_determinate_and_flips_with_body_orientation() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let seed = brep.faces.keys().next().unwrap();
    let sense = brep
        .component_material_sense(seed)
        .expect("cube sense determinate");

    let mut reversed = brep.clone();
    reversed.reverse_body_orientation();
    let seed_r = reversed.faces.keys().next().unwrap();
    let sense_r = reversed
        .component_material_sense(seed_r)
        .expect("reversed cube sense determinate");
    assert_ne!(
        sense, sense_r,
        "reversing the body orientation must flip the material sense"
    );
}

/// The sense does not depend on where the body sits: a thin bore far
/// from the origin reads as it does at the origin (the coarse grid's
/// residual used to be multiplied by the placement offset).
#[test]
fn material_sense_is_placement_invariant_for_a_thin_bore() {
    use crate::geometry::transform::Transform3;
    let tol = Tolerance::default();
    let bore = BRep::cylinder(0.04, 0.28, tol).unwrap();
    let seed = bore.faces.keys().next().unwrap();
    let home = bore.component_material_sense(seed).expect("determinate");
    let mut flipped = Vec::new();
    for i in -8..=8 {
        for j in -8..=8 {
            let (x, y) = (0.61 * f64::from(i) / 4.0, 0.31 * f64::from(j) / 4.0);
            let placed = bore
                .transform(&Transform3::translation(x, y, -0.1))
                .unwrap();
            let seed = placed.faces.keys().next().unwrap();
            if placed.component_material_sense(seed) != Some(home) {
                flipped.push((x, y));
            }
        }
    }
    assert!(flipped.is_empty(), "sense changed at {flipped:?}");
}
