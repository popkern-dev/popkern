use super::*;
use crate::classify::classify_faces;
use crate::inner_skin::build_inner_skin;
use crate::surface_offset::compute_offset_surfaces;
use crate::types::ShellDirection;
use popkern_core::tolerance::Tolerance;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn plan_for(brep: &BRep, direction: ShellDirection, thickness: f64, style: JoinStyle) -> JoinPlan {
    let classification = classify_faces(brep, thickness, &[], direction, &[]).unwrap();
    let offset = compute_offset_surfaces(brep, &classification, true, &tol()).unwrap();
    let (skin, skin_map) = build_inner_skin(brep, &[]).unwrap();
    classify_joins(&skin, &skin_map, &offset.surfaces, &tol(), style).unwrap()
}

#[test]
fn cube_outward_classifies_all_twelve_edges_as_gap() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let plan = plan_for(&brep, ShellDirection::Outward, 0.1, JoinStyle::Arc);
    assert_eq!(plan.gap_edges.len(), 12, "convex cube edges under outward");
    assert_eq!(plan.effective, JoinStyle::Arc);
    assert_eq!(plan.gap_corners.len(), 8, "all-gap trihedral corners");
    assert!(plan.warnings.is_empty());
    for gap in plan.gap_edges.values() {
        assert!((gap.ball_radius - 0.1).abs() < 1e-9, "{}", gap.ball_radius);
    }
}

#[test]
fn cube_inward_has_no_gap_edges() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let plan = plan_for(&brep, ShellDirection::Inward, 0.1, JoinStyle::Arc);
    assert!(
        plan.gap_edges.is_empty(),
        "erosion of a convex body keeps sharp edges: {:?}",
        plan.gap_edges
    );
    assert!(plan.warnings.is_empty());
}

#[test]
fn l_bracket_outward_downgrades_on_mixed_corners() {
    // L-shaped prism: one concave (reflex) vertical edge. Outward
    // offset: convex edges gap, the reflex edge is resolvable — its
    // endpoints are mixed corners, an unserved vertex class (PK-326).
    let profile = [
        popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
        popkern_core::geometry::Point3::new(2.0, 0.0, 0.0),
        popkern_core::geometry::Point3::new(2.0, 1.0, 0.0),
        popkern_core::geometry::Point3::new(1.0, 1.0, 0.0),
        popkern_core::geometry::Point3::new(1.0, 2.0, 0.0),
        popkern_core::geometry::Point3::new(0.0, 2.0, 0.0),
    ];
    let brep = BRep::prism(
        &profile,
        popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        1.0,
        tol(),
    )
    .unwrap();
    let plan = plan_for(&brep, ShellDirection::Outward, 0.1, JoinStyle::Arc);
    assert!(!plan.gap_edges.is_empty());
    assert_eq!(
        plan.effective,
        JoinStyle::Intersection,
        "mixed corners must contain to the sharp join"
    );
    assert!(
        plan.warnings
            .iter()
            .any(|w| matches!(w, ShellWarning::JoinStyleDowngraded { .. })),
        "downgrade must be named"
    );
    assert!(plan.gap_corners.is_empty());
}

#[test]
fn l_bracket_inward_gap_sits_on_the_reflex_edge() {
    // Erosion rounds reflex features: inward, the single concave edge
    // is the gap; its endpoints are mixed corners → downgrade.
    let profile = [
        popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
        popkern_core::geometry::Point3::new(2.0, 0.0, 0.0),
        popkern_core::geometry::Point3::new(2.0, 1.0, 0.0),
        popkern_core::geometry::Point3::new(1.0, 1.0, 0.0),
        popkern_core::geometry::Point3::new(1.0, 2.0, 0.0),
        popkern_core::geometry::Point3::new(0.0, 2.0, 0.0),
    ];
    let brep = BRep::prism(
        &profile,
        popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        1.0,
        tol(),
    )
    .unwrap();
    let plan = plan_for(&brep, ShellDirection::Inward, 0.1, JoinStyle::Arc);
    assert_eq!(
        plan.gap_edges.len(),
        1,
        "exactly the reflex vertical edge rounds under erosion"
    );
    assert_eq!(plan.effective, JoinStyle::Intersection);
}

#[test]
fn intersection_style_skips_corner_planning_but_keeps_gap_edges() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let plan = plan_for(&brep, ShellDirection::Outward, 0.1, JoinStyle::Intersection);
    assert_eq!(plan.gap_edges.len(), 12);
    assert_eq!(plan.effective, JoinStyle::Intersection);
    assert!(plan.gap_corners.is_empty());
    assert!(plan.arc_skip_edges().is_empty());
}

#[test]
fn cone_outward_classifies_base_rim_as_gap() {
    let brep = BRep::cone(1.0, 2.0, tol()).unwrap();
    let plan = plan_for(&brep, ShellDirection::Outward, 0.15, JoinStyle::Arc);
    assert_eq!(
        plan.gap_edges.len(),
        1,
        "the closed base rim is the cone's single gap edge: {:?}",
        plan.gap_edges
    );
}

#[test]
fn cone_outward_warns_sharp_apex_under_arc() {
    let brep = BRep::cone(1.0, 2.0, tol()).unwrap();
    let plan = plan_for(&brep, ShellDirection::Outward, 0.15, JoinStyle::Arc);
    assert!(
        plan.warnings
            .iter()
            .any(|w| matches!(w, ShellWarning::ConeApexSharp { .. })),
        "dilated cone apex deviation must be recorded: {:?}",
        plan.warnings
    );
    // Erosion keeps the apex exact — no warning.
    let inward = plan_for(&brep, ShellDirection::Inward, 0.15, JoinStyle::Arc);
    assert!(
        !inward
            .warnings
            .iter()
            .any(|w| matches!(w, ShellWarning::ConeApexSharp { .. })),
        "eroded apex is exact: {:?}",
        inward.warnings
    );
}
