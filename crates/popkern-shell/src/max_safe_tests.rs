use super::*;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn cube_inward_binds_at_half_extent() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let report = max_safe_offset(&brep, ShellDirection::Inward, tol()).unwrap();
    // True bound 0.5 (opposite walls meet / edges collapse), reported
    // at the safe fraction.
    assert!(
        (report.max_distance - 0.5 * SAFE_FRACTION).abs() < 1e-6,
        "{report}"
    );
    assert!(
        matches!(
            report.binding,
            OffsetConstraint::EdgeCollapse { .. } | OffsetConstraint::ThinWall { .. }
        ),
        "{report}"
    );
}

#[test]
fn cube_outward_is_unbounded() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let report = max_safe_offset(&brep, ShellDirection::Outward, tol()).unwrap();
    assert!(report.max_distance.is_infinite(), "{report}");
    assert!(matches!(report.binding, OffsetConstraint::Unbounded));
}

#[test]
fn cube_midplane_doubles_the_inward_bound() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let report = max_safe_offset(&brep, ShellDirection::Midplane, tol()).unwrap();
    assert!(
        (report.max_distance - 1.0 * SAFE_FRACTION).abs() < 1e-6,
        "{report}"
    );
}

#[test]
fn sphere_inward_binds_at_radius() {
    let brep = BRep::sphere(1.0, tol()).unwrap();
    let report = max_safe_offset(&brep, ShellDirection::Inward, tol()).unwrap();
    // Unit sphere: surface vanish at r = 1.
    assert!(
        (report.max_distance - 1.0 * SAFE_FRACTION).abs() < 1e-6,
        "{report}"
    );
    assert!(
        matches!(report.binding, OffsetConstraint::SurfaceVanish { .. }),
        "{report}"
    );
}

#[test]
fn report_display_carries_value_and_entity() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let report = max_safe_offset(&brep, ShellDirection::Inward, tol()).unwrap();
    let text = format!("{report}");
    assert!(text.contains("maximum safe offset"), "{text}");
    assert!(
        text.contains("Edge") || text.contains("Face") || text.contains("face"),
        "binding entity must be named: {text}"
    );
}

#[test]
fn smallest_positive_root_arms() {
    // Pure linear: 1 − d = 0.
    assert!((smallest_positive_root(1.0, -1.0, 0.0).unwrap() - 1.0).abs() < 1e-12);
    // No positive root.
    assert!(smallest_positive_root(1.0, 1.0, 1.0).is_none());
    // Quadratic with two positive roots picks the smaller.
    let r = smallest_positive_root(2.0, -3.0, 1.0).unwrap();
    assert!((r - 1.0).abs() < 1e-12);
}
