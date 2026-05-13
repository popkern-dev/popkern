use crate::a2_validate_helpers::{a2_fallback_reasons, build_under_a2};
use crate::helpers::measure_cap_outer_boundary;

/// Native cube-mixed A2 Gregory cap-outer gate. This is the first
/// repeated-apex Gregory fixture promoted out of containment: it must
/// complete natively and satisfy the measured cap/support seam quality
/// bar without relying on the legacy fallback route.
#[test]
fn a2_5_cap_outer_g0_cube_mixed_under_a2() {
    let result = build_under_a2(false, &[0.08, 0.10, 0.12]);
    let reasons = a2_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "cube_mixed cap-outer gate must complete through native A2; got:\n{}",
        reasons.join("\n"),
    );
    let cap_faces = result
        .corner_patches
        .iter()
        .map(|p| p.face)
        .collect::<Vec<_>>();
    let m = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
    eprintln!(
        "[A2.5 cap_outer cube_mixed] g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured,
    );
    assert!(m.edges_measured >= 1, "cube_mixed measured no cap seams");
    assert!(
        m.worst_g0_gap <= 5.0e-3,
        "cube_mixed cap_outer_g0 {:.3e} exceeds 5.0e-3",
        m.worst_g0_gap,
    );
    assert!(
        m.worst_g1_coplanarity <= 5.0e-2,
        "cube_mixed cap_outer_g1 {:.3e} exceeds 5.0e-2",
        m.worst_g1_coplanarity,
    );
}
