use crate::a2_validate_helpers::{a2_fallback_reasons, build_under_a2};
use crate::helpers::measure_cap_outer_boundary;

/// Native cube-mixed A2 Gregory cap-outer status. This fixture must
/// complete natively and preserve cap/support G0 without relying on
/// the legacy fallback route. Cap/support G1 remains a recorded
/// construction-quality blocker; do not hide it with local support-frame
/// rewrites that break ring continuity or route ownership.
#[test]
fn a2_5_cap_outer_cube_mixed_records_remaining_g1_blocker() {
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
    assert!(m.edges_measured >= 1, "cube_mixed measured no cap seams");
    assert!(
        m.worst_g0_gap <= 5.0e-3,
        "cube_mixed cap_outer_g0 {:.3e} exceeds 5.0e-3",
        m.worst_g0_gap,
    );
    assert!(
        m.worst_g1_coplanarity > 5.0e-2,
        "cube_mixed should continue exposing the cap_outer_g1 blocker until construction-side support continuation is repaired; measured {:.3e}",
        m.worst_g1_coplanarity,
    );
}
