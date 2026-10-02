//! PK-574 row 878 and E6: negative controls for the variational corner
//! fit — a corrupted cap is rejected on its certificate, and a folded
//! fit declines at the fold gate before certification.

#[path = "support/setback_fit_fixtures.rs"]
mod setback_fit_fixtures;

use popkern_blend::dispatch_overrides::{with_fit_cap_rim_tilted, with_kp20_cap_forced_to_decline};
use popkern_blend::{
    ContinuityLevel, CornerMethod, SetbackFitBorder, SetbackFitMode, SetbackFitOutcome,
};
use setback_fit_fixtures::*;

/// Row 878 negative control: a cap whose cross-boundary tangent is
/// corrupted after the fit (rim position intact) must not be served.
/// The certifier reads the committed topology, not the fit's report —
/// which still says the targets were met — so every attempt is
/// rejected on its certificate with a G1 bound past the gate, and the
/// corner falls through to the vertex-blend path.
#[test]
fn row_878_a_tilted_cap_is_rejected_on_its_certificate() {
    let result = with_fit_cap_rim_tilted(|| cube_corner(&[0.05, 0.2, 0.3]));
    let record = only_record("tilted cube", &result);
    assert_eq!(record.served, None, "{record:#?}");
    // Curvature-weighted, tangent, and the corner-cell tangent rescue,
    // on the corners border and then on the bridged border (PK-618), then
    // the two tangent modes on the bowed-bridge border (PK-619).
    assert_eq!(record.attempts.len(), 8, "{record:#?}");
    for attempt in &record.attempts {
        let SetbackFitOutcome::Rejected { certified, reason } = &attempt.outcome else {
            panic!("tilted cap must be rejected on its certificate: {attempt:#?}");
        };
        assert_eq!(certified.level, ContinuityLevel::G0, "{reason}");
        assert!(
            certified.rim[0] <= G0_GATE,
            "rim position intact: {:?}",
            certified.rim
        );
        assert!(
            certified.rim[1] > G1_GATE,
            "tilt measured: {:?}",
            certified.rim
        );
    }
    assert!(
        result
            .corner_patches
            .iter()
            .all(|p| p.method != CornerMethod::PlateMin),
        "no fitted sector survives the rejection"
    );
    result
        .brep
        .validate()
        .expect("the fallthrough body validates");
}

/// E6 negative control: on the steep variable square apex (h3, radii
/// 0.02/0.04/0.06/0.04, KP20 forced off) the corners border's fits fold
/// or pinch: the curvature-weighted fit at regularity 1.0e-2, the
/// tangent fit at 1.8e-2, both under the 2e-2 served floor and bending
/// past the fairness bound. The fold gate declines each typed, before
/// certification. The bridged border (PK-618) then serves the corner,
/// so the gate's decline is read from the record, not a fallthrough.
#[test]
fn e6_a_folded_fit_declines_at_the_fold_gate() {
    let body = popkern_core::topology::BRep::n_gonal_pyramid(
        4,
        1.0,
        3.0,
        popkern_core::tolerance::Tolerance::default(),
    )
    .expect("pyramid");
    let apex = body
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(v, _)| v)
        .expect("apex");
    let result = with_kp20_cap_forced_to_decline(|| {
        fillet(&body, apex, &[0.02, 0.04, 0.06, 0.04], ContinuityLevel::G2)
    });
    let record = only_record("steep square", &result);
    for mode in [SetbackFitMode::CurvatureWeighted, SetbackFitMode::Tangent] {
        let attempt = record
            .attempts
            .iter()
            .find(|a| a.mode == mode && a.border == SetbackFitBorder::Corners)
            .unwrap_or_else(|| panic!("{mode:?} on the corners border: {record:#?}"));
        let SetbackFitOutcome::Declined(reason) = &attempt.outcome else {
            panic!("{mode:?} declines before certification: {attempt:#?}");
        };
        assert!(
            reason.starts_with("variational fit: variational cap: sector")
                && (reason.contains("folds near") || reason.contains("is thin near")),
            "{mode:?} declines at the fold gate: {reason}"
        );
    }
}
