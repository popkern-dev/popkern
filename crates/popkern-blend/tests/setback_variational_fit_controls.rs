//! PK-574 row 878 and E6: negative controls for the variational corner
//! fit — a corrupted cap is rejected on its certificate, and a folded
//! fit declines at the fold gate before certification.

#[path = "support/setback_fit_fixtures.rs"]
mod setback_fit_fixtures;

use popkern_blend::dispatch_overrides::{with_fit_cap_rim_tilted, with_kp20_cap_forced_to_decline};
use popkern_blend::{ContinuityLevel, CornerMethod, SetbackFitMode, SetbackFitOutcome};
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
    assert_eq!(record.attempts.len(), 2, "{record:#?}");
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

/// E6 negative control: the curvature-weighted fit folds on the square
/// apex (measured regularity 0 inside a sector). The fold gate declines
/// it typed, before certification, and the tangent fit serves.
#[test]
fn e6_a_folded_fit_declines_at_the_fold_gate() {
    let result =
        with_kp20_cap_forced_to_decline(|| pyramid_apex(4, 0.8, &[0.16], ContinuityLevel::G2));
    let record = only_record("pyramid4 h0.8", &result);
    let first = &record.attempts[0];
    assert_eq!(first.mode, SetbackFitMode::CurvatureWeighted);
    let SetbackFitOutcome::Declined(reason) = &first.outcome else {
        panic!("the folded fit must decline before certification: {first:#?}");
    };
    assert!(reason.contains("folds"), "{reason}");
    assert_record_serves(
        "pyramid4 h0.8",
        &record,
        SetbackFitMode::Tangent,
        ContinuityLevel::G1,
    );
}
