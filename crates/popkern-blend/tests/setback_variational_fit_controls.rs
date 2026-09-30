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

/// E6 negative control: the curvature-weighted fit folds on the steep
/// octagonal apex (measured regularity 0: a flipped normal inside a
/// sector). The fold gate declines it typed, before certification, and
/// the tangent fit serves. (The square apex h0.8 this control first
/// used sat between 5e-3 and the old 2e-2 floor, not folded; since
/// PK-594 it clears the gate and the rim curvature witness declines it.)
#[test]
fn e6_a_folded_fit_declines_at_the_fold_gate() {
    let radii = [0.02, 0.04, 0.06, 0.04, 0.02, 0.05, 0.03, 0.06];
    let result =
        with_kp20_cap_forced_to_decline(|| pyramid_apex(8, 1.2, &radii, ContinuityLevel::G2));
    let record = only_record("pyramid8 h1.2", &result);
    let first = &record.attempts[0];
    assert_eq!(first.mode, SetbackFitMode::CurvatureWeighted);
    let SetbackFitOutcome::Declined(reason) = &first.outcome else {
        panic!("the folded fit must decline before certification: {first:#?}");
    };
    assert!(
        reason.contains("folds") && reason.contains("min regularity 0.000e0"),
        "{reason}"
    );
    // The tangent cap is thin (served regularity 8.4e-3, under the
    // PK-574 fixture set's 0.1 ratchet) and certified G1 on every seam.
    assert_eq!(record.served, Some(ContinuityLevel::G1), "{record:#?}");
    let last = record.attempts.last().expect("an attempt");
    assert_eq!(last.mode, SetbackFitMode::Tangent, "{record:#?}");
    let SetbackFitOutcome::Served(certified) = &last.outcome else {
        panic!("the tangent fit serves: {record:#?}");
    };
    assert_eq!(certified.level, ContinuityLevel::G1);
    assert!(result.brep.validate().is_ok());
}
