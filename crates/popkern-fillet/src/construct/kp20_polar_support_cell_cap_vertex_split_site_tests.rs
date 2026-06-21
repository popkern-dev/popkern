use super::*;

#[test]
fn within_cap_decision_routes_through_use_candidate() {
    // Displacement comfortably under the cap → caller proceeds.
    assert_eq!(
        decide_from_displacement(0.5, 0.05, 0.1),
        SplitParameterDecision::UseCandidate(0.5)
    );
}

#[test]
fn at_cap_boundary_decision_still_uses_candidate() {
    // Equality is accepted — `displacement <= cap`. Mirrors the
    // boundary-fit's bound test at `MAX_BOUNDARY_ROW_DISPLACEMENT`.
    assert_eq!(
        decide_from_displacement(0.42, 0.1, 0.1),
        SplitParameterDecision::UseCandidate(0.42)
    );
}

#[test]
fn above_cap_decision_routes_through_skip_unreachable() {
    // Displacement above cap → caller falls through to the
    // completion ladder's next candidate (Stage B1 default; Stage B2
    // refines into UseProjected when search succeeds).
    assert_eq!(
        decide_from_displacement(0.5, 0.3, 0.1),
        SplitParameterDecision::SkipUnreachable {
            candidate_t: 0.5,
            best_displacement: 0.3,
        }
    );
}

#[test]
fn skip_unreachable_carries_observed_displacement() {
    // The observed displacement is preserved so downstream
    // diagnostics can record the failure magnitude.
    let SplitParameterDecision::SkipUnreachable {
        candidate_t,
        best_displacement,
    } = decide_from_displacement(0.25, 0.817, 0.1)
    else {
        panic!("expected SkipUnreachable for above-cap displacement");
    };
    assert_eq!(candidate_t, 0.25);
    assert!((best_displacement - 0.817).abs() < 1.0e-12);
}
