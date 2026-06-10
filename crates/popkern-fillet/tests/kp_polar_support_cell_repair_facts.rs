#![allow(dead_code)]

#[derive(Debug, Clone, Copy)]
struct CandidateFact<'a> {
    candidate: &'a str,
    after_worst_g0_gap: f64,
}

#[derive(Debug, Clone, Copy)]
struct RepairFact<'a> {
    selected_candidate: &'a str,
    after_worst_g0_gap: f64,
}

pub(crate) fn assert_selected_repair_facts_are_consistent(
    fixture: &str,
    support_cell_details: &[String],
) {
    let selected_candidates = support_cell_details
        .iter()
        .filter_map(|detail| {
            detail
                .contains("candidate_status=accepted_selected")
                .then(|| candidate_fact(detail))
                .flatten()
        })
        .collect::<Vec<_>>();
    let repair_facts = support_cell_details
        .iter()
        .filter_map(|detail| {
            detail
                .contains("repair_status=accepted_non_regressing_final_repair")
                .then(|| repair_fact(detail))
                .flatten()
        })
        .collect::<Vec<_>>();

    assert!(
        !selected_candidates.is_empty(),
        "{fixture} should emit selected candidate facts; support_cell_details={support_cell_details:?}"
    );
    assert!(
        !repair_facts.is_empty(),
        "{fixture} should emit final repair facts; support_cell_details={support_cell_details:?}"
    );
    for repair in &repair_facts {
        assert!(
            selected_candidates.iter().any(|candidate| {
                candidate.candidate == repair.selected_candidate
                    && same_metric(candidate.after_worst_g0_gap, repair.after_worst_g0_gap)
            }),
            "{fixture} final repair should match a selected candidate score; repair={repair:?} selected_candidates={selected_candidates:?}"
        );
    }
    assert_final_quality_matches_last_repair(fixture, support_cell_details, &repair_facts);
}

fn assert_final_quality_matches_last_repair(
    fixture: &str,
    support_cell_details: &[String],
    repair_facts: &[RepairFact<'_>],
) {
    let Some(final_quality_g0) = support_cell_details
        .iter()
        .rev()
        .find(|detail| detail.contains("measurement_scope=final_post_pipeline"))
        .and_then(|detail| parse_f64_field(detail, "worst_g0_gap"))
    else {
        panic!(
            "{fixture} should emit final support-cell shared-trim quality; support_cell_details={support_cell_details:?}"
        );
    };
    let last_repair = repair_facts
        .last()
        .expect("repair facts checked nonempty before final quality correlation");
    assert!(
        same_metric(final_quality_g0, last_repair.after_worst_g0_gap),
        "{fixture} final quality G0 should match the last accepted repair; final_g0={final_quality_g0:.15e} last_repair={last_repair:?}"
    );
}

fn candidate_fact(detail: &str) -> Option<CandidateFact<'_>> {
    Some(CandidateFact {
        candidate: parse_str_field(detail, "candidate")?,
        after_worst_g0_gap: parse_f64_field(detail, "after_worst_g0_gap")?,
    })
}

fn repair_fact(detail: &str) -> Option<RepairFact<'_>> {
    Some(RepairFact {
        selected_candidate: parse_str_field(detail, "selected_candidate")?,
        after_worst_g0_gap: parse_f64_field(detail, "after_worst_g0_gap")?,
    })
}

fn parse_str_field<'a>(detail: &'a str, key: &str) -> Option<&'a str> {
    let prefix = format!("{key}=");
    detail
        .split(';')
        .find_map(|field| field.strip_prefix(&prefix))
}

fn parse_f64_field(detail: &str, key: &str) -> Option<f64> {
    parse_str_field(detail, key)?.parse().ok()
}

fn same_metric(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-12_f64.max(a.abs().max(b.abs()) * 1.0e-9)
}
