//! PK-328 row 481 gates: `max_safe_offset` is interactive (<100 ms per
//! query) on every corpus body, names its binding constraint and
//! entity, and is SLIDER-CONSISTENT — the full operation serves at the
//! reported maximum and declines or trims at 1.05× it.

use std::time::{Duration, Instant};

use popkern_core::tolerance::Tolerance;
use popkern_examples::offset_corpus::corpus;
use popkern_shell::offset::offset_solid;
use popkern_shell::{
    OffsetConstraint, ShellDirection, ShellOptions, ShellWarning, max_safe_offset, shell,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Row 481 timing: every query on every corpus body, both directions,
/// answers inside the interactive budget. The budget is generous
/// against scheduler noise; the design target is 100 ms.
#[test]
fn max_safe_offset_is_interactive_on_the_corpus() {
    let budget = Duration::from_millis(100);
    for f in corpus(tol()).expect("corpus builds") {
        for direction in [ShellDirection::Inward, ShellDirection::Outward] {
            let t0 = Instant::now();
            let report = max_safe_offset(&f.brep, direction, tol());
            let elapsed = t0.elapsed();
            assert!(
                elapsed < budget,
                "{}/{direction:?}: {elapsed:?} exceeds the interactive budget",
                f.name
            );
            match report {
                Ok(rep) => {
                    assert!(
                        rep.max_distance > 0.0,
                        "{}/{direction:?}: non-positive max {rep}",
                        f.name
                    );
                    if rep.max_distance.is_finite() {
                        // A finite bound must name its mechanism+entity.
                        assert!(
                            !matches!(rep.binding, OffsetConstraint::Unbounded),
                            "{}/{direction:?}: finite bound with Unbounded binding",
                            f.name
                        );
                        let text = format!("{rep}");
                        assert!(
                            text.contains("maximum safe offset"),
                            "{}: display: {text}",
                            f.name
                        );
                    }
                }
                Err(e) => {
                    // The open sheet has no resolvable material side —
                    // a structured decline is the honest answer there.
                    assert!(
                        !f.watertight,
                        "{}/{direction:?}: closed body declined: {e}",
                        f.name
                    );
                }
            }
        }
    }
}

/// The operation a slider drives for each direction.
fn run_operation(
    brep: &popkern_core::topology::BRep,
    direction: ShellDirection,
    distance: f64,
) -> Result<popkern_shell::ShellResult, popkern_shell::ShellError> {
    match direction {
        ShellDirection::Inward => shell(
            brep,
            distance,
            ShellDirection::Inward,
            &[],
            &[],
            &ShellOptions::default(),
            tol(),
            None,
        ),
        _ => offset_solid(brep, distance, &ShellOptions::default(), tol()),
    }
}

/// True when a served result crossed the plain-offset boundary: the
/// trim stage removed material or a feature collapsed.
fn served_past_boundary(warnings: &[ShellWarning]) -> bool {
    warnings.iter().any(|w| {
        matches!(
            w,
            ShellWarning::RegionTrimmed { .. } | ShellWarning::FaceCollapsed { .. }
        )
    })
}

/// Row 481 slider consistency: for every closed corpus body and
/// direction with a finite bound, the operation SERVES at the reported
/// max and declines or trims at 1.05× it.
///
/// Recorded containment: fixture×direction pairs whose join stage
/// downgrades to `JoinStyle::Intersection` (double-seam closed-rim
/// anchors, PK-344) may TYPED-decline (`JoinFailed`) at the geometric
/// max instead of serving — never wrong geometry. The set is pinned
/// EXACTLY so a new decline pair fails this gate.
#[test]
fn slider_consistency_serve_at_max_decline_or_trim_past_it() {
    use popkern_shell::ShellError;
    let containment_declines = [
        ("bracket_bore_countersink", ShellDirection::Outward),
        // Same class: the rod's cylinder×sphere junction circle is a
        // double-seam closed-rim anchor; its skin downgrades and the
        // Intersection fallback cannot span the junction gap at the
        // (large) cylinder-vanish max.
        ("sphere_capped_rod", ShellDirection::Inward),
    ];
    let mut checked = 0usize;
    let mut declined: Vec<(String, ShellDirection)> = Vec::new();
    for f in corpus(tol()).expect("corpus builds") {
        if !f.watertight {
            continue;
        }
        for direction in [ShellDirection::Inward, ShellDirection::Outward] {
            let Ok(report) = max_safe_offset(&f.brep, direction, tol()) else {
                continue;
            };
            if !report.max_distance.is_finite() {
                continue;
            }
            let at_max = run_operation(&f.brep, direction, report.max_distance);
            match &at_max {
                Ok(result) => assert!(
                    !served_past_boundary(&result.warnings),
                    "{}/{direction:?}: served at max {} but past the plain-offset \
                     boundary (trim/collapse fired): {report}",
                    f.name,
                    report.max_distance
                ),
                Err(ShellError::JoinFailed { .. }) => {
                    declined.push((f.name.to_string(), direction));
                }
                Err(e) => panic!(
                    "{}/{direction:?}: must serve at reported max {} ({report}); got {e}",
                    f.name, report.max_distance
                ),
            }
            let past = run_operation(&f.brep, direction, report.max_distance * 1.05);
            match past {
                Err(_) => {}
                Ok(result) => assert!(
                    served_past_boundary(&result.warnings),
                    "{}/{direction:?}: served plain at 1.05×max ({}) — the reported \
                     bound is not the slider boundary: {report}",
                    f.name,
                    report.max_distance * 1.05
                ),
            }
            checked += 1;
        }
    }
    let expected: Vec<(String, ShellDirection)> = containment_declines
        .iter()
        .map(|(n, d)| (n.to_string(), *d))
        .collect();
    assert_eq!(
        declined, expected,
        "typed-decline pairs must match the recorded containment set (PK-344)"
    );
    assert!(
        checked >= 6,
        "slider pin must exercise a real spread of fixtures, got {checked}"
    );
}
