//! PK-594: the wide-apex census rows that used to fail one strip with
//! `PartialFailure` and never reach a corner cap.
//!
//! Measured 2026-09-30 (G2 request, whole-blend evidence):
//! - Served since PK-591/593/596/597: pyr8 h0.4/h0.8 variable and wide.
//!   PK-618: pyr5 h1.2 wide and pyr8 h1.2 variable serve fits on the
//!   bridged border inside the fairness bound (their fans read 8.5 and
//!   15.3 / r_min). PK-613's faired guide serves
//!   pyr6 h1.2 variable (its unfaired fan pinched at 4.5e5 / r_min) and
//!   the wide octagons as G2 fans inside the fairness bound.
//! - Class B, pyr8 h1.2 variable and wide: the KP20 fan and the
//!   curvature-weighted fit folded; PK-594 served tangent fits at
//!   regularity 8.4e-3 / 1.6e-2 that PK-612 measured as near-fold
//!   pinches (2572 / 2850 x 1/r_min). With the faired guide (PK-613)
//!   the fan serves wide at G2. Variable: the corners border's fits
//!   still fold; the bridged border's tangent fit serves (PK-618).
//! - Class A, served by PK-611: pyr5 h0.4 variable and wide, pyr5 h0.8
//!   wide, pyr6 h0.4 wide. A thin strip's setback arc, run between the
//!   corners its wide neighbours fix, crosses its own trimmed boundary
//!   mid-arc, and the corners border's carve declines. The surgery
//!   retries on the bridged border (each strip at its own station, a
//!   chord bridge on every support face), where the tangent fit
//!   certifies G1 on every seam and passes `validate_all`.

use std::collections::BTreeSet;

use popkern_blend::blend_continuity_evidence::{
    ContinuityGates, measure_blend_continuity_evidence,
};
use popkern_blend::{
    ApexConstructionDiagnostic, BlendError, BlendResult, ContinuityLevel, CornerMethod,
    FilletOptions, FilletRadius, FilletSpec, SetbackFitBorder, SetbackFitMode, SetbackFitOutcome,
    fillet_edges,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

const VARIABLE: [f64; 8] = [0.02, 0.04, 0.06, 0.04, 0.02, 0.05, 0.03, 0.06];
const WIDE: [f64; 8] = [0.05, 0.15, 0.08, 0.2, 0.06, 0.12, 0.1, 0.18];

fn apex_fillet(n: usize, h: f64, radii: &[f64]) -> (BRep, Result<BlendResult, BlendError>) {
    let tol = Tolerance::default();
    let body = BRep::n_gonal_pyramid(n, 1.0, h, tol).expect("pyramid");
    let apex = body
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(v, _)| v)
        .expect("apex");
    let mut edges: Vec<_> = body
        .half_edges
        .values()
        .filter(|he| he.origin == apex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii)
        .map(|(&edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    let result = fillet_edges(&body, &specs, &options, tol, None);
    (body, result)
}

/// The served row's certified whole-blend level.
fn served(n: usize, h: f64, radii: &[f64]) -> Option<ContinuityLevel> {
    let label = format!("pyr{n} h{h} {radii:?}");
    let (body, result) = apex_fillet(n, h, radii);
    let result = result.unwrap_or_else(|e| panic!("{label}: {e}"));
    assert!(
        !result.corner_patches.is_empty() && result.corner_patches.iter().all(|p| p.setback_cap),
        "{label}: the setback rung serves"
    );
    assert!(result.brep.validate().is_ok(), "{label}");
    let pre: BTreeSet<FaceHandle> = body.faces.keys().collect();
    let blend: Vec<FaceHandle> = result
        .brep
        .faces
        .keys()
        .filter(|f| !pre.contains(f))
        .collect();
    measure_blend_continuity_evidence(
        &result.brep,
        &blend,
        ContinuityLevel::G2,
        &ContinuityGates::commercial(),
        33,
    )
    .certified_achieved
}

/// The served row's construction, the whole-blend evidence level (the
/// KP20 fans' seams sample their level but not all certify), and the
/// fan's recorded interior reading when it was past the fairness bound.
fn outcome(n: usize, h: f64, radii: &[f64]) -> (CornerMethod, ContinuityLevel, Option<f64>) {
    let (method, level, fairness, _) = outcome_and_fit(n, h, radii);
    (method, level, fairness)
}

/// [`outcome`] and the served fit attempt, from one build.
fn outcome_and_fit(
    n: usize,
    h: f64,
    radii: &[f64],
) -> (
    CornerMethod,
    ContinuityLevel,
    Option<f64>,
    Option<ServedFit>,
) {
    let label = format!("pyr{n} h{h} {radii:?}");
    let (body, result) = apex_fillet(n, h, radii);
    let result = result.unwrap_or_else(|e| panic!("{label}: {e}"));
    assert!(
        !result.corner_patches.is_empty() && result.corner_patches.iter().all(|p| p.setback_cap),
        "{label}: the setback rung serves"
    );
    assert!(result.brep.validate_all().is_empty(), "{label}");
    let pre: BTreeSet<FaceHandle> = body.faces.keys().collect();
    let blend: Vec<FaceHandle> = result
        .brep
        .faces
        .keys()
        .filter(|f| !pre.contains(f))
        .collect();
    let evidence = measure_blend_continuity_evidence(
        &result.brep,
        &blend,
        ContinuityLevel::G2,
        &ContinuityGates::commercial(),
        33,
    );
    let fairness = result
        .repair_report
        .apex_diagnostics
        .iter()
        .find_map(|d| match d {
            ApexConstructionDiagnostic::SetbackFit(record) => record.fairness.map(|f| f.kp20),
            _ => None,
        });
    (
        result.corner_patches[0].method,
        evidence.achieved,
        fairness,
        served_fit(&result),
    )
}

#[test]
fn octagonal_steep_apexes_serve_inside_the_fairness_bound() {
    // PK-613: the guide-faired KP20 fan serves the wide apex at G2
    // inside the fairness bound. The variable apex's fan reaches G1
    // only and reads 15.3 / r_min; on the corners border every fit
    // folds or pinches, and on the bridged border (PK-618) the tangent
    // fit serves at 2.40 / r_min, its sector beside the shortest strip
    // thin (regularity 1.8e-2) but fair.
    let (method, level, fairness) = outcome(8, 1.2, &WIDE);
    assert_eq!(
        (method, level, fairness),
        (CornerMethod::Gregory, ContinuityLevel::G2, None)
    );
    let (method, level, fairness, fit) = outcome_and_fit(8, 1.2, &VARIABLE);
    assert_eq!(
        (method, level),
        (CornerMethod::PlateMin, ContinuityLevel::G1)
    );
    assert!(fairness.is_some_and(|f| f > 4.0), "the fan's {fairness:?}");
    let fit = fit.expect("a served fit");
    assert_eq!(
        (fit.mode, fit.border),
        (SetbackFitMode::Tangent, SetbackFitBorder::Bridged)
    );
    assert!(fit.interior < 4.0, "{fit:?}");
    assert!(
        (5e-3..2e-2).contains(&fit.regularity),
        "a thin sector served inside the bound: {fit:?}"
    );
}

#[test]
fn octagonal_low_apexes_serve() {
    // Wide radii: the faired fan serves G2 inside the fairness bound
    // (PK-613; the tangent fit served these at G1 before). h 0.4
    // variable: the fan reads 5.5 / r_min and the corner-cell fit
    // (PK-595) serves, certified G2. h 0.8 variable: the fan reads 10.6
    // (a G1 downgrade); the tangent fit read 18.3 on the midpoint
    // layout and 2.76 with the spokes' feet nearest its centre, and
    // serves at G1 (PK-618).
    for h in [0.4, 0.8] {
        let (method, level, fairness) = outcome(8, h, &WIDE);
        assert_eq!(
            (method, level, fairness),
            (CornerMethod::Gregory, ContinuityLevel::G2, None),
            "h{h} wide"
        );
    }
    assert_eq!(served(8, 0.4, &VARIABLE), Some(ContinuityLevel::G2));
    let (method, level, fairness) = outcome(8, 0.8, &VARIABLE);
    assert_eq!(
        (method, level),
        (CornerMethod::PlateMin, ContinuityLevel::G1)
    );
    assert!(fairness.is_some_and(|f| f > 4.0), "{fairness:?}");
}

#[test]
fn steep_pentagon_and_hexagon_rows_serve_inside_the_fairness_bound() {
    // pyr5 h1.2 wide: the fan's interior reads 8.53 / r_min and the
    // corners border's fits 41 and 277 (a small corner between two
    // short sides gets a thin, bent sector); the bridged border's
    // curvature-weighted fit serves G2 at 3.89 (PK-618). pyr6 h1.2
    // variable: the unfaired fan pinched (4.5e5 / r_min); the faired
    // guide (PK-613) serves it at G2 inside the bound.
    let (method, level, fairness, fit) = outcome_and_fit(5, 1.2, &WIDE[..5]);
    assert_eq!(
        (method, level),
        (CornerMethod::PlateMin, ContinuityLevel::G2)
    );
    assert!(fairness.is_some_and(|f| f > 4.0), "the fan's {fairness:?}");
    let fit = fit.expect("a served fit");
    assert_eq!(
        (fit.mode, fit.border),
        (SetbackFitMode::CurvatureWeighted, SetbackFitBorder::Bridged)
    );
    assert!(fit.interior < 4.0 && fit.regularity >= 2e-2, "{fit:?}");
    let (method, level, fairness) = outcome(6, 1.2, &VARIABLE[..6]);
    assert_eq!(
        (method, level, fairness),
        (CornerMethod::Gregory, ContinuityLevel::G2, None)
    );
}

/// The served fit attempt of a row: mode, border, interior reading and
/// the certificate's regularity.
#[derive(Debug)]
struct ServedFit {
    mode: SetbackFitMode,
    border: SetbackFitBorder,
    interior: f64,
    regularity: f64,
}

fn served_fit(result: &BlendResult) -> Option<ServedFit> {
    result
        .repair_report
        .apex_diagnostics
        .iter()
        .find_map(|d| match d {
            ApexConstructionDiagnostic::SetbackFit(record) => {
                record.attempts.iter().find_map(|a| match &a.outcome {
                    SetbackFitOutcome::Served(c) => Some(ServedFit {
                        mode: a.mode,
                        border: a.border,
                        interior: a.interior.unwrap_or(f64::INFINITY),
                        regularity: c.min_regularity,
                    }),
                    _ => None,
                })
            }
            _ => None,
        })
}

/// A class A row serves the corner-cell fit (PK-595) on the bridged
/// border at G2, and records why the corners border and
/// the KP20 fan did not serve.
fn assert_bridged_fit_serve(n: usize, h: f64, radii: &[f64]) {
    let label = format!("pyr{n} h{h} {radii:?}");
    // The whole-blend certificate decides G2 on the pentagons; on the
    // hexagon some seams stay undecided at its budget (the fit's own
    // certificate reached G2 on every seam).
    let certified = served(n, h, radii);
    assert!(
        certified == Some(ContinuityLevel::G2) || (n == 6 && certified.is_none()),
        "{label}: certified {certified:?}"
    );
    let (_, result) = apex_fillet(n, h, radii);
    let result = result.unwrap_or_else(|e| panic!("{label}: {e}"));
    assert!(result.brep.validate_all().is_empty(), "{label}");
    assert!(
        result.corner_patches.len() == 2 * n
            && result
                .corner_patches
                .iter()
                .all(|p| p.method == CornerMethod::PlateMin),
        "{label}: one fitted sector per side of the 2n-sided border: {:?}",
        result.corner_patches
    );
    let record = result
        .repair_report
        .apex_diagnostics
        .iter()
        .find_map(|d| match d {
            ApexConstructionDiagnostic::SetbackFit(record) => Some(record),
            _ => None,
        })
        .unwrap_or_else(|| panic!("{label}: the fit's record"));
    assert!(
        record.kp20_decline.contains("arc/boundary crossings")
            && record.kp20_decline.contains("bridged border"),
        "{label}: {}",
        record.kp20_decline
    );
}

#[test]
fn a_narrow_strip_overrun_serves_the_bridged_border() {
    assert_bridged_fit_serve(5, 0.4, &VARIABLE[..5]);
    assert_bridged_fit_serve(5, 0.4, &WIDE[..5]);
}

#[test]
fn wide_overruns_serve_the_bridged_border() {
    assert_bridged_fit_serve(5, 0.8, &WIDE[..5]);
    assert_bridged_fit_serve(6, 0.4, &WIDE[..6]);
}
