//! PK-594: the wide-apex census rows that used to fail one strip with
//! `PartialFailure` and never reach a corner cap.
//!
//! Measured 2026-09-30 (G2 request, whole-blend evidence):
//! - Served since PK-591/593/596/597: pyr5 h1.2 wide (setback fan, G2)
//!   and pyr8 h0.4/h0.8 variable and wide (tangent fit, certified G1).
//!   pyr6 h1.2 variable's fan pinches (interior 4.5e5 / r_min at
//!   regularity 2.4e-3); since PK-612 the fairer tangent fit (2.09)
//!   serves it at certified G1.
//! - Class B, pyr8 h1.2 variable and wide: declined. The KP20 fan
//!   folds and the curvature-weighted fit folds. PK-594 served the
//!   tangent fit's caps at served regularity 8.4e-3 and 1.6e-2, but
//!   PK-612 measured their interiors at 2572 and 2850 times 1/r_min:
//!   near-fold pinches the seam certificates cannot see. With the fold
//!   floor back at 2e-2 the corner declines typed.
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
    ApexConstructionDiagnostic, BlendError, BlendFailure, BlendResult, ContinuityLevel,
    CornerMethod, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
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

#[test]
fn octagonal_steep_apexes_decline_where_every_cap_nearly_folds() {
    for radii in [VARIABLE, WIDE] {
        let (_, result) = apex_fillet(8, 1.2, &radii);
        let Err(BlendError::PartialFailure { failures, .. }) = result else {
            panic!("{radii:?}: the corner declines typed");
        };
        assert!(
            failures
                .iter()
                .all(|(_, f)| matches!(f, BlendFailure::VertexBlendFailed { .. })),
            "{radii:?}: {failures:?}"
        );
    }
}

#[test]
fn octagonal_low_apexes_serve_a_certified_g1_fit() {
    for h in [0.4, 0.8] {
        for radii in [VARIABLE, WIDE] {
            assert_eq!(
                served(8, h, &radii),
                Some(ContinuityLevel::G1),
                "h{h} {radii:?}"
            );
        }
    }
}

#[test]
fn steep_pentagon_keeps_its_fan_and_the_hexagon_serves_the_fairer_fit() {
    // pyr5 h1.2 wide: the fan's interior reads 8.53 / r_min, past the
    // fairness bound, but its fit nearly folds (declined at the 2e-2
    // floor), so the fan serves G2 with both readings recorded.
    let label = "pyr5 h1.2 wide";
    let (_, result) = apex_fillet(5, 1.2, &WIDE[..5]);
    let result = result.unwrap_or_else(|e| panic!("{label}: {e}"));
    assert!(
        result
            .corner_patches
            .iter()
            .all(|p| p.setback_cap && p.achieved_continuity == ContinuityLevel::G2),
        "{label}: {:?}",
        result.corner_patches
    );
    let fairness = result
        .repair_report
        .apex_diagnostics
        .iter()
        .find_map(|d| match d {
            ApexConstructionDiagnostic::SetbackFit(record) => record.fairness,
            _ => None,
        })
        .unwrap_or_else(|| panic!("{label}: the fairness record"));
    assert!(
        fairness.kp20 > fairness.bound && fairness.fit.is_none(),
        "{label}: {fairness:?}"
    );
    // pyr6 h1.2 variable: the pinched fan yields to the tangent fit.
    assert_eq!(
        served(6, 1.2, &VARIABLE[..6]),
        Some(ContinuityLevel::G1),
        "pyr6 h1.2 variable"
    );
}

/// A class A row serves the tangent fit on the bridged border,
/// certified G1 on every seam, and records why the corners border and
/// the KP20 fan did not serve.
fn assert_bridged_fit_serve(n: usize, h: f64, radii: &[f64]) {
    let label = format!("pyr{n} h{h} {radii:?}");
    assert_eq!(
        served(n, h, radii),
        Some(ContinuityLevel::G1),
        "{label}: certified"
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
