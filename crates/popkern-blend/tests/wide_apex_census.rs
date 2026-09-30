//! PK-594: the wide-apex census rows that used to fail one strip with
//! `PartialFailure` and never reach a corner cap.
//!
//! Measured 2026-09-30 (G2 request, whole-blend evidence):
//! - Served since PK-591/593/596/597: pyr5 h1.2 wide and pyr6 h1.2
//!   variable (setback fan, G2); pyr8 h0.4/h0.8 variable and wide
//!   (tangent fit, certified G1).
//! - Class B, served by PK-594: pyr8 h1.2 variable and wide. The KP20
//!   fan folds and the curvature-weighted fit folds; the tangent fit's
//!   caps sat at served regularity 8.4e-3 and 1.6e-2, under a 2e-2
//!   floor the E6 data never supported (near-folds sat below 5e-3).
//!   They certify G1 on every seam (g1 ~1e-6) and pass `validate_all`.
//! - Class A, declined by name: pyr5 h0.4 variable and wide, pyr5 h0.8
//!   wide, pyr6 h0.4 wide. A thin strip's setback arc, run between the
//!   corners its wide neighbours fix, crosses its own trimmed boundary
//!   mid-arc; no setback carve exists on the corners-at-crossings
//!   border, and the fallbacks miss their supports. The convex bridge
//!   border that would serve them is PK-611.

use std::collections::BTreeSet;

use popkern_blend::blend_continuity_evidence::{
    ContinuityGates, measure_blend_continuity_evidence,
};
use popkern_blend::{
    BlendError, BlendFailure, BlendResult, ContinuityLevel, FilletOptions, FilletRadius,
    FilletSpec, fillet_edges,
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
fn octagonal_steep_apexes_serve_a_certified_g1_fit() {
    for radii in [VARIABLE, WIDE] {
        assert_eq!(
            served(8, 1.2, &radii),
            Some(ContinuityLevel::G1),
            "{radii:?}"
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
fn steep_pentagon_and_hexagon_rows_serve_the_setback_fan() {
    for (n, radii) in [(5, &WIDE[..5]), (6, &VARIABLE[..6])] {
        let label = format!("pyr{n} h1.2 {radii:?}");
        let (_, result) = apex_fillet(n, 1.2, radii);
        let result = result.unwrap_or_else(|e| panic!("{label}: {e}"));
        assert!(
            result
                .corner_patches
                .iter()
                .all(|p| p.setback_cap && p.achieved_continuity == ContinuityLevel::G2),
            "{label}: {:?}",
            result.corner_patches
        );
    }
}

#[test]
fn a_narrow_strip_overrun_declines_by_name() {
    for (n, h, radii) in [
        (5, 0.4, &VARIABLE[..5]),
        (5, 0.4, &WIDE[..5]),
        (5, 0.8, &WIDE[..5]),
        (6, 0.4, &WIDE[..6]),
    ] {
        let label = format!("pyr{n} h{h} {radii:?}");
        let (_, result) = apex_fillet(n, h, radii);
        let Err(BlendError::PartialFailure { failures, .. }) = result else {
            panic!("{label}: the corner declines typed");
        };
        assert!(
            failures.iter().any(|(_, failure)| matches!(
                failure,
                BlendFailure::VertexBlendFailed { reason, .. }
                    if reason.contains("setback arcs overrun a narrow strip")
            )),
            "{label}: the decline names the class: {failures:?}"
        );
    }
}
