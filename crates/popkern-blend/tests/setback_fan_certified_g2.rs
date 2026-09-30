//! PK-596: the KP20 setback fan's G2 claim agrees with the certified
//! evidence.
//!
//! Before PK-596 the fan's seam-G2 postcondition read its own sampling
//! (end margin, finite-difference normals) and passed on pyramid-4
//! apexes at r = 0.05 whose spokes certifiably jump 8.7e-2 near the
//! centre: every corner piece reported G2 while the whole-blend
//! evidence (PK-499 certificates) read G1. The rung now reads the
//! fan's verdict from the same evidence, hands a failing fan to the
//! curvature-weighted fit, and serves the fan as a recorded G1
//! downgrade when the fit declines. Each row pins the corner record,
//! the served construction and the whole-blend verdict, and that the
//! record never claims more than the evidence.

use std::collections::BTreeSet;

use popkern_blend::blend_continuity_evidence::{
    ContinuityGates, measure_blend_continuity_evidence,
};
use popkern_blend::dispatch_overrides::with_cap_fairness_ladder_disabled;
use popkern_blend::{
    ApexConstructionDiagnostic, ContinuityLevel, CornerMethod, FilletOptions, FilletRadius,
    FilletSpec, fillet_edges,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// `(corner record, served method, whole-blend evidence)` for the
/// pyramid-4 apex of height `h` at r = 0.05 under a G2 request.
fn serve(h: f64) -> (ContinuityLevel, CornerMethod, ContinuityLevel) {
    let tol = Tolerance::default();
    let body = BRep::n_gonal_pyramid(4, 1.0, h, tol).expect("pyramid");
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
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.05),
        })
        .collect();
    let pre: BTreeSet<FaceHandle> = body.faces.keys().collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    let result = fillet_edges(&body, &specs, &options, tol, None).expect("fillet");
    assert!(result.brep.validate().is_ok(), "h{h}");
    let fairness = result
        .repair_report
        .apex_diagnostics
        .iter()
        .find_map(|d| match d {
            ApexConstructionDiagnostic::SetbackFit(record) => record.fairness,
            _ => None,
        });
    if let Some(f) = fairness {
        // Whichever cap served, the record names both readings.
        assert!(f.kp20 > f.bound, "h{h}: {f:?}");
    }
    let patches = &result.corner_patches;
    assert!(
        !patches.is_empty() && patches.iter().all(|p| p.setback_cap),
        "h{h}: the setback rung serves: {patches:?}"
    );
    let record = patches
        .iter()
        .map(|p| p.achieved_continuity)
        .min()
        .expect("patches");
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
    assert!(
        record <= evidence.achieved,
        "h{h}: corner record {record:?} claims more than the evidence {:?}",
        evidence.achieved
    );
    (record, patches[0].method, evidence.achieved)
}

#[test]
fn a_fan_certified_below_g2_serves_a_recorded_g1_downgrade() {
    // Measured 2026-09-30: h 3's spokes certify a g2 lower bound
    // 8.7e-2; the curvature-weighted fit declines and the tangent fit's
    // seam g2 (19) is far worse than the fan's. The fan's interior also
    // bends past the fairness bound (4.22 / r_min, PK-612), so the
    // downgrade path is pinned with the fairness ladder off.
    assert_eq!(
        with_cap_fairness_ladder_disabled(|| serve(3.0)),
        (
            ContinuityLevel::G1,
            CornerMethod::Gregory,
            ContinuityLevel::G1
        )
    );
}

#[test]
fn an_unfair_fan_below_g2_yields_to_the_fairer_tangent_fit() {
    // PK-612, measured 2026-09-30: h 3's fan reads interior curvature
    // 4.22 / r_min, the tangent fit 1.53 with cleaner rims (zebra
    // cap-rim normal p95 1.3e-2 against the fan's 4.1e-2). Both are G1.
    assert_eq!(
        serve(3.0),
        (
            ContinuityLevel::G1,
            CornerMethod::PlateMin,
            ContinuityLevel::G1
        )
    );
}

#[test]
fn the_faired_fan_certifies_g2_at_h_0_8() {
    // PK-613: before the guide fairing, h 0.8's fan spokes certified a
    // g2 lower bound 5.3e-2 (a G1 downgrade) and its interior read
    // 4.35 / r_min. The faired guide's fan certifies G2 on every seam.
    assert_eq!(
        serve(0.8),
        (
            ContinuityLevel::G2,
            CornerMethod::Gregory,
            ContinuityLevel::G2
        )
    );
}

#[test]
fn a_fan_below_g2_hands_the_corner_to_a_g2_fit() {
    // Measured 2026-09-30: the fan reads g2 2.9e-1; the
    // curvature-weighted fit certifies G2 (whole-blend g2 2.8e-2).
    assert_eq!(
        serve(0.25),
        (
            ContinuityLevel::G2,
            CornerMethod::PlateMin,
            ContinuityLevel::G2
        )
    );
}

#[test]
fn a_fan_certified_at_g2_keeps_its_g2_claim() {
    // Negative control: the gate does not demote a genuine G2 fan
    // (measured whole-blend g2 1.6e-2).
    assert_eq!(
        serve(1.2),
        (
            ContinuityLevel::G2,
            CornerMethod::Gregory,
            ContinuityLevel::G2
        )
    );
}
