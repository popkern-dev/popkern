//! PK-212 — rim (cap ↔ strip) tangency gate for the setback fan,
//! split from `kp_polar_commercial_gates.rs` for the 500-line
//! structural gate. Same fixtures, same instrument.

#[path = "support/kp_polar_fixtures.rs"]
mod fixtures;
mod kp_gate_helpers;

use fixtures::fillet_n_gonal_pyramid_result_with_continuity;

/// PK-212 — rim (cap ↔ strip) tangency gate. The setback arcs' cross
/// tangency was NEVER measured before the PK-188 evidence substrate
/// (every commercial gate is cap↔cap): the border derivation's
/// point-reflection prolongation tilted the cap's cross-jets out of
/// the strip tangent planes by the chord sagitta κℓ/2 — a UNIFORM
/// 2.6°..6.9° along every arc. Fixed by the tangent-plane
/// prolongation in the border extraction plus the committed-fan rim
/// cross-jet reconciliation (shared-corner clamped-linear deltas,
/// curvature-invariant depth application). Gated at the commercial
/// 5e-2 bar on every fixture: the blunt pentagon G2 request's 6.7e-2
/// mid-arc bump (once pinned at 8e-2 after every k-varying interior
/// correction was falsified) went with PK-398's field-sampled cross
/// rows (measured 2.9e-2; its corners are gapped — the G2 sweeps'
/// contact curves miss by 3.4e-3 — and serve with the χ̃-calibrated
/// containment ends).
#[test]
fn rim_tangency_meets_commercial_bar() {
    for (n, h, continuity, fixture, bar) in [
        (
            5usize,
            3.0,
            popkern_blend::ContinuityLevel::G1,
            "pentagon5_blunt/G1",
            5.0e-2,
        ),
        (
            5,
            3.0,
            popkern_blend::ContinuityLevel::G2,
            "pentagon5_blunt/G2",
            5.0e-2,
        ),
        (
            6,
            3.0,
            popkern_blend::ContinuityLevel::G2,
            "hexagon6_blunt/G2",
            5.0e-2,
        ),
        (
            5,
            0.8,
            popkern_blend::ContinuityLevel::G2,
            "pentagon5_sharp/G2",
            5.0e-2,
        ),
        (
            5,
            0.8,
            popkern_blend::ContinuityLevel::G1,
            "pentagon5_sharp/G1",
            5.0e-2,
        ),
    ] {
        let radii = vec![0.05; n];
        let result = fillet_n_gonal_pyramid_result_with_continuity(n, 1.0, h, &radii, continuity);
        let caps: Vec<_> = result
            .corner_patches
            .iter()
            .filter(|p| p.n_sides == n)
            .map(|p| p.face)
            .collect();
        assert!(!caps.is_empty(), "[{fixture}] no cap faces");
        let ev = popkern_blend::measure_blend_continuity_evidence(
            &result.brep,
            &caps,
            continuity,
            &popkern_blend::ContinuityGates::commercial(),
            33,
        );
        let (mut count, mut worst) = (0usize, 0.0_f64);
        for s in ev
            .seams
            .iter()
            .filter(|s| matches!(s.kind, popkern_blend::BlendSeamKind::BlendToSupport) && s.gated)
        {
            count += 1;
            worst = worst.max(s.worst_g1_normal);
        }
        println!("[{fixture}] rim seams={count} worst_g1={worst:.3e} (bar {bar:.1e})");
        assert!(count > 0, "[{fixture}] rim seams must measure");
        assert!(
            worst < bar,
            "[{fixture}] rim tangency {worst:.3e} exceeds the {bar:.1e} bar"
        );
    }
}
