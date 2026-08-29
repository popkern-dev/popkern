//! Low-valence PnS cap zebra evidence, split out of
//! `pns_low_valence_gates.rs` for the 500-line structural gate.
//!
//! The PK-192 campaign's ignored diagnosis probes that lived here
//! (loop anatomy, pcurve attribution, pipeline-end fidelity,
//! rim-vs-spoke deviation) were retired with the route decision;
//! their findings are recorded on the ticket and in the live-gate
//! doc comments. What remains measures through
//! `run_pyramid_apex_outcome`, the same harness the gates use.

mod kp_gate_helpers;
#[path = "support/pns_apex_outcome.rs"]
mod pns_apex_outcome;

use pns_apex_outcome::run_pyramid_apex_outcome;

/// PK-197 criterion: zebra evidence across the internal cap seams of
/// the n=3 and n=4 mixed apexes — the stripe-phase form of the zebra
/// render, measured where a render would be eyeballed. A creased seam
/// reads a large stripe-phase jump; a G1 seam reads near zero.
///
/// Asserted (not just printed): with `harmonize_mismatched_spoke_nets`
/// plus `cross_spoke_g1` the caps must stay below a phase jump of 0.25
/// stripes on every internal seam. Pre-fix, the n=3 crease (89.9 deg)
/// saturated the measure.
#[test]
fn zebra_refutes_the_low_valence_cap_crease() {
    use popkern_core::geometry::Vec3;
    for (n, radii) in [
        (3_usize, vec![0.04, 0.06, 0.05]),
        (4, vec![0.04, 0.06, 0.05, 0.055]),
    ] {
        // n = 3 measures the CUBE corner — the fixture the
        // exact-boundary cross-spoke + harmonization machinery
        // serves. The old n=3 row measured a pyramid BASE corner
        // (the valence-3 pick was ambiguous on a pyramid3) whose
        // defective direct serve now declines (PK-260); the true
        // pyramid apex serves the rational-Gregory fan, which is
        // creased by construction and never carried the machinery
        // this probe refutes the crease FOR.
        let o = if n == 3 {
            pns_apex_outcome::run_cube_corner_outcome(&radii)
        } else {
            run_pyramid_apex_outcome(n, &radii)
        };
        assert!(o.metrics.seam_count > 0, "n={n}: no internal seams");
        let brep = &o.brep;
        let caps: std::collections::BTreeSet<_> = o.cap_face_handles.iter().copied().collect();
        let mut measured = 0usize;
        let mut worst = 0.0_f64;
        for (_eh, edge) in brep.edges.iter() {
            let Some(he) = brep.half_edges.get(edge.half_edge) else {
                continue;
            };
            let Some(twin) = brep.half_edges.get(he.twin) else {
                continue;
            };
            let (Some(la), Some(lb)) = (brep.loops.get(he.loop_), brep.loops.get(twin.loop_))
            else {
                continue;
            };
            let (fa, fb) = (la.face, lb.face);
            if fa == fb || !caps.contains(&fa) || !caps.contains(&fb) {
                continue;
            }
            if let Some(kink) = popkern_blend::zebra::zebra_boundary_kink(
                brep,
                fa,
                fb,
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(1.0, 0.0, 0.0),
                0.02,
                24,
                1e-6,
            ) {
                measured += 1;
                worst = worst.max(kink);
            }
        }
        println!(
            "[pk197-zebra] n={n} seams_measured={measured} worst_stripe_phase_jump={worst:.4}"
        );
        assert!(measured > 0, "n={n}: zebra measured no internal seams");
        assert!(
            worst < 0.25,
            "n={n}: internal cap seam shows a zebra stripe-phase jump of {worst:.4} \
             (creased); the cross-spoke + harmonization fix should have refuted this"
        );
    }
}
