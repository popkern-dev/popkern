//! Certified claim for a variational-fit cap (PK-574): every seam that
//! touches a cap face is bounded by the PK-499 certifier, and the claim
//! is the highest level whose commercial gates every bound meets. A seam
//! whose bounds do not cover it counts as failing: no level is claimed
//! from samples.

use popkern_core::geometry::certify::SeamCertifyOptions;
use popkern_core::topology::{BRep, FaceHandle};

use crate::blend_continuity_certification::evidence_certify_options;
use crate::blend_continuity_evidence::ContinuityGates;
use crate::construct::setback_decline::CapCommitError;
use crate::types::{ContinuityLevel, SetbackFitCertified};

/// Normal-sine resolution of the recorded bounds: three decades inside
/// the commercial G1 gate.
const RECORD_SINE: f64 = 1e-5;

/// Intervals per seam. The corner-cell fit's degree-9 spokes are G2 by
/// construction (every certified lower bound 0) but their curvature
/// enclosures widen beside the valence-n centre: on the cube corner
/// 0.05/0.2/0.3 one spoke stopped at the evidence's 256 intervals with
/// an upper bound 7.2e-2 over the 5e-2 gate (PK-595).
const FIT_CERTIFY_BUDGET: usize = 1024;

/// Certify the cap's seams and derive its level.
pub(super) fn certify_cap(
    brep: &BRep,
    cap_faces: &[FaceHandle],
    gates: &ContinuityGates,
    report: &popkern_core::geometry::variational_cap::VariationalCapReport,
) -> Result<SetbackFitCertified, CapCommitError> {
    let caps: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();
    // The claim needs bounds decided against the commercial gates; the
    // record also resolves the normal sine to `RECORD_SINE`. Pure
    // decision mode left recorded bounds just under the G1 gate
    // (pyramid3 h3: rim g1 bound 2.8e-2 against 9.7e-7 refined); full
    // refinement ran spoke curvature enclosures to the interval budget
    // (certification 1.2 s mean, 3.8 s max per corner, E4). Deciding
    // g0 and g2 at the gates and g1 at `RECORD_SINE` keeps both.
    let options = SeamCertifyOptions {
        targets: popkern_core::topology::edge_seam_certificate::seam_certify_options(brep).targets,
        decision_gates: Some([gates.g0_gap, RECORD_SINE, gates.g2_curvature]),
        budget: FIT_CERTIFY_BUDGET,
        ..evidence_certify_options(brep, gates)
    };
    // Every seam touching the cap: (spoke?, edge). Each certificate is
    // a pure read of the BRep, so the seams certify in parallel (PK-598:
    // octagon8 mixed, 16 seams, 1.8 s serial); the fold is a max, so
    // the claim does not depend on order.
    let targets: Vec<(bool, popkern_core::topology::EdgeHandle)> = brep
        .edges
        .iter()
        .filter_map(|(edge, data)| {
            let he = brep.half_edges.get(data.half_edge)?;
            let tw = brep.half_edges.get(he.twin)?;
            let (fa, fb) = (brep.loops[he.loop_].face, brep.loops[tw.loop_].face);
            let (ca, cb) = (caps.contains(&fa), caps.contains(&fb));
            ((ca || cb) && fa != fb).then_some((ca && cb, edge))
        })
        .collect();
    let certify = |&(is_spoke, edge): &(bool, popkern_core::topology::EdgeHandle)| {
        let cert =
            popkern_core::topology::edge_seam_certificate::certify_edge_seam(brep, edge, &options);
        let bounds = if cert.bounds_cover_seam() {
            [cert.g0_upper, cert.g1_upper, cert.g2_upper]
        } else {
            [f64::INFINITY; 3]
        };
        (is_spoke, bounds)
    };
    #[cfg(feature = "parallel")]
    let certified: Vec<(bool, [f64; 3])> = {
        use rayon::prelude::*;
        targets.par_iter().map(certify).collect()
    };
    #[cfg(not(feature = "parallel"))]
    let certified: Vec<(bool, [f64; 3])> = targets.iter().map(certify).collect();
    let (mut rim, mut spoke) = ([0.0f64; 3], [0.0f64; 3]);
    for (is_spoke, bounds) in &certified {
        let slot = if *is_spoke { &mut spoke } else { &mut rim };
        for (s, b) in slot.iter_mut().zip(bounds) {
            *s = s.max(*b);
        }
    }
    let seams = certified.len();
    if seams == 0 {
        return Err("variational fit cap has no seams to certify".into());
    }
    let worst = [
        rim[0].max(spoke[0]),
        rim[1].max(spoke[1]),
        rim[2].max(spoke[2]),
    ];
    let level = if worst[0] <= gates.g0_gap && worst[1] <= gates.g1_normal {
        if worst[2] <= gates.g2_curvature {
            ContinuityLevel::G2
        } else {
            ContinuityLevel::G1
        }
    } else {
        ContinuityLevel::G0
    };
    Ok(SetbackFitCertified {
        level,
        rim,
        spoke,
        refinement_levels: report.refinement_levels,
        unknowns: report.unknowns,
        min_regularity: report.min_regularity,
    })
}
