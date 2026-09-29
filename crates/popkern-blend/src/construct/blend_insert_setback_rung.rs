//! The setback vertex-blend rung of `insert_vertex_blend` (PK-094
//! surgery entry, PK-225 N=3 production gate, PK-574 fallback ladder).
//! See the call site for the rung ordering contract and
//! `docs/architecture/variational-corner-fit.md` for the ladder.

use super::InsertVertexBlendOutput;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};

use crate::construct::setback_decline::{SetbackCapStrategy, SetbackDecline, SetbackStage};
use crate::construct::setback_vertex_topology::{
    SETBACK_TOPOLOGY_LAST, SetbackTopologyOutcome, try_setback_boundary_surgery,
};
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{
    ApexConstructionDiagnostic, ContinuityLevel, CornerMethod, CornerPatchInfo,
    CornerSubPatchFamily, CornerSubPatchProvenance, CornerSubPatchRole, SetbackFitAttempt,
    SetbackFitDiagnostic, SetbackFitMode, SetbackFitOutcome,
};
use crate::vertex_classifier::VertexBlendClassification;

/// The rung's result: `served` is `Some` when the rung decided (served
/// or hard-errored), `None` to fall through to the vertex-blend corner
/// machinery; `diagnostics` carry the variational-fit record whenever
/// the fit ran, served or not.
pub(super) type SetbackRungOutcome = (
    Option<Result<InsertVertexBlendOutput, String>>,
    Vec<ApexConstructionDiagnostic>,
);

/// Attempt the setback rung: the KP20 cap, its G1 downgrade on a
/// seam-G2-only failure, then the variational fit (curvature-weighted
/// under a G2 request, then tangent), each transactional on the same
/// carve.
#[allow(clippy::too_many_arguments)]
pub(super) fn try_setback_rung(
    brep: &mut BRep,
    vertex: VertexHandle,
    edge_indices: &[usize],
    fillet_surfaces_pre: &[FilletSurfaceResult],
    tolerance: &Tolerance,
    requested_continuity: ContinuityLevel,
    classification: VertexBlendClassification,
    setback_min_n: usize,
    cap_presolve: Option<&crate::construct::setback_cap_presolve::SetbackCapPresolveEntry>,
) -> SetbackRungOutcome {
    if !(crate::construct::setback_vertex_topology::setback_topology_commit_enabled()
        && edge_indices.len() >= setback_min_n)
    {
        return (None, Vec::new());
    }
    // On any surgery error the snapshot restore inside the surgery has
    // already run; the next attempt starts from the same model.
    let commit_snapshot = brep.clone();
    let gate_seam_g2 = requested_continuity == ContinuityLevel::G2;
    let surgery = |brep: &mut BRep, strategy| -> Result<SetbackTopologyOutcome, SetbackDecline> {
        if crate::dispatch_overrides::setback_decline_forced() {
            return Err("forced decline (dispatch override)".into());
        }
        try_setback_boundary_surgery(
            brep,
            &commit_snapshot,
            vertex,
            fillet_surfaces_pre,
            tolerance,
            strategy,
            cap_presolve,
        )
    };
    let info = |face, method, achieved, provenance| CornerPatchInfo {
        vertex,
        n_sides: edge_indices.len(),
        method,
        achieved_continuity: achieved,
        classification,
        face,
        bridge_fallbacks: Vec::new(),
        overlap_remediated: false,
        relaxed_continuity_coplanarity: None,
        cap_sphere: None,
        sub_patch_provenance: provenance,
        setback_cap: true,
        exact_boundary_cap: false,
    };
    // 1. The KP20 cap. It serves both continuity levels: it is the
    //    proven watertight tiling and a G2 cap satisfies G1.
    let kp20 = surgery(brep, SetbackCapStrategy::Kp20 { gate_seam_g2 });
    let mut kp20_reason = match kp20 {
        Ok(outcome) => {
            // A KP20 (G2-interior) cap on G1 sweeps yields a G1 blend.
            let achieved = match requested_continuity {
                ContinuityLevel::G2 => ContinuityLevel::G2,
                _ => ContinuityLevel::G1,
            };
            return (
                Some(commit(brep, &commit_snapshot, &outcome, |f, _| {
                    info(f, CornerMethod::Gregory, achieved, None)
                })),
                Vec::new(),
            );
        }
        Err(d) if d.stage == SetbackStage::Carve => {
            record_decline(vertex, &d.reason);
            return (None, Vec::new());
        }
        Err(d) => d.reason,
    };
    // 2. A fan that fails ONLY the seam-G2 postcondition still beats
    //    every other fallback measured on that class (rim+verticals box:
    //    the declined fan's seams at g1 2.9e-3 / g2 5.4e-2 against the
    //    Gregory fallback's 27° crease); serve it as a measured G1
    //    downgrade.
    if gate_seam_g2 && kp20_reason.starts_with("fan seam G2 postcondition failed") {
        match surgery(
            brep,
            SetbackCapStrategy::Kp20 {
                gate_seam_g2: false,
            },
        ) {
            Ok(outcome) => {
                return (
                    Some(commit(brep, &commit_snapshot, &outcome, |f, _| {
                        info(f, CornerMethod::Gregory, ContinuityLevel::G1, None)
                    })),
                    Vec::new(),
                );
            }
            Err(d) => kp20_reason = format!("{kp20_reason}; G1 downgrade: {}", d.reason),
        }
    }
    // 3. The variational fit on the same carve (PK-574).
    let mut record = SetbackFitDiagnostic {
        vertex,
        kp20_decline: kp20_reason,
        attempts: Vec::new(),
        served: None,
    };
    let modes: &[SetbackFitMode] = if gate_seam_g2 {
        &[SetbackFitMode::CurvatureWeighted, SetbackFitMode::Tangent]
    } else {
        &[SetbackFitMode::Tangent]
    };
    for &mode in modes {
        let attempt = surgery(brep, SetbackCapStrategy::VariationalFit(mode));
        match attempt {
            Ok(outcome) => {
                let Some(serve) = outcome.fit.clone() else {
                    *brep = commit_snapshot.clone();
                    record.attempts.push(SetbackFitAttempt {
                        mode,
                        outcome: SetbackFitOutcome::Declined(
                            "fit surgery returned no fit record".into(),
                        ),
                    });
                    continue;
                };
                if let Err(e) = brep.validate() {
                    *brep = commit_snapshot.clone();
                    record.attempts.push(SetbackFitAttempt {
                        mode,
                        outcome: SetbackFitOutcome::Rejected {
                            certified: serve.certified,
                            reason: format!("validation: {e}"),
                        },
                    });
                    continue;
                }
                let achieved = match mode {
                    SetbackFitMode::CurvatureWeighted => ContinuityLevel::G2,
                    SetbackFitMode::Tangent => ContinuityLevel::G1,
                };
                record.attempts.push(SetbackFitAttempt {
                    mode,
                    outcome: SetbackFitOutcome::Served(serve.certified),
                });
                record.served = Some(achieved);
                let n = outcome.cap_faces.len();
                let infos = outcome
                    .cap_faces
                    .iter()
                    .enumerate()
                    .map(|(s, &f)| {
                        info(
                            f,
                            CornerMethod::PlateMin,
                            achieved,
                            Some(CornerSubPatchProvenance {
                                family: CornerSubPatchFamily::VariationalFit,
                                sector_index: s,
                                sector_count: n,
                                piece_index: 0,
                                piece_count: 1,
                                role: CornerSubPatchRole::Sector,
                            }),
                        )
                    })
                    .collect();
                let diagnostics = vec![ApexConstructionDiagnostic::SetbackFit(record)];
                return (Some(Ok((infos, None, diagnostics))), Vec::new());
            }
            Err(d) => {
                if popkern_core::trace::trace_enabled("setback-fit") {
                    eprintln!(
                        "  [setback-fit] {mode:?} declined at {vertex:?}: {}",
                        d.reason
                    );
                }
                record.attempts.push(SetbackFitAttempt {
                    mode,
                    outcome: match d.fit_certified {
                        Some(certified) => SetbackFitOutcome::Rejected {
                            certified,
                            reason: d.reason,
                        },
                        None => SetbackFitOutcome::Declined(d.reason),
                    },
                });
            }
        }
    }
    record_decline(vertex, &record.kp20_decline);
    (None, vec![ApexConstructionDiagnostic::SetbackFit(record)])
}

/// Validate a served KP20 surgery and build its corner records; a
/// validation failure restores the snapshot and is a hard error.
fn commit(
    brep: &mut BRep,
    snapshot: &BRep,
    outcome: &SetbackTopologyOutcome,
    info: impl Fn(popkern_core::topology::FaceHandle, usize) -> CornerPatchInfo,
) -> Result<InsertVertexBlendOutput, String> {
    if let Err(e) = brep.validate() {
        *brep = snapshot.clone();
        *SETBACK_TOPOLOGY_LAST
            .lock()
            .unwrap_or_else(|p| p.into_inner()) =
            Some(format!("setback commit failed validation: {e}"));
        return Err(format!("setback vertex blend failed validation: {e}"));
    }
    let infos = outcome
        .cap_faces
        .iter()
        .enumerate()
        .map(|(s, &face)| info(face, s))
        .collect();
    Ok((infos, None, Vec::new()))
}

fn record_decline(vertex: VertexHandle, reason: &str) {
    *SETBACK_TOPOLOGY_LAST
        .lock()
        .unwrap_or_else(|p| p.into_inner()) = Some(format!("setback commit declined: {reason}"));
    if popkern_core::trace::trace_enabled("corner-carve") {
        eprintln!("  [setback] declined at {vertex:?}: {reason}");
    }
}
