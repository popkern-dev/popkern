//! Transactional PnS dispatch acceptance for vertex-blend insertion.

use popkern_core::error::ValidationError;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::blend::BlendBoundary;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo};
use crate::vertex_classifier::VertexBlendClassification;

use super::blend_insert_diagnostics::pns_dispatch_declined_diagnostic;

pub(super) struct PnsDispatchAcceptanceInput<'a> {
    pub brep: &'a mut BRep,
    pub blend_face: FaceHandle,
    pub boundary: &'a BlendBoundary,
    pub fillet_surfaces: &'a [FilletSurfaceResult],
    pub vertex: VertexHandle,
    pub classification: VertexBlendClassification,
    pub n_sides: usize,
    pub diagnostics: &'a mut Vec<ApexConstructionDiagnostic>,
}

pub(super) fn try_pns_dispatch_with_acceptance(
    input: PnsDispatchAcceptanceInput<'_>,
) -> Option<Vec<CornerPatchInfo>> {
    let PnsDispatchAcceptanceInput {
        brep,
        blend_face,
        boundary,
        fillet_surfaces,
        vertex,
        classification,
        n_sides,
        diagnostics,
    } = input;
    let pns_snapshot = brep.clone();
    let patches = match crate::blend::pns_promote::try_pns_dispatch(
        brep,
        blend_face,
        boundary,
        fillet_surfaces,
        vertex,
        classification,
        boundary.bridge_fallbacks.clone(),
    ) {
        Ok(patches) => patches,
        Err(error) => {
            // Restore before returning: a partially-inserted cap must
            // not survive a decline. Then record WHY — this path used
            // to return None silently, so every construction failure
            // looked identical to "PnS was never tried".
            *brep = pns_snapshot;
            diagnostics.push(pns_dispatch_declined_diagnostic(
                n_sides,
                "construction_failed",
                Some(error.to_string()),
            ));
            return None;
        }
    };

    let cap_faces = patches.iter().map(|patch| patch.face).collect::<Vec<_>>();
    // Measure pcurve fidelity on the cap faces and their neighbours —
    // the faces the acceptance gate judges — with the validated
    // instrument. Read-only; reported through the non-blocking
    // channel so it cannot change which route serves the corner.
    let mut audited: Vec<popkern_core::topology::FaceHandle> = cap_faces.clone();
    for &cap in &cap_faces {
        let Some(face_data) = brep.faces.get(cap) else {
            continue;
        };
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut he = start;
        for _ in 0..32 {
            let Some(hd) = brep.half_edges.get(he) else {
                break;
            };
            if let Some(twin) = brep.half_edges.get(hd.twin)
                && let Some(l) = brep.loops.get(twin.loop_)
                && !audited.contains(&l.face)
            {
                audited.push(l.face);
            }
            he = hd.next;
            if he == start {
                break;
            }
        }
    }
    // PK-192: the acceptance gate below judges the cap against its
    // neighbours, but PnS dispatch runs BEFORE the pipeline's
    // pcurve-completion passes — so those neighbours are unfinished at
    // this moment. Measured on the same class of face: 3.277 here
    // against 4.140e-7 at pipeline end, a factor of eight million. The
    // cap was being refused for its neighbours' unfinished state.
    //
    // Complete them first, so the gate measures geometry rather than
    // timing.
    // Verified: this succeeds on every audited face, including the one
    // the acceptance gate then reports against — so a silent recompute
    // failure is NOT the explanation for the two measurements
    // disagreeing (see the open inconsistency on PK-192).
    let mut recompute_failures = 0usize;
    for &face in &audited {
        if brep.recompute_face_pcurves(face).is_err() {
            recompute_failures += 1;
        }
    }
    let fidelity = crate::pcurve_fidelity::measure_pcurve_fidelity(brep, &audited);
    // A bi-4 patch has domain [0,1]^2 and EXTRAPOLATES violently
    // outside it, so a large fidelity reading on a cap face is
    // ambiguous until the domain is checked: a pcurve straying out of
    // domain is a real defect, an instrument sampling out of domain is
    // not. Report the UV extent so the two are distinguishable.
    let mut uv_min = (f64::MAX, f64::MAX);
    let mut uv_max = (f64::MIN, f64::MIN);
    for &cap in &cap_faces {
        let Some(fd) = brep.faces.get(cap) else {
            continue;
        };
        let Some(ld) = brep.loops.get(fd.outer_loop) else {
            continue;
        };
        let start = ld.half_edge;
        let mut he = start;
        for _ in 0..32 {
            let Some(hd) = brep.half_edges.get(he) else {
                break;
            };
            if let Some(pc) = hd.pcurve.as_ref() {
                for k in 0..=10 {
                    let t = pc.t_start + f64::from(k) / 10.0 * (pc.t_end - pc.t_start);
                    let uv = pc.point_at(t);
                    uv_min = (uv_min.0.min(uv.x), uv_min.1.min(uv.y));
                    uv_max = (uv_max.0.max(uv.x), uv_max.1.max(uv.y));
                }
            }
            he = hd.next;
            if he == start {
                break;
            }
        }
    }
    diagnostics.push(pns_dispatch_declined_diagnostic(
        n_sides,
        "pcurve_fidelity_audit",
        Some(format!(
            "{fidelity:?};cap_uv_extent=({:.4},{:.4})..({:.4},{:.4});recompute_failures={}",
            uv_min.0, uv_min.1, uv_max.0, uv_max.1, recompute_failures
        )),
    ));
    crate::construct::corner_patch_support_completion::complete_external_cap_support_pcurves(
        brep,
        &patches,
        diagnostics,
    );
    let topology_error = brep.validate_collect_fast().into_iter().find(|error| {
        matches!(
            error,
            ValidationError::BrokenCycle { .. }
                | ValidationError::ShellGap { .. }
                | ValidationError::OpenLoop { .. }
        )
    });
    let continuity_error = crate::construct::cap_outer_acceptance::cap_support_acceptance_error(
        brep, &cap_faces, 5.0e-3, 5.0e-2,
    );
    if topology_error.is_none() && continuity_error.is_none() {
        return Some(patches);
    }

    if let Some(diagnostic) =
        crate::construct::cap_outer_acceptance::cap_support_continuity_diagnostic(
            brep, &cap_faces, 5.0e-3, 5.0e-2,
        )
    {
        diagnostics.push(diagnostic);
    }
    let detail_code = topology_error
        .map(|error| format!("topology_validation_failed={error}"))
        .or(continuity_error);
    diagnostics.push(pns_dispatch_declined_diagnostic(
        n_sides,
        "cap_support_acceptance",
        detail_code,
    ));
    *brep = pns_snapshot;
    None
}
