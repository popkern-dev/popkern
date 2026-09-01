//! The setback vertex-blend rung of `insert_vertex_blend` (PK-094
//! surgery entry, PK-225 N=3 production gate). Extracted verbatim
//! from `blend_insert.rs` for the structural file budget; see the
//! call site for the rung ordering contract.

use super::InsertVertexBlendOutput;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};

use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod, CornerPatchInfo};
use crate::vertex_classifier::VertexBlendClassification;

/// Attempt the setback boundary surgery. `Some(result)` when the rung
/// decided (served or hard-errored); `None` falls through to the
/// vertex-blend corner machinery with the decline recorded.
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
) -> Option<Result<InsertVertexBlendOutput, String>> {
    if crate::construct::setback_vertex_topology::setback_topology_commit_enabled()
        && edge_indices.len() >= setback_min_n
    {
        // PK-094 commit mode: the setback vertex blend replaces the
        // vertex-blend corner machinery entirely for eligible apexes. On
        // any surgery error the snapshot restore inside the surgery
        // has already run; fall through to the vertex-blend path with a
        // recorded reason.
        let commit_snapshot = brep.clone();
        // The KP20 cap serves BOTH continuity levels: it is the
        // proven watertight tiling and a G2 cap satisfies G1 (the
        // PnS bi-4 family integrated at its historical 5e-2 seam class —
        // its patches overlap by construction, see PK-097).
        let gate_seam_g2 = requested_continuity == ContinuityLevel::G2;
        let mut fan_downgraded_to_g1 = false;
        let mut surgery = crate::construct::setback_vertex_topology::try_setback_boundary_surgery(
            brep,
            &commit_snapshot,
            vertex,
            fillet_surfaces_pre,
            tolerance,
            gate_seam_g2,
            cap_presolve,
        );
        // A fan that fails ONLY the seam-G2 postcondition still beats
        // every fallback by orders of magnitude (measured on the
        // rim+verticals box class: the declined fan's seams sat at
        // g1 2.9e-3 / g2 5.4e-2 while the Gregory fallback served a
        // 27° crease). Serve the fan as a measured G1 downgrade: the
        // corner patch records achieved G1 and the quality pass emits
        // ContinuityDowngraded.
        if gate_seam_g2
            && matches!(&surgery, Err(e) if e.starts_with("fan seam G2 postcondition failed"))
        {
            fan_downgraded_to_g1 = true;
            surgery = crate::construct::setback_vertex_topology::try_setback_boundary_surgery(
                brep,
                &commit_snapshot,
                vertex,
                fillet_surfaces_pre,
                tolerance,
                false,
                cap_presolve,
            );
        }
        match surgery {
            Ok(outcome) => {
                if let Err(e) = brep.validate() {
                    *brep = commit_snapshot;
                    *crate::construct::setback_vertex_topology::SETBACK_TOPOLOGY_LAST
                        .lock()
                        .unwrap_or_else(|p| p.into_inner()) =
                        Some(format!("setback commit failed validation: {e}"));
                    return Some(Err(format!("setback vertex blend failed validation: {e}")));
                }
                {
                    let infos: Vec<CornerPatchInfo> = outcome
                        .cap_faces
                        .iter()
                        .map(|&face| CornerPatchInfo {
                            vertex,
                            n_sides: edge_indices.len(),
                            method: CornerMethod::Gregory,
                            // The blend's continuity is capped by the
                            // strips: a KP20 (G2-interior) cap on G1
                            // sweeps yields a G1 blend and must say so.
                            achieved_continuity: match requested_continuity {
                                ContinuityLevel::G2 if !fan_downgraded_to_g1 => {
                                    crate::types::ContinuityLevel::G2
                                }
                                _ => crate::types::ContinuityLevel::G1,
                            },
                            classification,
                            face,
                            bridge_fallbacks: Vec::new(),
                            overlap_remediated: false,
                            relaxed_continuity_coplanarity: None,
                            cap_sphere: None,
                            sub_patch_provenance: None,
                            setback_cap: true,
                            exact_boundary_cap: false,
                        })
                        .collect();
                    return Some(Ok((infos, None, Vec::new())));
                }
            }
            Err(e) => {
                // Fall through to the vertex-blend path.
                *crate::construct::setback_vertex_topology::SETBACK_TOPOLOGY_LAST
                    .lock()
                    .unwrap_or_else(|p| p.into_inner()) =
                    Some(format!("setback commit declined: {e}"));
                if popkern_core::trace::trace_enabled("corner-carve") {
                    eprintln!("  [setback] declined at {vertex:?}: {e}");
                }
            }
        }
    }
    None
}
