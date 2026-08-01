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
/// legacy corner machinery with the decline recorded.
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
) -> Option<Result<InsertVertexBlendOutput, String>> {
    if crate::construct::setback_vertex_topology::setback_topology_commit_enabled()
        && edge_indices.len() >= setback_min_n
    {
        // PK-094 commit mode: the setback vertex blend replaces the
        // legacy corner machinery entirely for eligible apexes. On
        // any surgery error the snapshot restore inside the surgery
        // has already run; fall through to the legacy path with a
        // recorded reason.
        let commit_snapshot = brep.clone();
        // The KP20 cap is the default for BOTH continuity levels: it
        // is the proven watertight tiling and a G2 cap satisfies G1.
        // The PnS bi-4 family stays reachable as a diagnostic /
        // comparison path (POPKERN_G1_PNS_CAP=1); it integrates at
        // the legacy PnS quality class (5e-2 seams — its patches
        // overlap by construction, see PK-097).
        let family = match requested_continuity {
            ContinuityLevel::G2 => crate::construct::setback_vertex_topology::CapFamily::Kp20G2,
            _ if std::env::var_os("POPKERN_G1_PNS_CAP").is_some() => {
                crate::construct::setback_vertex_topology::CapFamily::PnsG1
            }
            _ => crate::construct::setback_vertex_topology::CapFamily::Kp20G2,
        };
        match crate::construct::setback_vertex_topology::try_setback_boundary_surgery(
            brep,
            vertex,
            fillet_surfaces_pre,
            tolerance,
            family,
            requested_continuity == ContinuityLevel::G2,
        ) {
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
                            method: match family {
                                crate::construct::setback_vertex_topology::CapFamily::Kp20G2 => {
                                    CornerMethod::Gregory
                                }
                                crate::construct::setback_vertex_topology::CapFamily::PnsG1 => {
                                    CornerMethod::PolyhedralNetSpline
                                }
                            },
                            // The blend's continuity is capped by the
                            // strips: a KP20 (G2-interior) cap on G1
                            // sweeps yields a G1 blend and must say so.
                            achieved_continuity: match (family, requested_continuity) {
                                (
                                    crate::construct::setback_vertex_topology::CapFamily::Kp20G2,
                                    ContinuityLevel::G2,
                                ) => crate::types::ContinuityLevel::G2,
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
                // Fall through to the legacy path.
                *crate::construct::setback_vertex_topology::SETBACK_TOPOLOGY_LAST
                    .lock()
                    .unwrap_or_else(|p| p.into_inner()) =
                    Some(format!("setback commit declined: {e}"));
                if std::env::var_os("PK_CARVE_TRACE").is_some() {
                    eprintln!("  [setback] declined at {vertex:?}: {e}");
                }
            }
        }
    }
    None
}
