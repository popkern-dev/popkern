//! Read-only ownership classification for BRep export-fidelity blockers.

use crate::topology::brep_edge_surface_compatibility_diagnostics::surface_kind;
use crate::topology::{BRep, BrepEdgeSurfaceKind, EdgeHandle, FaceHandle, TopologyRole};

/// Ownership classification for pcurve/edge-surface fidelity blockers.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BrepFidelityBlockerOwnershipReport {
    /// Number of inspected blocker edges.
    pub blockers: usize,
    /// Blocker edges with production edge ownership and at least one production
    /// face-side blocker.
    pub production_owned: usize,
    /// Blocker edges involving construction-residual edge or face ownership.
    pub construction_residual_owned: usize,
    /// Per-edge ownership facts.
    pub edges: Vec<BrepFidelityBlockerOwnership>,
}

/// Per-edge ownership fact for an export-fidelity blocker.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepFidelityBlockerOwnership {
    /// Edge carrying the blocker.
    pub edge: EdgeHandle,
    /// Edge semantic ownership role.
    pub edge_role: TopologyRole,
    /// Maximum edge/pcurve deviation for this edge.
    pub max_pcurve_deviation: f64,
    /// Face-side measurements contributing to this blocker.
    pub half_edges: Vec<BrepFidelityBlockerSideOwnership>,
}

/// Per-half-edge face-side ownership fact for an export-fidelity blocker.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepFidelityBlockerSideOwnership {
    /// Owning face.
    pub face: FaceHandle,
    /// Face semantic ownership role.
    pub face_role: TopologyRole,
    /// Surface kind on the owning face.
    pub surface_kind: Option<BrepEdgeSurfaceKind>,
    /// Maximum edge/pcurve deviation on this side.
    pub max_pcurve_deviation: f64,
}

/// Classify export-fidelity blockers by semantic ownership.
///
/// This read-only diagnostic separates production geometry defects from
/// construction-residual leakage so repair and export decisions can act on
/// model-owned topology instead of debug trace text.
pub fn classify_brep_fidelity_blocker_ownership(
    brep: &BRep,
    tolerance: f64,
) -> BrepFidelityBlockerOwnershipReport {
    let mut report = BrepFidelityBlockerOwnershipReport::default();
    for edge in brep.edges.keys() {
        let Some(deviation) = brep.measure_edge_pcurve_deviation(edge) else {
            continue;
        };
        if deviation.max_deviation <= tolerance {
            continue;
        }
        let Some(edge_data) = brep.edges.get(edge) else {
            continue;
        };
        let half_edges = deviation
            .half_edges
            .iter()
            .filter_map(|side| {
                let face = brep.faces.get(side.face)?;
                Some(BrepFidelityBlockerSideOwnership {
                    face: side.face,
                    face_role: face.role,
                    surface_kind: face
                        .surface
                        .as_ref()
                        .map(|surface| surface_kind(&surface.geometry)),
                    max_pcurve_deviation: side.max_deviation,
                })
            })
            .collect::<Vec<_>>();
        let production_side = half_edges
            .iter()
            .any(|side| side.face_role == TopologyRole::Production);
        let residual_side = half_edges
            .iter()
            .any(|side| side.face_role == TopologyRole::ConstructionResidual);
        if edge_data.role == TopologyRole::Production && production_side {
            report.production_owned += 1;
        }
        if edge_data.role == TopologyRole::ConstructionResidual || residual_side {
            report.construction_residual_owned += 1;
        }
        report.blockers += 1;
        report.edges.push(BrepFidelityBlockerOwnership {
            edge,
            edge_role: edge_data.role,
            max_pcurve_deviation: deviation.max_deviation,
            half_edges,
        });
    }
    report.edges.sort_by(|a, b| {
        b.max_pcurve_deviation
            .total_cmp(&a.max_pcurve_deviation)
            .then_with(|| format!("{:?}", a.edge).cmp(&format!("{:?}", b.edge)))
    });
    report
}
