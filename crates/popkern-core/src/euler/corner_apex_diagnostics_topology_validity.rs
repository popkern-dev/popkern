//! Topology-validity facts for corner/apex diagnostics.

use super::types::*;
use crate::geometry::surface::SurfaceGeometry;
use crate::topology::{BRep, EdgeHandle, FaceHandle, TopologyRole};

pub(super) fn extract_topology_validity(
    brep: &BRep,
    incident_edges: &[EdgeHandle],
    incident_faces: &[FaceHandle],
) -> CornerTopologyValidityFact {
    CornerTopologyValidityFact {
        body_shell: body_shell_membership(brep, incident_faces),
        boundary: boundary_membership(brep, incident_edges),
        radial: radial_topology(brep, incident_edges),
        generated_context: generated_context(brep, incident_edges, incident_faces),
    }
}

pub(super) fn missing_topology_validity() -> CornerTopologyValidityFact {
    CornerTopologyValidityFact::default()
}

fn body_shell_membership(brep: &BRep, faces: &[FaceHandle]) -> CornerBodyShellMembership {
    if faces.is_empty() {
        return CornerBodyShellMembership::UnknownInputFact;
    }
    let mut shells = std::collections::BTreeSet::new();
    for face in faces {
        let shell = brep.faces[*face].shell;
        let Some(shell_record) = brep.shells.get(shell) else {
            return CornerBodyShellMembership::MissingShell;
        };
        if !brep.solids.contains_key(shell_record.solid) {
            return CornerBodyShellMembership::MissingSolid;
        }
        shells.insert(shell);
    }
    if shells.len() == 1 {
        CornerBodyShellMembership::SameBodyShell
    } else {
        CornerBodyShellMembership::MultipleShells
    }
}

fn boundary_membership(brep: &BRep, edges: &[EdgeHandle]) -> CornerBoundaryMembership {
    if edges.is_empty() {
        return CornerBoundaryMembership::WireOrSheet;
    }
    if edges.iter().any(|edge| brep.boundary_edges.contains(edge)) {
        CornerBoundaryMembership::Boundary
    } else {
        CornerBoundaryMembership::Interior
    }
}

fn radial_topology(brep: &BRep, edges: &[EdgeHandle]) -> CornerRadialTopology {
    if edges.is_empty() {
        return CornerRadialTopology::UnknownInputFact;
    }
    if edges.iter().any(|edge| {
        let half_edge = brep.edges[*edge].half_edge;
        brep.nonmanifold_radial.contains_key(&half_edge)
            || brep
                .half_edges
                .get(half_edge)
                .is_some_and(|he| brep.nonmanifold_radial.contains_key(&he.twin))
    }) {
        CornerRadialTopology::NonManifoldRadial
    } else {
        CornerRadialTopology::Manifold
    }
}

fn generated_context(
    brep: &BRep,
    edges: &[EdgeHandle],
    faces: &[FaceHandle],
) -> CornerGeneratedTopologyContext {
    if edges
        .iter()
        .any(|edge| brep.edges[*edge].role == TopologyRole::ConstructionResidual)
        || faces
            .iter()
            .any(|face| brep.faces[*face].role == TopologyRole::ConstructionResidual)
    {
        return CornerGeneratedTopologyContext::ConstructionResidual;
    }
    if faces.iter().any(|face| {
        brep.faces[*face]
            .surface
            .as_ref()
            .is_some_and(|surface| matches!(surface.geometry, SurfaceGeometry::Gregory(_)))
    }) {
        CornerGeneratedTopologyContext::ExistingBlend
    } else {
        CornerGeneratedTopologyContext::Production
    }
}
