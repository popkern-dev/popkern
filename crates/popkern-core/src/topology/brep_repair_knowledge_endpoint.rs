//! Shared endpoint membership facts for BRep repair knowledge.

use crate::topology::{BrepLoopKind, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

/// Shared facts describing where a repair endpoint vertex appears in loops.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BrepEndpointMembershipFacts {
    /// Endpoint vertex.
    pub vertex: VertexHandle,
    /// Loop memberships containing the endpoint vertex.
    pub containing_loops: Vec<BrepEndpointLoopMembership>,
}

impl BrepEndpointMembershipFacts {
    /// Loop memberships containing the endpoint vertex.
    pub fn containing_loops(&self) -> &[BrepEndpointLoopMembership] {
        &self.containing_loops
    }

    /// Migration-compatible compact diagnostic label (pre-render-key scheme).
    pub fn diagnostic_label(&self) -> String {
        format!(
            "{:?} in [{}]",
            self.vertex,
            self.containing_loops
                .iter()
                .map(BrepEndpointLoopMembership::diagnostic_label)
                .collect::<Vec<_>>()
                .join(", ")
        )
    }
}

/// Shared facts for one loop occurrence of an endpoint vertex.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BrepEndpointLoopMembership {
    /// Face containing the endpoint.
    pub face: FaceHandle,
    /// Loop containing the endpoint.
    pub loop_handle: LoopHandle,
    /// Half-edge whose origin is the endpoint vertex.
    pub half_edge: HalfEdgeHandle,
    /// Loop role on the owning face.
    pub kind: BrepLoopKind,
}

impl BrepEndpointLoopMembership {
    /// Face containing the endpoint.
    pub fn face(&self) -> FaceHandle {
        self.face
    }

    /// Loop containing the endpoint.
    pub fn loop_handle(&self) -> LoopHandle {
        self.loop_handle
    }

    /// Half-edge whose origin is the endpoint vertex.
    pub fn half_edge(&self) -> HalfEdgeHandle {
        self.half_edge
    }

    /// Loop role on the owning face.
    pub fn kind(&self) -> BrepLoopKind {
        self.kind
    }

    /// Migration-compatible compact diagnostic label (pre-render-key scheme).
    pub fn diagnostic_label(&self) -> String {
        format!(
            "{:?}/{:?}/{:?}@{:?}",
            self.face, self.loop_handle, self.kind, self.half_edge
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn endpoint_membership_label_contains_legacy_fields() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let face = brep.faces.keys().next().expect("cube has face");
        let loop_handle = brep.faces[face].outer_loop;
        let half_edge = brep.loops[loop_handle].half_edge;
        let vertex = brep.half_edges[half_edge].origin;

        let facts = BrepEndpointMembershipFacts {
            vertex,
            containing_loops: vec![BrepEndpointLoopMembership {
                face,
                loop_handle,
                half_edge,
                kind: BrepLoopKind::Outer,
            }],
        };

        let label = facts.diagnostic_label();
        assert!(label.contains(" in ["), "{label}");
        assert!(label.contains("Outer"), "{label}");
    }
}
