//! Read-only side-support region planner for repeated-apex trim rewrites.

use popkern_core::topology::BRep;

use super::repeated_apex_region_candidates::{
    CandidateRejection, CandidateSetSnapshot, CandidateStatus, RegionOperationKind,
};
use super::repeated_apex_region_graph::RegionGraph;
use super::repeated_apex_topology_knowledge::TopologyKnowledgeGraph;
use super::repeated_apex_trim_subcomplex::TrimRewriteSubcomplexPlan;

#[derive(Clone, Debug)]
pub(super) struct SideSupportRegionDiagnosticSnapshot {
    knowledge: TopologyKnowledgeGraph,
    region_graph: RegionGraph,
    candidates: CandidateSetSnapshot,
    selection: SideSupportRegionPlanSelection,
}

#[derive(Clone, Debug)]
pub(super) enum SideSupportRegionPlanSelection {
    Selected(SelectedSideSupportRegionPlan),
    Refused(SideSupportRegionRefusal),
}

#[derive(Clone, Debug)]
pub(super) struct SelectedSideSupportRegionPlan {
    operations: Vec<PlannedTopologyOperation>,
}

#[derive(Clone, Debug)]
pub(super) enum PlannedTopologyOperation {
    SplitFaceWithEdge { fillet_index: usize },
}

#[derive(Clone, Debug)]
pub(super) struct SideSupportRegionRefusal {
    stage: PlannerStage,
    reason: SideSupportRegionRefusalReason,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum PlannerStage {
    SelectPlan,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum SideSupportRegionRefusalReason {
    UnknownEndpointMembership,
    MissingValidatedTopologyPrimitive,
    NoAcceptedOperationCandidate,
}

impl SideSupportRegionDiagnosticSnapshot {
    pub(super) fn build_read_only(brep: &BRep, plan: &TrimRewriteSubcomplexPlan) -> Self {
        let knowledge = TopologyKnowledgeGraph::build_local(brep, plan);
        let region_graph = RegionGraph::build(&knowledge);
        let candidates = CandidateSetSnapshot::enumerate(&region_graph);
        let selection = SideSupportRegionPlanSelection::select(&candidates);
        Self {
            knowledge,
            region_graph,
            candidates,
            selection,
        }
    }

    pub(super) fn validate_fresh(&self, brep: &BRep) -> bool {
        self.knowledge.validate_fresh(brep)
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "{}; {}; {}; {}",
            self.knowledge.diagnostic_label(),
            self.region_graph.diagnostic_label(),
            self.candidates.diagnostic_label(),
            self.selection.diagnostic_label()
        )
    }
}

impl SideSupportRegionPlanSelection {
    fn select(candidates: &CandidateSetSnapshot) -> Self {
        let mut operations = Vec::new();
        let mut saw_unknown_endpoint = false;
        let mut saw_missing_primitive = false;
        for candidate in &candidates.candidates {
            match (&candidate.kind, &candidate.status) {
                (RegionOperationKind::SplitFaceWithEdge, CandidateStatus::Accepted) => {
                    operations.push(PlannedTopologyOperation::SplitFaceWithEdge {
                        fillet_index: candidate.fillet_index,
                    });
                }
                (
                    _,
                    CandidateStatus::Rejected {
                        reason: CandidateRejection::EndpointMembershipUnknown,
                    },
                ) => saw_unknown_endpoint = true,
                (
                    RegionOperationKind::RequiresSupportRegionBridge,
                    CandidateStatus::Rejected {
                        reason: CandidateRejection::MissingValidatedBridgePrimitive,
                    },
                ) => saw_missing_primitive = true,
                _ => {}
            }
        }
        if !operations.is_empty() && !saw_unknown_endpoint && !saw_missing_primitive {
            return SideSupportRegionPlanSelection::Selected(SelectedSideSupportRegionPlan {
                operations,
            });
        }
        let reason = if saw_unknown_endpoint {
            SideSupportRegionRefusalReason::UnknownEndpointMembership
        } else if saw_missing_primitive {
            SideSupportRegionRefusalReason::MissingValidatedTopologyPrimitive
        } else {
            SideSupportRegionRefusalReason::NoAcceptedOperationCandidate
        };
        SideSupportRegionPlanSelection::Refused(SideSupportRegionRefusal {
            stage: PlannerStage::SelectPlan,
            reason,
        })
    }

    fn diagnostic_label(&self) -> String {
        match self {
            SideSupportRegionPlanSelection::Selected(plan) => plan.diagnostic_label(),
            SideSupportRegionPlanSelection::Refused(refusal) => refusal.diagnostic_label(),
        }
    }
}

impl SelectedSideSupportRegionPlan {
    fn diagnostic_label(&self) -> String {
        format!(
            "selected_plan operations=[{}]",
            self.operations
                .iter()
                .map(PlannedTopologyOperation::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }
}

impl PlannedTopologyOperation {
    fn diagnostic_label(&self) -> String {
        match self {
            PlannedTopologyOperation::SplitFaceWithEdge { fillet_index } => {
                format!("FilletEnd[{fillet_index}] split-face-with-edge")
            }
        }
    }
}

impl SideSupportRegionRefusal {
    fn diagnostic_label(&self) -> String {
        format!("refusal stage={:?} reason={:?}", self.stage, self.reason)
    }
}

#[cfg(test)]
mod tests {
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

    use super::super::repeated_apex_region_candidates::{
        CandidateRejection, CandidateStatus, RegionOperationKind,
    };
    use super::super::repeated_apex_trim_subcomplex::{
        TopologyCountDelta, TrimEdgeConstruction, TrimOwnerPair, TrimRewriteEdgePlan,
        TrimRewriteSubcomplexPlan, trim_rewrite_support_groups,
    };
    use super::*;

    #[test]
    fn read_only_snapshot_accepts_existing_support_loop_split_candidate() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let faces = brep.faces.keys().collect::<Vec<_>>();
        let support_face = faces[0];
        let cylinder_face = faces[1];
        let support_vertices = outer_loop_vertices(&brep, support_face);
        let plan = trim_plan(vec![TrimRewriteEdgePlan {
            fillet_index: 0,
            owner_pair: TrimOwnerPair::CylinderSupport {
                cylinder_face,
                support_face,
            },
            construction: TrimEdgeConstruction::NewExactGeneratorTrim,
            start: support_vertices[0],
            end: support_vertices[1],
        }]);

        let snapshot = SideSupportRegionDiagnosticSnapshot::build_read_only(&brep, &plan);

        assert!(snapshot.validate_fresh(&brep));
        assert_eq!(snapshot.candidates.candidates.len(), 1);
        let candidate = &snapshot.candidates.candidates[0];
        assert_eq!(candidate.kind, RegionOperationKind::SplitFaceWithEdge);
        assert!(matches!(candidate.status, CandidateStatus::Accepted));
        assert!(matches!(
            snapshot.selection,
            SideSupportRegionPlanSelection::Selected(_)
        ));
        assert!(snapshot.diagnostic_label().contains("selected_plan"));
    }

    #[test]
    fn read_only_snapshot_rejects_generator_trim_that_leaves_support_family() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let faces = brep.faces.keys().collect::<Vec<_>>();
        let support_face = faces[0];
        let cylinder_face = faces[1];
        let cap_face = face_with_vertex_outside_support_loop(&brep, support_face);
        let support_vertices = outer_loop_vertices(&brep, support_face);
        let cap_vertices = outer_loop_vertices(&brep, cap_face);
        let outside_vertex = cap_vertices
            .into_iter()
            .find(|vertex| !support_vertices.contains(vertex))
            .expect("fixture has a cap vertex outside support loop");
        let plan = trim_plan(vec![
            TrimRewriteEdgePlan {
                fillet_index: 0,
                owner_pair: TrimOwnerPair::CylinderCap {
                    cylinder_face,
                    cap_face,
                },
                construction: TrimEdgeConstruction::ExistingSeam,
                start: outside_vertex,
                end: outside_vertex,
            },
            TrimRewriteEdgePlan {
                fillet_index: 0,
                owner_pair: TrimOwnerPair::CylinderSupport {
                    cylinder_face,
                    support_face,
                },
                construction: TrimEdgeConstruction::NewExactGeneratorTrim,
                start: support_vertices[0],
                end: outside_vertex,
            },
        ]);

        let snapshot = SideSupportRegionDiagnosticSnapshot::build_read_only(&brep, &plan);

        let labels = snapshot.diagnostic_label();
        assert!(labels.contains("outside-support-family"));
        assert!(labels.contains("missing-validated-bridge-primitive"));
        assert_eq!(snapshot.candidates.candidates.len(), 2);
        assert!(snapshot.candidates.candidates.iter().any(|candidate| {
            candidate.kind == RegionOperationKind::RequiresSupportRegionBridge
                && matches!(
                    candidate.status,
                    CandidateStatus::Rejected {
                        reason: CandidateRejection::MissingValidatedBridgePrimitive
                    }
                )
        }));
        assert!(matches!(
            snapshot.selection,
            SideSupportRegionPlanSelection::Refused(SideSupportRegionRefusal {
                reason: SideSupportRegionRefusalReason::MissingValidatedTopologyPrimitive,
                ..
            })
        ));
    }

    fn trim_plan(edge_plans: Vec<TrimRewriteEdgePlan>) -> TrimRewriteSubcomplexPlan {
        TrimRewriteSubcomplexPlan {
            support_groups: trim_rewrite_support_groups(&edge_plans),
            edge_plans,
            expected_edge_only_delta: TopologyCountDelta::zero(),
        }
    }

    fn face_with_vertex_outside_support_loop(brep: &BRep, support_face: FaceHandle) -> FaceHandle {
        let support_vertices = outer_loop_vertices(brep, support_face);
        brep.faces
            .keys()
            .find(|face| {
                *face != support_face
                    && outer_loop_vertices(brep, *face)
                        .iter()
                        .any(|vertex| !support_vertices.contains(vertex))
            })
            .expect("fixture has a face with a vertex outside support loop")
    }

    fn outer_loop_vertices(brep: &BRep, face: FaceHandle) -> Vec<VertexHandle> {
        let loop_handle = brep.faces[face].outer_loop;
        let start = brep.loops[loop_handle].half_edge;
        let mut vertices = Vec::new();
        let mut cursor = start;
        loop {
            let half_edge = &brep.half_edges[cursor];
            vertices.push(half_edge.origin);
            cursor = half_edge.next;
            if cursor == start {
                break;
            }
        }
        vertices
    }
}
