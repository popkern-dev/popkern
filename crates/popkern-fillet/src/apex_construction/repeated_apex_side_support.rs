//! Read-only side-support region planner for repeated-apex trim rewrites.

use popkern_core::topology::BRep;

use super::repeated_apex_region_candidates::CandidateSetSnapshot;
use super::repeated_apex_region_graph::RegionGraph;
use super::repeated_apex_side_support_execution::SideSupportRegionExecutionReport;
use super::repeated_apex_side_support_facts::SideSupportRegionFacts;
use super::repeated_apex_side_support_plan::SideSupportRegionPlanSelection;
use super::repeated_apex_topology_knowledge::TopologyKnowledgeGraph;
use super::repeated_apex_trim_subcomplex::TrimRewriteSubcomplexPlan;

#[derive(Clone, Debug)]
pub(super) struct SideSupportRegionDiagnosticSnapshot {
    knowledge: TopologyKnowledgeGraph,
    facts: SideSupportRegionFacts,
    region_graph: RegionGraph,
    candidates: CandidateSetSnapshot,
    selection: SideSupportRegionPlanSelection,
    execution: SideSupportRegionExecutionReport,
}

impl SideSupportRegionDiagnosticSnapshot {
    pub(super) fn build_read_only(brep: &BRep, plan: &TrimRewriteSubcomplexPlan) -> Self {
        let knowledge = TopologyKnowledgeGraph::build_local(brep, plan);
        let facts = SideSupportRegionFacts::classify(&knowledge);
        let region_graph = RegionGraph::build(&facts);
        let candidates = CandidateSetSnapshot::enumerate(&region_graph);
        let selection = SideSupportRegionPlanSelection::select(&candidates);
        let execution = SideSupportRegionExecutionReport::preflight(&selection, brep, true);
        Self {
            knowledge,
            facts,
            region_graph,
            candidates,
            selection,
            execution,
        }
    }

    pub(super) fn validate_fresh(&self, brep: &BRep) -> bool {
        self.knowledge.validate_fresh(brep)
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "{}; {}; {}; {}; {}; {}",
            self.knowledge.diagnostic_label(),
            self.facts.diagnostic_label(),
            self.region_graph.diagnostic_label(),
            self.candidates.diagnostic_label(),
            self.selection.diagnostic_label(),
            self.execution.diagnostic_label()
        )
    }
}

#[cfg(test)]
mod tests {
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

    use super::super::repeated_apex_region_candidates::{CandidateStatus, RegionOperationKind};
    use super::super::repeated_apex_side_support_plan::{
        SideSupportRegionPlanSelection, SideSupportRegionRefusal, SideSupportRegionRefusalReason,
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
        assert!(snapshot.diagnostic_label().contains("ready_for_executor"));
        assert!(
            snapshot
                .diagnostic_label()
                .contains("total_delta=dV=0 dE=1 dF=1 dL=1")
        );
    }

    #[test]
    fn read_only_snapshot_blocks_bridge_plan_that_leaves_support_family() {
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
            TrimRewriteEdgePlan {
                fillet_index: 1,
                owner_pair: TrimOwnerPair::CylinderSupport {
                    cylinder_face,
                    support_face,
                },
                construction: TrimEdgeConstruction::NewExactGeneratorTrim,
                start: outside_vertex,
                end: support_vertices[1],
            },
        ]);

        let snapshot = SideSupportRegionDiagnosticSnapshot::build_read_only(&brep, &plan);

        let labels = snapshot.diagnostic_label();
        assert!(labels.contains("outside-support-family"));
        assert!(labels.contains("bridge-support-region"));
        assert!(labels.contains("total_delta=dV=0 dE=2 dF=2 dL=2"));
        assert!(labels.contains("dV=0 dE=2 dF=2 dL=2"));
        assert!(labels.contains("one_endpoint_on_support_loop=true"));
        assert!(labels.contains("MissingValidatedFaceFamilyBridgePrimitive"));
        assert!(labels.contains("face-family-bridge primitive_ready=true"));
        assert!(labels.contains("MissingFaceFamilyBridgeTopologyOperator"));
        assert!(labels.contains("topology_operator_ready=false"));
        assert!(labels.contains("boundary_cycles=2"));
        assert!(labels.contains("missing-cross-face-bridge-euler-operator"));
        assert!(labels.contains("execution refused stage=Preflight"));
        assert_eq!(snapshot.candidates.candidates.len(), 4);
        assert!(snapshot.candidates.candidates.iter().any(|candidate| {
            candidate.kind == RegionOperationKind::BridgeSupportRegion
                && matches!(candidate.status, CandidateStatus::Accepted)
        }));
        assert!(matches!(
            snapshot.selection,
            SideSupportRegionPlanSelection::Blocked(
                _,
                SideSupportRegionRefusal {
                    reason:
                        SideSupportRegionRefusalReason::MissingValidatedFaceFamilyBridgePrimitive,
                    ..
                }
            )
        ));
    }

    #[test]
    fn read_only_snapshot_refuses_bridge_region_with_one_trim() {
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
        let plan = trim_plan(vec![TrimRewriteEdgePlan {
            fillet_index: 0,
            owner_pair: TrimOwnerPair::CylinderSupport {
                cylinder_face,
                support_face,
            },
            construction: TrimEdgeConstruction::NewExactGeneratorTrim,
            start: support_vertices[0],
            end: outside_vertex,
        }]);

        let snapshot = SideSupportRegionDiagnosticSnapshot::build_read_only(&brep, &plan);

        let labels = snapshot.diagnostic_label();
        assert!(labels.contains("UnsupportedBridgeTrimCount"));
        assert!(!labels.contains("MissingValidatedFaceFamilyBridgePrimitive"));
        assert_eq!(snapshot.candidates.candidates.len(), 2);
        assert!(matches!(
            snapshot.selection,
            SideSupportRegionPlanSelection::Refused(SideSupportRegionRefusal {
                reason: SideSupportRegionRefusalReason::UnsupportedBridgeTrimCount,
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
