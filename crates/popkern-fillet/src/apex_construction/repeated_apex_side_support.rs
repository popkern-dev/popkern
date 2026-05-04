//! Read-only side-support region planner for repeated-apex trim rewrites.

use popkern_core::topology::BRep;

use super::repeated_apex_region_candidates::CandidateSetSnapshot;
use super::repeated_apex_region_graph::RegionGraph;
use super::repeated_apex_topology_knowledge::TopologyKnowledgeGraph;
use super::repeated_apex_trim_subcomplex::TrimRewriteSubcomplexPlan;

#[derive(Clone, Debug)]
pub(super) struct SideSupportRegionDiagnosticSnapshot {
    knowledge: TopologyKnowledgeGraph,
    region_graph: RegionGraph,
    candidates: CandidateSetSnapshot,
}

impl SideSupportRegionDiagnosticSnapshot {
    pub(super) fn build_read_only(brep: &BRep, plan: &TrimRewriteSubcomplexPlan) -> Self {
        let knowledge = TopologyKnowledgeGraph::build_local(brep, plan);
        let region_graph = RegionGraph::build(&knowledge);
        let candidates = CandidateSetSnapshot::enumerate(&region_graph);
        Self {
            knowledge,
            region_graph,
            candidates,
        }
    }

    pub(super) fn validate_fresh(&self, brep: &BRep) -> bool {
        self.knowledge.validate_fresh(brep)
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "{}; {}; {}",
            self.knowledge.diagnostic_label(),
            self.region_graph.diagnostic_label(),
            self.candidates.diagnostic_label()
        )
    }
}
