//! Read-only side-support region planner for repeated-apex trim rewrites.

use popkern_core::topology::BRep;

use super::error::ApexBuildError;
use super::repeated_apex_region_candidates::CandidateSetSnapshot;
use super::repeated_apex_region_graph::RegionGraph;
use super::repeated_apex_side_support_bridge::FaceFamilyBridgePrimitiveDerivation;
use super::repeated_apex_side_support_execution::SideSupportRegionExecutionReport;
use super::repeated_apex_side_support_facts::SideSupportRegionFacts;
use super::repeated_apex_side_support_merge::SupportRegionMergePrimitiveDerivation;
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

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum SideSupportRegionContinuationOutcome {
    Applied,
    NotReady,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum SideSupportRegionContinuationFailure {
    StaleKnowledgeGraph,
    BridgeTopologyOperator { reason: String },
    MergeTopologyAndGeometry { reason: String },
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

    fn apply_executable_continuation(
        &self,
        brep: &mut BRep,
    ) -> Result<SideSupportRegionContinuationOutcome, SideSupportRegionContinuationFailure> {
        if !self.validate_fresh(brep) {
            return Err(SideSupportRegionContinuationFailure::StaleKnowledgeGraph);
        }
        let bridge_primitive = FaceFamilyBridgePrimitiveDerivation::build(&self.selection, brep);
        let topology_report = bridge_primitive.topology_operator_report(brep);
        if bridge_primitive.is_complete() && topology_report.is_ready() {
            let rollback = brep.clone();
            if let Err(err) = bridge_primitive.apply_topology_operator(brep) {
                *brep = rollback;
                return Err(
                    SideSupportRegionContinuationFailure::BridgeTopologyOperator { reason: err },
                );
            }
            return Ok(SideSupportRegionContinuationOutcome::Applied);
        }
        let merge_primitive = SupportRegionMergePrimitiveDerivation::build(&self.selection, brep);
        if !merge_primitive.exact_geometry_ready() {
            return Ok(SideSupportRegionContinuationOutcome::NotReady);
        }
        let rollback = brep.clone();
        if let Err(err) = merge_primitive.apply_topology_and_geometry(brep) {
            *brep = rollback;
            return Err(
                SideSupportRegionContinuationFailure::MergeTopologyAndGeometry { reason: err },
            );
        }
        Ok(SideSupportRegionContinuationOutcome::Applied)
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

pub(super) fn apply_selected_side_support_continuation(
    brep: &mut BRep,
    subcomplex_plan: &TrimRewriteSubcomplexPlan,
) -> Result<SideSupportRegionContinuationOutcome, ApexBuildError> {
    let diagnostics = SideSupportRegionDiagnosticSnapshot::build_read_only(brep, subcomplex_plan);
    diagnostics
        .apply_executable_continuation(brep)
        .map_err(|failure| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition target-cylinder trim rewrite could not apply \
                 executable side-support continuation: {}. Side-support diagnostics: {}",
                failure.diagnostic_label(),
                diagnostics.diagnostic_label()
            ),
        })
}

impl SideSupportRegionContinuationOutcome {
    pub(super) fn diagnostic_label(&self) -> &'static str {
        match self {
            SideSupportRegionContinuationOutcome::Applied => "applied",
            SideSupportRegionContinuationOutcome::NotReady => "not-ready",
        }
    }
}

impl SideSupportRegionContinuationFailure {
    pub(super) fn diagnostic_label(&self) -> String {
        match self {
            SideSupportRegionContinuationFailure::StaleKnowledgeGraph => {
                "stale-knowledge-graph".to_string()
            }
            SideSupportRegionContinuationFailure::BridgeTopologyOperator { reason } => {
                format!("bridge-topology-operator-failed reason={reason}")
            }
            SideSupportRegionContinuationFailure::MergeTopologyAndGeometry { reason } => {
                format!("merge-topology-and-geometry-failed reason={reason}")
            }
        }
    }
}

#[cfg(test)]
#[path = "repeated_apex_side_support_tests.rs"]
mod tests;
