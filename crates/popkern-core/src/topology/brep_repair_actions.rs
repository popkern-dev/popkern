//! Generic repair-action diagnostics shared by BRep repair pipelines.

use crate::topology::{BrepRepairFactId, BrepRepairPipelineStage, BrepRepairTopologyHandle};
use serde::{Deserialize, Serialize};

/// Aggregate repair-action diagnostics for API and GUI callers.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepRepairActionReport {
    /// Number of repair actions considered.
    pub attempted_actions: usize,
    /// Number of actions that mutated the repaired clone.
    pub authored_actions: usize,
    /// Number of actions declined before mutation.
    pub declined_actions: usize,
    /// Number of transaction rollbacks.
    pub rolled_back_actions: usize,
    /// Maximum vertex or endpoint motion reported by any action.
    pub max_motion: f64,
    /// Maximum sampled geometric deviation reported by any action.
    pub max_deviation: f64,
    /// Per-action decisions.
    pub decisions: Vec<BrepRepairActionDecision>,
}

impl Default for BrepRepairActionReport {
    fn default() -> Self {
        Self {
            attempted_actions: 0,
            authored_actions: 0,
            declined_actions: 0,
            rolled_back_actions: 0,
            max_motion: 0.0,
            max_deviation: 0.0,
            decisions: Vec::new(),
        }
    }
}

impl BrepRepairActionReport {
    /// Append one action decision and update aggregate counters.
    pub(crate) fn push(&mut self, decision: BrepRepairActionDecision) {
        self.attempted_actions += 1;
        match decision.outcome {
            BrepRepairActionOutcome::Authored => self.authored_actions += 1,
            BrepRepairActionOutcome::Declined => self.declined_actions += 1,
            BrepRepairActionOutcome::RolledBack => self.rolled_back_actions += 1,
        }
        self.max_motion = self.max_motion.max(decision.max_motion);
        self.max_deviation = self.max_deviation.max(decision.max_deviation);
        self.decisions.push(decision);
    }
}

/// Serialization-ready repair-action artifact for API, UI, and debug tooling.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BrepRepairActionArtifact {
    /// Aggregate action counters and worst observed distances.
    pub summary: BrepRepairActionSummary,
    /// Per-action records with stable stage, outcome, and reason labels.
    pub actions: Vec<BrepRepairActionRecord>,
}

/// Serialized aggregate repair-action counters.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BrepRepairActionSummary {
    /// Number of repair actions considered.
    pub attempted_actions: usize,
    /// Number of actions that mutated the repaired clone.
    pub authored_actions: usize,
    /// Number of actions declined before mutation.
    pub declined_actions: usize,
    /// Number of transaction rollbacks.
    pub rolled_back_actions: usize,
    /// Maximum vertex or endpoint motion reported by any action.
    pub max_motion: f64,
    /// Maximum sampled geometric deviation reported by any action.
    pub max_deviation: f64,
}

/// Serialized per-action repair decision.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BrepRepairActionRecord {
    /// Stable pipeline stage label.
    pub stage: String,
    /// Stable outcome label.
    pub outcome: String,
    /// Stable subject label, usually a provenance label or entity handle.
    pub subject: String,
    /// Stage-specific stable reason label for declined or rolled-back actions.
    pub reason: Option<String>,
    /// Vertex or endpoint motion associated with this action.
    pub max_motion: f64,
    /// Sampled geometric deviation associated with this action.
    pub max_deviation: f64,
    /// Topology handles directly affected or considered by this action.
    pub affected_handles: Vec<String>,
    /// Shared repair knowledge facts referenced by this action.
    pub shared_fact_ids: Vec<usize>,
}

/// Build a serialization-ready action artifact from a generic report.
pub fn build_brep_repair_action_artifact(
    report: &BrepRepairActionReport,
) -> BrepRepairActionArtifact {
    BrepRepairActionArtifact {
        summary: BrepRepairActionSummary {
            attempted_actions: report.attempted_actions,
            authored_actions: report.authored_actions,
            declined_actions: report.declined_actions,
            rolled_back_actions: report.rolled_back_actions,
            max_motion: report.max_motion,
            max_deviation: report.max_deviation,
        },
        actions: report
            .decisions
            .iter()
            .map(|decision| BrepRepairActionRecord {
                stage: brep_repair_action_stage_label(decision.stage).to_string(),
                outcome: brep_repair_action_outcome_label(decision.outcome).to_string(),
                subject: decision.subject.clone(),
                reason: decision.reason.clone(),
                max_motion: decision.max_motion,
                max_deviation: decision.max_deviation,
                affected_handles: decision
                    .affected_handles
                    .iter()
                    .map(|handle| repair_action_handle_label(*handle))
                    .collect(),
                shared_fact_ids: decision
                    .shared_fact_ids
                    .iter()
                    .map(|id| id.0)
                    .collect::<Vec<_>>(),
            })
            .collect(),
    }
}

/// Stable serialized label for a repair pipeline stage.
pub fn brep_repair_action_stage_label(stage: BrepRepairPipelineStage) -> &'static str {
    match stage {
        BrepRepairPipelineStage::PolicyPreflight => "policy_preflight",
        BrepRepairPipelineStage::Classification => "classification",
        BrepRepairPipelineStage::SnapSewing => "snap_sewing",
        BrepRepairPipelineStage::EdgeCurveCompletion => "edge_curve_completion",
        BrepRepairPipelineStage::CurveTrimReconstruction => "curve_trim_reconstruction",
        BrepRepairPipelineStage::EdgeSurfaceCompatibility => "edge_surface_compatibility",
        BrepRepairPipelineStage::TrimRepair => "trim_repair",
        BrepRepairPipelineStage::LoopShellDiagnostics => "loop_shell_diagnostics",
    }
}

/// Stable serialized label for a repair action outcome.
pub fn brep_repair_action_outcome_label(outcome: BrepRepairActionOutcome) -> &'static str {
    match outcome {
        BrepRepairActionOutcome::Authored => "authored",
        BrepRepairActionOutcome::Declined => "declined",
        BrepRepairActionOutcome::RolledBack => "rolled_back",
    }
}

/// Stable serialized label for a repair-action topology handle.
pub fn repair_action_handle_label(handle: BrepRepairTopologyHandle) -> String {
    match handle {
        BrepRepairTopologyHandle::Vertex(handle) => format!("vertex:{handle:?}"),
        BrepRepairTopologyHandle::HalfEdge(handle) => format!("half_edge:{handle:?}"),
        BrepRepairTopologyHandle::Edge(handle) => format!("edge:{handle:?}"),
        BrepRepairTopologyHandle::Loop(handle) => format!("loop:{handle:?}"),
        BrepRepairTopologyHandle::Face(handle) => format!("face:{handle:?}"),
        BrepRepairTopologyHandle::Shell(handle) => format!("shell:{handle:?}"),
        BrepRepairTopologyHandle::Solid(handle) => format!("solid:{handle:?}"),
    }
}

/// One generic repair-action decision.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepRepairActionDecision {
    /// Pipeline stage that produced the decision.
    pub stage: BrepRepairPipelineStage,
    /// Machine-readable outcome class.
    pub outcome: BrepRepairActionOutcome,
    /// Stable subject label, usually a provenance label or entity handle.
    pub subject: String,
    /// Stage-specific reason text for declined or rolled-back actions.
    pub reason: Option<String>,
    /// Vertex or endpoint motion associated with this action.
    pub max_motion: f64,
    /// Sampled geometric deviation associated with this action.
    pub max_deviation: f64,
    /// Topology handles directly affected or considered by this action.
    pub affected_handles: Vec<BrepRepairTopologyHandle>,
    /// Shared repair knowledge facts referenced by this action.
    pub shared_fact_ids: Vec<BrepRepairFactId>,
}

/// Generic action outcome class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepRepairActionOutcome {
    /// The action authored topology or geometry into the repaired clone.
    Authored,
    /// The action was declined before mutation.
    Declined,
    /// The action mutated a clone but the transaction rolled back.
    RolledBack,
}
