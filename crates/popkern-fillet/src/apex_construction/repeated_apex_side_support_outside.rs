//! Outside/outside support-continuation primitive derivation.

use super::repeated_apex_side_support_plan::{
    OutsideOutsideSupportContinuationPlan, PlannedTopologyOperation, SelectedSideSupportRegionPlan,
    SideSupportRegionPlanSelection,
};

#[derive(Clone, Debug)]
pub(super) struct OutsideOutsideSupportPrimitiveDerivation {
    regions: usize,
    ready_regions: usize,
    details: Vec<OutsideOutsideSupportRegionDerivation>,
}

#[derive(Clone, Debug)]
struct OutsideOutsideSupportRegionDerivation {
    region: OutsideOutsideSupportContinuationPlan,
    reason: OutsideOutsideSupportRefusal,
}

#[derive(Clone, Debug)]
enum OutsideOutsideSupportRefusal {
    MissingTopologyOperator,
}

impl OutsideOutsideSupportPrimitiveDerivation {
    pub(super) fn build(selection: &SideSupportRegionPlanSelection) -> Self {
        match selection {
            SideSupportRegionPlanSelection::Selected(plan)
            | SideSupportRegionPlanSelection::Blocked(plan, _) => Self::from_selected_plan(plan),
            SideSupportRegionPlanSelection::Refused(_) => Self {
                regions: 0,
                ready_regions: 0,
                details: Vec::new(),
            },
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "outside-outside-support primitive_ready={} regions={} ready_regions={} details=[{}]",
            self.is_ready(),
            self.regions,
            self.ready_regions,
            self.details
                .iter()
                .map(OutsideOutsideSupportRegionDerivation::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }

    fn is_ready(&self) -> bool {
        self.regions > 0 && self.ready_regions == self.regions
    }

    fn from_selected_plan(plan: &SelectedSideSupportRegionPlan) -> Self {
        let details = plan
            .operations()
            .filter_map(|operation| match operation {
                PlannedTopologyOperation::OutsideOutsideSupportContinuation(region) => {
                    Some(OutsideOutsideSupportRegionDerivation {
                        region: region.clone(),
                        reason: OutsideOutsideSupportRefusal::MissingTopologyOperator,
                    })
                }
                PlannedTopologyOperation::SplitFaceWithEdge { .. }
                | PlannedTopologyOperation::BridgeSupportRegion(_) => None,
            })
            .collect::<Vec<_>>();
        Self {
            regions: details.len(),
            ready_regions: 0,
            details,
        }
    }
}

impl OutsideOutsideSupportRegionDerivation {
    fn diagnostic_label(&self) -> String {
        format!(
            "support={:?} trims={} {}",
            self.region.support_face(),
            self.region.trims().len(),
            self.reason.diagnostic_label()
        )
    }
}

impl OutsideOutsideSupportRefusal {
    fn diagnostic_label(&self) -> &'static str {
        match self {
            OutsideOutsideSupportRefusal::MissingTopologyOperator => "missing-topology-operator",
        }
    }
}
