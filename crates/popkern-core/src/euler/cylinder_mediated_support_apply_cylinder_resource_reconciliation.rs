//! Reconciliation between outside bridge cycles and cylinder-side resources.

use std::collections::{BTreeMap, BTreeSet};

use super::super::apply_bridge_cycles::SegmentBridgeCycleOperandSet;
use super::{CylinderBridgeResourcePartition, CylinderBridgeResourceSegment};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(in crate::euler::cylinder_mediated_support_apply) struct BridgeResourceReconciliation {
    outside_segments: usize,
    cylinder_segments: usize,
    exact_matches: usize,
    missing_cylinder_matches: usize,
    extra_cylinder_resources: usize,
    outside_request_sets: Vec<Vec<usize>>,
    cylinder_request_sets: Vec<Vec<usize>>,
    missing_request_sets: Vec<Vec<usize>>,
    actions: Vec<BridgeResourceConstructionAction>,
}

impl BridgeResourceReconciliation {
    pub(in crate::euler::cylinder_mediated_support_apply) fn status(
        &self,
    ) -> BridgeResourceReconciliationStatus {
        if self
            .actions
            .iter()
            .any(|action| action.kind == BridgeResourceConstructionActionKind::Unsupported)
        {
            return BridgeResourceReconciliationStatus::IncompatibleResourceSets;
        }
        if self.missing_request_sets.is_empty() && self.extra_cylinder_resources == 0 {
            return BridgeResourceReconciliationStatus::ExactResourcePairing;
        }
        if self.has_split_cylinder_resource_transition() {
            return BridgeResourceReconciliationStatus::SplitCylinderResourcesNeedMerge;
        }
        if self.extra_cylinder_resources > 0 {
            return BridgeResourceReconciliationStatus::IncompatibleResourceSets;
        }
        if self.actions.iter().any(|action| {
            action.kind == BridgeResourceConstructionActionKind::SharedOutsideTransition
        }) && self
            .actions
            .iter()
            .filter(|action| action.kind != BridgeResourceConstructionActionKind::ExactPair)
            .all(|action| {
                action.kind == BridgeResourceConstructionActionKind::SharedOutsideTransition
            })
        {
            return BridgeResourceReconciliationStatus::SharedOutsideTransitionNeedsResource;
        }
        if self.actions.iter().any(|action| {
            action.kind
                == BridgeResourceConstructionActionKind::SingleSupportCoveredBySharedCylinder
        }) && self
            .actions
            .iter()
            .filter(|action| action.kind != BridgeResourceConstructionActionKind::ExactPair)
            .all(|action| {
                action.kind
                    == BridgeResourceConstructionActionKind::SingleSupportCoveredBySharedCylinder
            })
        {
            return BridgeResourceReconciliationStatus::SingleSupportCoveredBySharedCylinderResource;
        }
        BridgeResourceReconciliationStatus::IncompatibleResourceSets
    }

    fn has_split_cylinder_resource_transition(&self) -> bool {
        duplicate_request_set_count(&self.cylinder_request_sets) > 0
            && !self.missing_request_sets.is_empty()
            && self.missing_request_sets.iter().all(|missing| {
                cylinder_request_sets_cover_missing_set(missing, &self.cylinder_request_sets)
            })
    }

    pub(in crate::euler::cylinder_mediated_support_apply) fn diagnostic_label(&self) -> String {
        format!(
            "bridge_resource_reconciliation=outside_segments={} cylinder_segments={} exact_matches={} missing_cylinder_matches={} extra_cylinder_resources={} duplicate_cylinder_request_sets={} resource_status={} construction_actions={} action_sequence=[{}] outside_request_sets=[{}] cylinder_request_sets=[{}] missing_request_sets=[{}] {}",
            self.outside_segments,
            self.cylinder_segments,
            self.exact_matches,
            self.missing_cylinder_matches,
            self.extra_cylinder_resources,
            duplicate_request_set_count(&self.cylinder_request_sets),
            self.status().diagnostic_label(),
            construction_action_summary(&self.actions),
            construction_action_sequence(&self.actions),
            request_set_summary(&self.outside_request_sets),
            request_set_summary(&self.cylinder_request_sets),
            request_set_summary(&self.missing_request_sets),
            split_resource_merge_summary(&self.missing_request_sets, &self.cylinder_request_sets)
                .diagnostic_label()
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(in crate::euler::cylinder_mediated_support_apply) enum BridgeResourceReconciliationStatus {
    ExactResourcePairing,
    SharedOutsideTransitionNeedsResource,
    SingleSupportCoveredBySharedCylinderResource,
    SplitCylinderResourcesNeedMerge,
    IncompatibleResourceSets,
}

impl BridgeResourceReconciliationStatus {
    pub(in crate::euler::cylinder_mediated_support_apply) fn diagnostic_label(
        self,
    ) -> &'static str {
        match self {
            Self::ExactResourcePairing => "exact-resource-pairing",
            Self::SharedOutsideTransitionNeedsResource => {
                "shared-outside-transition-needs-resource"
            }
            Self::SingleSupportCoveredBySharedCylinderResource => {
                "single-support-covered-by-shared-cylinder-resource"
            }
            Self::SplitCylinderResourcesNeedMerge => "split-cylinder-resources-need-merge",
            Self::IncompatibleResourceSets => "incompatible-resource-sets",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BridgeResourceConstructionAction {
    kind: BridgeResourceConstructionActionKind,
    outside_segment_index: usize,
    cylinder_segment_index: Option<usize>,
    request_indices: Vec<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BridgeResourceConstructionActionKind {
    ExactPair,
    SharedOutsideTransition,
    SingleSupportCoveredBySharedCylinder,
    Unsupported,
}

impl CylinderBridgeResourcePartition {
    pub(in crate::euler::cylinder_mediated_support_apply) fn reconcile(
        &self,
        bridge_cycles: &SegmentBridgeCycleOperandSet,
    ) -> BridgeResourceReconciliation {
        let outside_request_sets = bridge_cycles
            .request_index_sets()
            .map(<[usize]>::to_vec)
            .collect::<Vec<_>>();
        let mut matched_cylinder_segments = BTreeSet::new();
        let mut exact_matches = 0usize;
        let mut missing_request_sets = Vec::new();
        let mut actions = Vec::new();
        for (outside_segment_index, outside_request_set) in outside_request_sets.iter().enumerate()
        {
            let Some((index, _)) = self.segments.iter().enumerate().find(|(index, cylinder)| {
                !matched_cylinder_segments.contains(index)
                    && &cylinder.request_indices == outside_request_set
            }) else {
                missing_request_sets.push(outside_request_set.clone());
                actions.push(unmatched_construction_action(
                    &self.segments,
                    outside_segment_index,
                    outside_request_set,
                ));
                continue;
            };
            matched_cylinder_segments.insert(index);
            exact_matches += 1;
            actions.push(BridgeResourceConstructionAction {
                kind: BridgeResourceConstructionActionKind::ExactPair,
                outside_segment_index,
                cylinder_segment_index: Some(index),
                request_indices: outside_request_set.clone(),
            });
        }
        BridgeResourceReconciliation {
            outside_segments: outside_request_sets.len(),
            cylinder_segments: self.segments.len(),
            exact_matches,
            missing_cylinder_matches: missing_request_sets.len(),
            extra_cylinder_resources: self.segments.len() - matched_cylinder_segments.len(),
            outside_request_sets,
            cylinder_request_sets: self
                .segments
                .iter()
                .map(|segment| segment.request_indices.clone())
                .collect(),
            missing_request_sets,
            actions,
        }
    }
}

fn unmatched_construction_action(
    cylinder_segments: &[CylinderBridgeResourceSegment],
    outside_segment_index: usize,
    request_set: &[usize],
) -> BridgeResourceConstructionAction {
    if request_set.len() > 1 {
        return BridgeResourceConstructionAction {
            kind: BridgeResourceConstructionActionKind::SharedOutsideTransition,
            outside_segment_index,
            cylinder_segment_index: None,
            request_indices: request_set.to_vec(),
        };
    }
    let covering_cylinder_segment = cylinder_segments
        .iter()
        .enumerate()
        .find(|(_, cylinder_segment)| {
            cylinder_segment.request_indices.len() > 1
                && request_set
                    .iter()
                    .all(|request_index| cylinder_segment.request_indices.contains(request_index))
        })
        .map(|(index, _)| index);
    if request_set.len() == 1 && covering_cylinder_segment.is_some() {
        return BridgeResourceConstructionAction {
            kind: BridgeResourceConstructionActionKind::SingleSupportCoveredBySharedCylinder,
            outside_segment_index,
            cylinder_segment_index: covering_cylinder_segment,
            request_indices: request_set.to_vec(),
        };
    }
    BridgeResourceConstructionAction {
        kind: BridgeResourceConstructionActionKind::Unsupported,
        outside_segment_index,
        cylinder_segment_index: None,
        request_indices: request_set.to_vec(),
    }
}

fn construction_action_summary(actions: &[BridgeResourceConstructionAction]) -> String {
    let exact_pairs =
        construction_action_count(actions, BridgeResourceConstructionActionKind::ExactPair);
    let shared_outside = construction_action_count(
        actions,
        BridgeResourceConstructionActionKind::SharedOutsideTransition,
    );
    let single_covered = construction_action_count(
        actions,
        BridgeResourceConstructionActionKind::SingleSupportCoveredBySharedCylinder,
    );
    let unsupported =
        construction_action_count(actions, BridgeResourceConstructionActionKind::Unsupported);
    format!(
        "exact_pair:{} shared_outside:{} single_covered:{} unsupported:{}",
        exact_pairs, shared_outside, single_covered, unsupported
    )
}

fn construction_action_count(
    actions: &[BridgeResourceConstructionAction],
    kind: BridgeResourceConstructionActionKind,
) -> usize {
    actions.iter().filter(|action| action.kind == kind).count()
}

fn construction_action_sequence(actions: &[BridgeResourceConstructionAction]) -> String {
    actions
        .iter()
        .map(|action| {
            format!(
                "{}:{}:{}:{}",
                action.outside_segment_index,
                action.kind.diagnostic_label(),
                optional_index_label(action.cylinder_segment_index),
                request_set_label(&action.request_indices)
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

impl BridgeResourceConstructionActionKind {
    fn diagnostic_label(self) -> &'static str {
        match self {
            Self::ExactPair => "exact",
            Self::SharedOutsideTransition => "shared-outside",
            Self::SingleSupportCoveredBySharedCylinder => "single-covered",
            Self::Unsupported => "unsupported",
        }
    }
}

fn optional_index_label(index: Option<usize>) -> String {
    index.map_or_else(|| "-".to_string(), |index| index.to_string())
}

fn request_set_summary(request_sets: &[Vec<usize>]) -> String {
    request_sets
        .iter()
        .map(|set| request_set_label(set))
        .collect::<Vec<_>>()
        .join(",")
}

fn duplicate_request_set_count(request_sets: &[Vec<usize>]) -> usize {
    let mut counts = BTreeMap::<&[usize], usize>::new();
    for request_set in request_sets {
        *counts.entry(request_set).or_default() += 1;
    }
    counts.values().filter(|&&count| count > 1).count()
}

fn cylinder_request_sets_cover_missing_set(
    missing: &[usize],
    cylinder_sets: &[Vec<usize>],
) -> bool {
    let covered = cylinder_sets
        .iter()
        .filter(|set| cylinder_request_set_participates_in_missing_set(missing, set))
        .flat_map(|set| set.iter().copied())
        .collect::<BTreeSet<_>>();
    missing.iter().all(|index| covered.contains(index))
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct SplitResourceMergeSummary {
    missing_sets: usize,
    coverable_sets: usize,
    candidate_cylinder_runs: usize,
}

impl SplitResourceMergeSummary {
    fn diagnostic_label(self) -> String {
        format!(
            "split_resource_merge_candidates=missing_sets={} coverable_sets={} cylinder_runs={}",
            self.missing_sets, self.coverable_sets, self.candidate_cylinder_runs
        )
    }
}

fn split_resource_merge_summary(
    missing_sets: &[Vec<usize>],
    cylinder_sets: &[Vec<usize>],
) -> SplitResourceMergeSummary {
    let mut summary = SplitResourceMergeSummary {
        missing_sets: missing_sets.len(),
        ..SplitResourceMergeSummary::default()
    };
    for missing in missing_sets {
        let candidate_runs = cylinder_sets
            .iter()
            .filter(|set| cylinder_request_set_participates_in_missing_set(missing, set))
            .count();
        if candidate_runs > 0 && cylinder_request_sets_cover_missing_set(missing, cylinder_sets) {
            summary.coverable_sets += 1;
            summary.candidate_cylinder_runs += candidate_runs;
        }
    }
    summary
}

fn cylinder_request_set_participates_in_missing_set(missing: &[usize], set: &[usize]) -> bool {
    set.iter().all(|index| missing.contains(index))
        || missing.iter().all(|index| set.contains(index))
}

fn request_set_label(set: &[usize]) -> String {
    set.iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join("+")
}
