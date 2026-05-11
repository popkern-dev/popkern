//! Cylinder-side resource partitioning for mediated support continuation.

use std::collections::{BTreeMap, BTreeSet};

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle, LoopHandle};

use super::CylinderMediatedSupportContinuationBatchPlan;
use super::apply_bridge_cycles::SegmentBridgeCycleOperandSet;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct BridgeResourceReconciliation {
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
    pub(super) fn status(&self) -> BridgeResourceReconciliationStatus {
        if self
            .actions
            .iter()
            .any(|action| action.kind == BridgeResourceConstructionActionKind::Unsupported)
            || self.extra_cylinder_resources > 0
        {
            return BridgeResourceReconciliationStatus::IncompatibleResourceSets;
        }
        if self.missing_request_sets.is_empty() && self.extra_cylinder_resources == 0 {
            return BridgeResourceReconciliationStatus::ExactResourcePairing;
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

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "bridge_resource_reconciliation=outside_segments={} cylinder_segments={} exact_matches={} missing_cylinder_matches={} extra_cylinder_resources={} duplicate_cylinder_request_sets={} resource_status={} construction_actions={} action_sequence=[{}] outside_request_sets=[{}] cylinder_request_sets=[{}] missing_request_sets=[{}]",
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
            request_set_summary(&self.missing_request_sets)
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum BridgeResourceReconciliationStatus {
    ExactResourcePairing,
    SharedOutsideTransitionNeedsResource,
    SingleSupportCoveredBySharedCylinderResource,
    IncompatibleResourceSets,
}

impl BridgeResourceReconciliationStatus {
    pub(super) fn diagnostic_label(self) -> &'static str {
        match self {
            Self::ExactResourcePairing => "exact-resource-pairing",
            Self::SharedOutsideTransitionNeedsResource => {
                "shared-outside-transition-needs-resource"
            }
            Self::SingleSupportCoveredBySharedCylinderResource => {
                "single-support-covered-by-shared-cylinder-resource"
            }
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct CylinderBridgeResourcePartition {
    segments: Vec<CylinderBridgeResourceSegment>,
}

impl CylinderBridgeResourcePartition {
    pub(super) fn from_remapped_segments(
        segments: &[CylinderBridgeResourceSegment],
    ) -> CylinderBridgeResourcePartition {
        CylinderBridgeResourcePartition {
            segments: segments.to_vec(),
        }
    }

    pub(super) fn segments(&self) -> impl Iterator<Item = &CylinderBridgeResourceSegment> {
        self.segments.iter()
    }

    pub(super) fn reconcile(
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

    pub(super) fn diagnostic_label(&self) -> String {
        let shared_segments = self
            .segments
            .iter()
            .filter(|segment| segment.request_indices.len() > 1)
            .count();
        let max_request_owners = self
            .segments
            .iter()
            .map(|segment| segment.request_indices.len())
            .max()
            .unwrap_or_default();
        let half_edges = self
            .segments
            .iter()
            .map(|segment| segment.half_edges.len())
            .sum::<usize>();
        format!(
            "cylinder_resource_partition=segments={} half_edges={} shared_segments={} max_request_owners={}",
            self.segments.len(),
            half_edges,
            shared_segments,
            max_request_owners
        )
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct CylinderBridgeResourceSegment {
    pub(super) half_edges: Vec<HalfEdgeHandle>,
    pub(super) edges: Vec<EdgeHandle>,
    pub(super) request_indices: Vec<usize>,
}

pub(super) fn derive_cylinder_bridge_resource_partition(
    brep: &BRep,
    plan: &CylinderMediatedSupportContinuationBatchPlan,
) -> Result<CylinderBridgeResourcePartition, CoreError> {
    let mut selected_by_loop = BTreeMap::<_, BTreeMap<_, BTreeSet<usize>>>::new();
    for (request_index, request_plan) in plan.plans.iter().enumerate() {
        for &half_edge in &request_plan.request.cylinder_interval {
            let loop_ = brep
                .half_edges
                .get(half_edge)
                .map(|data| data.loop_)
                .ok_or_else(|| {
                    CoreError::InvalidHandle(format!(
                        "cylinder resource partition references stale half-edge {half_edge:?}"
                    ))
                })?;
            selected_by_loop
                .entry(loop_)
                .or_default()
                .entry(half_edge)
                .or_default()
                .insert(request_index);
        }
    }
    let mut segments = Vec::new();
    for (loop_, membership) in selected_by_loop {
        segments.extend(partition_cylinder_loop(brep, loop_, &membership)?);
    }
    Ok(CylinderBridgeResourcePartition { segments })
}

fn partition_cylinder_loop(
    brep: &BRep,
    loop_: LoopHandle,
    membership: &BTreeMap<HalfEdgeHandle, BTreeSet<usize>>,
) -> Result<Vec<CylinderBridgeResourceSegment>, CoreError> {
    let start = brep
        .loops
        .get(loop_)
        .map(|loop_data| loop_data.half_edge)
        .ok_or_else(|| {
            CoreError::InvalidHandle(format!("cylinder resource loop {loop_:?} is stale"))
        })?;
    let mut segments = Vec::new();
    let mut current_half_edges = Vec::new();
    let mut current_edges = Vec::new();
    let mut current_request_indices = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let data = brep.half_edges.get(cursor).ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "cylinder resource partition references stale loop half-edge {cursor:?}"
            ))
        })?;
        if data.loop_ != loop_ {
            return Err(CoreError::PreconditionFailed(format!(
                "cylinder resource partition traversal left loop {loop_:?} at {cursor:?}"
            )));
        }
        match membership.get(&cursor) {
            Some(request_indices) => {
                let request_indices = request_indices.iter().copied().collect::<Vec<_>>();
                if current_half_edges.is_empty() || current_request_indices == request_indices {
                    current_half_edges.push(cursor);
                    current_edges.push(data.edge);
                    current_request_indices = request_indices;
                } else {
                    push_cylinder_resource_segment(
                        &mut segments,
                        &mut current_half_edges,
                        &mut current_edges,
                        &mut current_request_indices,
                    );
                    current_half_edges.push(cursor);
                    current_edges.push(data.edge);
                    current_request_indices = request_indices;
                }
            }
            None => push_cylinder_resource_segment(
                &mut segments,
                &mut current_half_edges,
                &mut current_edges,
                &mut current_request_indices,
            ),
        }
        cursor = data.next;
        if cursor == start {
            break;
        }
    }
    if cursor != start {
        return Err(CoreError::PreconditionFailed(format!(
            "cylinder resource partition could not close loop {loop_:?}"
        )));
    }
    push_cylinder_resource_segment(
        &mut segments,
        &mut current_half_edges,
        &mut current_edges,
        &mut current_request_indices,
    );
    Ok(segments)
}

fn push_cylinder_resource_segment(
    segments: &mut Vec<CylinderBridgeResourceSegment>,
    current_half_edges: &mut Vec<HalfEdgeHandle>,
    current_edges: &mut Vec<EdgeHandle>,
    current_request_indices: &mut Vec<usize>,
) {
    if current_half_edges.is_empty() {
        return;
    }
    segments.push(CylinderBridgeResourceSegment {
        half_edges: std::mem::take(current_half_edges),
        edges: std::mem::take(current_edges),
        request_indices: std::mem::take(current_request_indices),
    });
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

fn request_set_label(set: &[usize]) -> String {
    set.iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join("+")
}
