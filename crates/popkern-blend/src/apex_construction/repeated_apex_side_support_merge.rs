//! Support-region merge primitive for adjacent repeated-apex residual faces.

use std::collections::BTreeMap;

use popkern_core::euler::{FaceMergeSplitRequest, FaceMergeSplitResult};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::repeated_apex_side_support_merge_geometry::SideGeneratorCurveAuthoring;
use super::repeated_apex_side_support_plan::{
    BridgeSupportRegionPlan, BridgeSupportTrimPlan, PlannedTopologyOperation,
    SelectedSideSupportRegionPlan, SideSupportRegionPlanSelection,
};

#[derive(Clone, Debug)]
pub(super) struct SupportRegionMergePrimitiveDerivation {
    merge_regions: usize,
    ready_regions: usize,
    details: Vec<SupportRegionMergeRegionDerivation>,
}

#[derive(Clone, Debug)]
struct SupportRegionMergeRegionDerivation {
    support_face: FaceHandle,
    primitive: Option<SupportRegionMergePrimitivePlan>,
    reason: Option<SupportRegionMergeRefusal>,
}

#[derive(Clone, Debug)]
struct SupportRegionMergePrimitivePlan {
    support_face: FaceHandle,
    outside_face: FaceHandle,
    shared_half_edge: HalfEdgeHandle,
    continuation: Option<SupportRegionMergeContinuationPlan>,
}

#[derive(Clone, Debug)]
struct SupportRegionMergeContinuationPlan {
    retained_face: FaceHandle,
    retained_loop: LoopHandle,
    trims: Vec<PostMergeSideTrimPlan>,
    expected_delta: SupportRegionMergeDelta,
}

#[derive(Clone, Debug)]
struct PostMergeSideTrimPlan {
    fillet_index: usize,
    cylinder_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    start_half_edge: HalfEdgeHandle,
    end_half_edge: HalfEdgeHandle,
    curve_authoring: SideGeneratorCurveAuthoring,
}

#[derive(Clone, Copy, Debug)]
struct SupportRegionMergeDelta {
    edges: i64,
    faces: i64,
    loops: i64,
}

#[derive(Clone, Debug)]
enum SupportRegionMergeRefusal {
    UnsupportedTrimCount,
    MixedOutsideFaces,
    NoSharedSupportOutsideEdge,
    MultipleSharedSupportOutsideEdges,
    CoreMergeRefused(String),
}

impl SupportRegionMergePrimitiveDerivation {
    pub(super) fn build(selection: &SideSupportRegionPlanSelection, brep: &BRep) -> Self {
        match selection {
            SideSupportRegionPlanSelection::Selected(plan)
            | SideSupportRegionPlanSelection::Blocked(plan, _) => {
                Self::from_selected_plan(plan, brep)
            }
            SideSupportRegionPlanSelection::Refused(_) => Self {
                merge_regions: 0,
                ready_regions: 0,
                details: Vec::new(),
            },
        }
    }

    pub(super) fn is_ready(&self) -> bool {
        self.merge_regions > 0 && self.ready_regions == self.merge_regions
    }

    pub(super) fn exact_geometry_ready(&self) -> bool {
        self.is_ready()
            && self
                .details
                .iter()
                .all(SupportRegionMergeRegionDerivation::exact_geometry_ready)
    }

    pub(super) fn has_regions(&self) -> bool {
        self.merge_regions > 0
    }

    pub(super) fn apply_topology_and_geometry(&self, brep: &mut BRep) -> Result<(), String> {
        if !self.exact_geometry_ready() {
            return Err("support-region merge exact geometry is not ready".to_string());
        }
        for detail in &self.details {
            detail.apply_topology_and_geometry(brep)?;
        }
        Ok(())
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "support-region-merge primitive_ready={} merge_regions={} ready_regions={} details=[{}]",
            self.is_ready(),
            self.merge_regions,
            self.ready_regions,
            self.details
                .iter()
                .map(SupportRegionMergeRegionDerivation::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }

    fn from_selected_plan(plan: &SelectedSideSupportRegionPlan, brep: &BRep) -> Self {
        let details = plan
            .operations()
            .filter_map(|operation| match operation {
                PlannedTopologyOperation::BridgeSupportRegion(region) => {
                    Some(SupportRegionMergeRegionDerivation::build(region, brep))
                }
                PlannedTopologyOperation::SplitFaceWithEdge { .. }
                | PlannedTopologyOperation::OutsideOutsideSupportContinuation(_) => None,
            })
            .collect::<Vec<_>>();
        let merge_regions = details.len();
        let ready_regions = details
            .iter()
            .filter(|detail| detail.reason.is_none())
            .count();
        Self {
            merge_regions,
            ready_regions,
            details,
        }
    }
}

impl SupportRegionMergeRegionDerivation {
    fn build(region: &BridgeSupportRegionPlan, brep: &BRep) -> Self {
        if region.trims().len() != 2 {
            return Self::refused(
                region.support_face(),
                SupportRegionMergeRefusal::UnsupportedTrimCount,
            );
        }
        let Some(outside_face) = single_outside_face(region.trims()) else {
            return Self::refused(
                region.support_face(),
                SupportRegionMergeRefusal::MixedOutsideFaces,
            );
        };
        let shared_half_edges = shared_support_half_edges(region, brep);
        let shared_half_edge = match shared_half_edges.as_slice() {
            [] => {
                return Self::refused(
                    region.support_face(),
                    SupportRegionMergeRefusal::NoSharedSupportOutsideEdge,
                );
            }
            [shared] => *shared,
            _ => {
                return Self::refused(
                    region.support_face(),
                    SupportRegionMergeRefusal::MultipleSharedSupportOutsideEdges,
                );
            }
        };
        let split_requests = region
            .trims()
            .iter()
            .map(|trim| FaceMergeSplitRequest {
                start: trim.start(),
                end: trim.end(),
            })
            .collect::<Vec<_>>();
        let mut candidate = brep.clone();
        let topology_result =
            match candidate.merge_faces_across_edge_then_split(shared_half_edge, &split_requests) {
                Ok(result) => result,
                Err(err) => {
                    return Self::refused(
                        region.support_face(),
                        SupportRegionMergeRefusal::CoreMergeRefused(err.to_string()),
                    );
                }
            };
        let continuation =
            SupportRegionMergeContinuationPlan::derive(&candidate, region, &topology_result);
        Self {
            support_face: region.support_face(),
            primitive: Some(SupportRegionMergePrimitivePlan {
                support_face: region.support_face(),
                outside_face,
                shared_half_edge,
                continuation: Some(continuation),
            }),
            reason: None,
        }
    }

    fn refused(support_face: FaceHandle, reason: SupportRegionMergeRefusal) -> Self {
        Self {
            support_face,
            primitive: None,
            reason: Some(reason),
        }
    }

    fn diagnostic_label(&self) -> String {
        match (&self.primitive, &self.reason) {
            (Some(primitive), None) => primitive.diagnostic_label(),
            (None, Some(reason)) => {
                format!(
                    "support={:?} {}",
                    self.support_face,
                    reason.diagnostic_label()
                )
            }
            _ => format!(
                "support={:?} inconsistent-support-region-merge-derivation",
                self.support_face
            ),
        }
    }

    fn exact_geometry_ready(&self) -> bool {
        self.primitive
            .as_ref()
            .is_some_and(SupportRegionMergePrimitivePlan::exact_geometry_ready)
    }

    fn apply_topology_and_geometry(&self, brep: &mut BRep) -> Result<(), String> {
        self.primitive
            .as_ref()
            .ok_or_else(|| "support-region merge primitive is incomplete".to_string())?
            .apply_topology_and_geometry(brep)
    }
}

impl SupportRegionMergePrimitivePlan {
    fn diagnostic_label(&self) -> String {
        format!(
            "support={:?} outside={:?} shared_half_edge={:?} core_merge_probe=ready delta=dV=0 dE=-1 dF=-1 dL=-1 continuation={}",
            self.support_face,
            self.outside_face,
            self.shared_half_edge,
            self.continuation
                .as_ref()
                .map(SupportRegionMergeContinuationPlan::diagnostic_label)
                .unwrap_or_else(|| "missing-post-merge-continuation".to_string())
        )
    }

    fn exact_geometry_ready(&self) -> bool {
        self.continuation
            .as_ref()
            .is_some_and(SupportRegionMergeContinuationPlan::exact_geometry_ready)
    }

    fn apply_topology_and_geometry(&self, brep: &mut BRep) -> Result<(), String> {
        let continuation = self
            .continuation
            .as_ref()
            .ok_or_else(|| "missing support-region merge continuation".to_string())?;
        let split_requests = continuation
            .trims
            .iter()
            .map(|trim| FaceMergeSplitRequest {
                start: trim.start,
                end: trim.end,
            })
            .collect::<Vec<_>>();
        let result = brep
            .merge_faces_across_edge_then_split(self.shared_half_edge, &split_requests)
            .map_err(|err| err.to_string())?;
        continuation.apply_geometry(brep, &result)
    }
}

impl SupportRegionMergeContinuationPlan {
    fn derive(
        applied: &BRep,
        region: &BridgeSupportRegionPlan,
        topology_result: &FaceMergeSplitResult,
    ) -> Self {
        let retained_face = topology_result.merge.retained_face;
        let retained_loop = applied.faces[retained_face].outer_loop;
        let mut trims = Vec::new();
        for (trim, split) in region.trims().iter().zip(&topology_result.splits) {
            trims.push(PostMergeSideTrimPlan {
                fillet_index: trim.fillet_index(),
                cylinder_face: trim.cylinder_face(),
                start: trim.start(),
                end: trim.end(),
                start_half_edge: split.half_edge,
                end_half_edge: split.half_edge_twin,
                curve_authoring: SideGeneratorCurveAuthoring::build(
                    applied,
                    trim.cylinder_face(),
                    retained_face,
                    trim.start(),
                    trim.end(),
                ),
            });
        }
        Self {
            retained_face,
            retained_loop,
            expected_delta: SupportRegionMergeDelta::for_trim_count(trims.len()),
            trims,
        }
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "post_merge_side_trims_ready=true retained={:?}/{:?} trims=[{}] combined_delta={}",
            self.retained_face,
            self.retained_loop,
            self.trims
                .iter()
                .map(PostMergeSideTrimPlan::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.expected_delta.diagnostic_label()
        )
    }

    fn exact_geometry_ready(&self) -> bool {
        self.trims
            .iter()
            .all(PostMergeSideTrimPlan::exact_geometry_ready)
    }

    fn apply_geometry(&self, brep: &mut BRep, result: &FaceMergeSplitResult) -> Result<(), String> {
        if result.splits.len() != self.trims.len() {
            return Err(format!(
                "support-region merge split count mismatch: expected {} got {}",
                self.trims.len(),
                result.splits.len()
            ));
        }
        for (trim, split) in self.trims.iter().zip(&result.splits) {
            trim.apply_geometry(brep, split.half_edge, split.half_edge_twin)?;
        }
        Ok(())
    }
}

impl PostMergeSideTrimPlan {
    fn diagnostic_label(&self) -> String {
        format!(
            "FilletEnd[{}] cylinder={:?} {:?}->{:?} half_edges={:?}->{:?} {}",
            self.fillet_index,
            self.cylinder_face,
            self.start,
            self.end,
            self.start_half_edge,
            self.end_half_edge,
            self.curve_authoring.diagnostic_label()
        )
    }

    fn exact_geometry_ready(&self) -> bool {
        self.curve_authoring.is_ready()
    }

    fn apply_geometry(
        &self,
        brep: &mut BRep,
        half_edge: HalfEdgeHandle,
        half_edge_twin: HalfEdgeHandle,
    ) -> Result<(), String> {
        self.curve_authoring
            .apply_to_support_split(brep, half_edge, half_edge_twin)
    }
}

impl SupportRegionMergeDelta {
    fn for_trim_count(trim_count: usize) -> Self {
        Self {
            edges: -1 + trim_count as i64,
            faces: -1 + trim_count as i64,
            loops: -1 + trim_count as i64,
        }
    }

    fn diagnostic_label(&self) -> String {
        format!("dV=0 dE={} dF={} dL={}", self.edges, self.faces, self.loops)
    }
}

impl SupportRegionMergeRefusal {
    fn diagnostic_label(&self) -> String {
        match self {
            SupportRegionMergeRefusal::UnsupportedTrimCount => "unsupported-trim-count".to_string(),
            SupportRegionMergeRefusal::MixedOutsideFaces => "mixed-outside-faces".to_string(),
            SupportRegionMergeRefusal::NoSharedSupportOutsideEdge => {
                "no-shared-support-outside-edge".to_string()
            }
            SupportRegionMergeRefusal::MultipleSharedSupportOutsideEdges => {
                "multiple-shared-support-outside-edges".to_string()
            }
            SupportRegionMergeRefusal::CoreMergeRefused(err) => {
                format!("core-merge-refused={err}")
            }
        }
    }
}

fn single_outside_face(trims: &[BridgeSupportTrimPlan]) -> Option<FaceHandle> {
    let first = trims.first()?.outside_face();
    trims
        .iter()
        .all(|trim| trim.outside_face() == first)
        .then_some(first)
}

fn shared_support_half_edges(region: &BridgeSupportRegionPlan, brep: &BRep) -> Vec<HalfEdgeHandle> {
    let Some(first_trim) = region.trims().first() else {
        return Vec::new();
    };
    let support_edges = loop_edges(brep, first_trim.support_half_edge());
    let mut shared = Vec::new();
    for trim in region.trims() {
        for (edge, _) in loop_edges(brep, trim.outside_half_edge()) {
            if let Some(support_half_edge) = support_edges.get(&edge) {
                shared.push(*support_half_edge);
            }
        }
    }
    shared.sort();
    shared.dedup();
    shared
}

fn loop_edges(
    brep: &BRep,
    start: HalfEdgeHandle,
) -> BTreeMap<popkern_core::topology::EdgeHandle, HalfEdgeHandle> {
    let mut edges = BTreeMap::new();
    let Some(first) = brep.half_edges.get(start) else {
        return edges;
    };
    let loop_handle = first.loop_;
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let Some(half_edge) = brep.half_edges.get(cursor) else {
            return BTreeMap::new();
        };
        if half_edge.loop_ != loop_handle {
            return BTreeMap::new();
        }
        edges.insert(half_edge.edge, cursor);
        cursor = half_edge.next;
        if cursor == start {
            return edges;
        }
    }
    BTreeMap::new()
}
