//! Trim-subcomplex planning diagnostics for repeated-apex partition rewrites.

use popkern_core::topology::{FaceHandle, VertexHandle};

use super::error::ApexBuildError;

#[path = "repeated_apex_trim_subcomplex_topology.rs"]
mod trim_subcomplex_topology;

pub(super) use trim_subcomplex_topology::{
    TopologyBalanceRequirement, TopologyCountDelta, topology_count_summary,
};

pub(super) struct TrimRewriteSubcomplexPlan {
    edge_plans: Vec<TrimRewriteEdgePlan>,
    expected_edge_only_delta: TopologyCountDelta,
    support_groups: Vec<TrimRewriteSupportGroup>,
}

pub(super) struct TrimRewriteSubcomplexBuilder {
    edge_plans: Vec<TrimRewriteEdgePlan>,
    expected_edge_only_delta: TopologyCountDelta,
}

#[derive(Clone, Copy)]
struct TrimRewriteEdgePlan {
    fillet_index: usize,
    owner_pair: TrimOwnerPair,
    construction: TrimEdgeConstruction,
    start: VertexHandle,
    end: VertexHandle,
}

pub(super) struct TrimRewriteSupportGroup {
    support_face: FaceHandle,
    generated_trim_count: usize,
    cylinder_faces: Vec<FaceHandle>,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct TrimRewriteIntentFacts {
    pub(super) fillet_index: usize,
    pub(super) is_new_generator_trim: bool,
    pub(super) cylinder_face: FaceHandle,
    pub(super) support_face: Option<FaceHandle>,
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
}

#[derive(Clone, Copy, Debug)]
enum TrimOwnerPair {
    CylinderCap {
        cylinder_face: FaceHandle,
        cap_face: FaceHandle,
    },
    CylinderSupport {
        cylinder_face: FaceHandle,
        support_face: FaceHandle,
    },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TrimEdgeConstruction {
    ExistingSeam,
    NewExactGeneratorTrim,
    NewExactSupportCylinderTrim,
}

fn trim_rewrite_support_groups(edge_plans: &[TrimRewriteEdgePlan]) -> Vec<TrimRewriteSupportGroup> {
    let mut groups: std::collections::BTreeMap<
        FaceHandle,
        (usize, std::collections::BTreeSet<FaceHandle>),
    > = std::collections::BTreeMap::new();
    for edge_plan in edge_plans {
        if !edge_plan.construction.is_generated() {
            continue;
        }
        let TrimOwnerPair::CylinderSupport {
            cylinder_face,
            support_face,
        } = edge_plan.owner_pair
        else {
            continue;
        };
        let (count, cylinder_faces) = groups.entry(support_face).or_default();
        *count += 1;
        cylinder_faces.insert(cylinder_face);
    }
    groups
        .into_iter()
        .map(
            |(support_face, (generated_trim_count, cylinder_faces))| TrimRewriteSupportGroup {
                support_face,
                generated_trim_count,
                cylinder_faces: cylinder_faces.into_iter().collect(),
            },
        )
        .collect()
}

impl TrimRewriteSubcomplexBuilder {
    pub(super) fn new() -> Self {
        Self {
            edge_plans: Vec::new(),
            expected_edge_only_delta: TopologyCountDelta::zero(),
        }
    }

    pub(super) fn add_topology_delta(&mut self, delta: TopologyCountDelta) {
        self.expected_edge_only_delta.add_assign(delta);
    }

    pub(super) fn push_existing_seam(
        &mut self,
        fillet_index: usize,
        cylinder_face: FaceHandle,
        cap_face: FaceHandle,
        start: VertexHandle,
        end: VertexHandle,
    ) {
        self.edge_plans.push(TrimRewriteEdgePlan::existing_seam(
            fillet_index,
            cylinder_face,
            cap_face,
            start,
            end,
        ));
    }

    pub(super) fn push_new_generator_trim(
        &mut self,
        fillet_index: usize,
        cylinder_face: FaceHandle,
        support_face: FaceHandle,
        start: VertexHandle,
        end: VertexHandle,
    ) {
        self.edge_plans
            .push(TrimRewriteEdgePlan::new_generator_trim(
                fillet_index,
                cylinder_face,
                support_face,
                start,
                end,
            ));
        self.add_topology_delta(TopologyCountDelta::from_components(0, 1, 0, 0, 0));
    }

    pub(super) fn push_new_support_cylinder_trim(
        &mut self,
        fillet_index: usize,
        cylinder_face: FaceHandle,
        support_face: FaceHandle,
        start: VertexHandle,
        end: VertexHandle,
    ) {
        self.edge_plans
            .push(TrimRewriteEdgePlan::new_support_cylinder_trim(
                fillet_index,
                cylinder_face,
                support_face,
                start,
                end,
            ));
    }

    pub(super) fn finish(self) -> TrimRewriteSubcomplexPlan {
        TrimRewriteSubcomplexPlan::from_edge_plans(self.edge_plans, self.expected_edge_only_delta)
    }
}

impl TrimRewriteSubcomplexPlan {
    fn from_edge_plans(
        edge_plans: Vec<TrimRewriteEdgePlan>,
        expected_edge_only_delta: TopologyCountDelta,
    ) -> Self {
        Self {
            support_groups: trim_rewrite_support_groups(&edge_plans),
            edge_plans,
            expected_edge_only_delta,
        }
    }

    fn edge_plans(&self) -> &[TrimRewriteEdgePlan] {
        &self.edge_plans
    }

    pub(super) fn support_groups(&self) -> &[TrimRewriteSupportGroup] {
        &self.support_groups
    }

    pub(super) fn seed_faces(&self) -> Vec<FaceHandle> {
        let mut seed_faces = std::collections::BTreeSet::new();
        for group in self.support_groups() {
            seed_faces.insert(group.support_face());
            seed_faces.extend(group.cylinder_faces().iter().copied());
        }
        for edge_plan in self.edge_plans() {
            seed_faces.extend(edge_plan.owner_faces());
        }
        seed_faces.into_iter().collect()
    }

    pub(super) fn trim_intents(&self) -> Vec<TrimRewriteIntentFacts> {
        self.edge_plans()
            .iter()
            .map(TrimRewriteEdgePlan::intent_facts)
            .collect()
    }

    pub(super) fn planned_generator_support_face(
        &self,
        cylinder_face: FaceHandle,
        fillet_index: usize,
        start: VertexHandle,
        end: VertexHandle,
        label: &'static str,
    ) -> Result<FaceHandle, ApexBuildError> {
        let candidates = self
            .edge_plans()
            .iter()
            .filter(|edge_plan| {
                edge_plan.fillet_index() == fillet_index
                    && edge_plan.start() == start
                    && edge_plan.end() == end
                    && edge_plan.is_new_generator_trim()
                    && edge_plan.cylinder_face() == cylinder_face
            })
            .filter_map(TrimRewriteEdgePlan::support_face)
            .collect::<Vec<_>>();
        if candidates.len() == 1 {
            return Ok(candidates[0]);
        }
        Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot apply \
                 generator trim {label} {start:?}->{end:?}: expected one planned \
                 cylinder/support owner in trim subcomplex for cylinder {cylinder_face:?}, got \
                 {:?}. {}",
                candidates,
                self.diagnostic_label()
            ),
        })
    }

    pub(super) fn diagnostic_label(&self) -> String {
        let generated_edges = self
            .edge_plans
            .iter()
            .filter(|edge| edge.construction.is_generated())
            .count();
        let existing_edges = self.edge_plans.len() - generated_edges;
        let support_groups = self
            .support_groups
            .iter()
            .map(TrimRewriteSupportGroup::diagnostic_label)
            .collect::<Vec<_>>()
            .join("; ");
        let by_owner = self
            .edge_plans
            .iter()
            .map(TrimRewriteEdgePlan::diagnostic_label)
            .collect::<Vec<_>>()
            .join("; ");
        format!(
            "{} planned trim side(s): {} generated, {} existing, edge-only delta {}; support \
             groups [{}]; [{}]",
            self.edge_plans.len(),
            generated_edges,
            existing_edges,
            self.expected_edge_only_delta.diagnostic_label(),
            support_groups,
            by_owner
        )
    }
}

impl TrimRewriteSupportGroup {
    pub(super) fn support_face(&self) -> FaceHandle {
        self.support_face
    }

    pub(super) fn generated_trim_count(&self) -> usize {
        self.generated_trim_count
    }

    pub(super) fn cylinder_faces(&self) -> &[FaceHandle] {
        &self.cylinder_faces
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "{:?}: {} generated trim(s) from cylinder face(s) {:?}",
            self.support_face, self.generated_trim_count, self.cylinder_faces
        )
    }
}

impl TrimRewriteEdgePlan {
    fn existing_seam(
        fillet_index: usize,
        cylinder_face: FaceHandle,
        cap_face: FaceHandle,
        start: VertexHandle,
        end: VertexHandle,
    ) -> Self {
        Self {
            fillet_index,
            owner_pair: TrimOwnerPair::CylinderCap {
                cylinder_face,
                cap_face,
            },
            construction: TrimEdgeConstruction::ExistingSeam,
            start,
            end,
        }
    }

    fn new_generator_trim(
        fillet_index: usize,
        cylinder_face: FaceHandle,
        support_face: FaceHandle,
        start: VertexHandle,
        end: VertexHandle,
    ) -> Self {
        Self {
            fillet_index,
            owner_pair: TrimOwnerPair::CylinderSupport {
                cylinder_face,
                support_face,
            },
            construction: TrimEdgeConstruction::NewExactGeneratorTrim,
            start,
            end,
        }
    }

    fn new_support_cylinder_trim(
        fillet_index: usize,
        cylinder_face: FaceHandle,
        support_face: FaceHandle,
        start: VertexHandle,
        end: VertexHandle,
    ) -> Self {
        Self {
            fillet_index,
            owner_pair: TrimOwnerPair::CylinderSupport {
                cylinder_face,
                support_face,
            },
            construction: TrimEdgeConstruction::NewExactSupportCylinderTrim,
            start,
            end,
        }
    }

    fn fillet_index(&self) -> usize {
        self.fillet_index
    }

    fn start(&self) -> VertexHandle {
        self.start
    }

    fn end(&self) -> VertexHandle {
        self.end
    }

    fn owner_faces(&self) -> [FaceHandle; 2] {
        self.owner_pair.owner_faces()
    }

    fn cylinder_face(&self) -> FaceHandle {
        self.owner_pair.cylinder_face()
    }

    fn support_face(&self) -> Option<FaceHandle> {
        self.owner_pair.support_face()
    }

    fn is_new_generator_trim(&self) -> bool {
        self.construction == TrimEdgeConstruction::NewExactGeneratorTrim
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "FilletEnd[{}] {} {} {:?}->{:?}",
            self.fillet_index,
            self.construction.diagnostic_label(),
            self.owner_pair.diagnostic_label(),
            self.start,
            self.end
        )
    }

    fn intent_facts(&self) -> TrimRewriteIntentFacts {
        TrimRewriteIntentFacts {
            fillet_index: self.fillet_index(),
            is_new_generator_trim: self.is_new_generator_trim(),
            cylinder_face: self.cylinder_face(),
            support_face: self.support_face(),
            start: self.start(),
            end: self.end(),
        }
    }
}

impl TrimOwnerPair {
    fn owner_faces(&self) -> [FaceHandle; 2] {
        match *self {
            TrimOwnerPair::CylinderCap {
                cylinder_face,
                cap_face,
            } => [cylinder_face, cap_face],
            TrimOwnerPair::CylinderSupport {
                cylinder_face,
                support_face,
            } => [cylinder_face, support_face],
        }
    }

    fn cylinder_face(&self) -> FaceHandle {
        match *self {
            TrimOwnerPair::CylinderCap { cylinder_face, .. }
            | TrimOwnerPair::CylinderSupport { cylinder_face, .. } => cylinder_face,
        }
    }

    fn support_face(&self) -> Option<FaceHandle> {
        match *self {
            TrimOwnerPair::CylinderCap { .. } => None,
            TrimOwnerPair::CylinderSupport { support_face, .. } => Some(support_face),
        }
    }

    fn diagnostic_label(&self) -> String {
        match self {
            TrimOwnerPair::CylinderCap {
                cylinder_face,
                cap_face,
            } => format!("owner=cylinder/cap({cylinder_face:?}/{cap_face:?})"),
            TrimOwnerPair::CylinderSupport {
                cylinder_face,
                support_face,
            } => format!("owner=cylinder/support({cylinder_face:?}/{support_face:?})"),
        }
    }
}

impl TrimEdgeConstruction {
    fn is_generated(&self) -> bool {
        !matches!(self, TrimEdgeConstruction::ExistingSeam)
    }

    fn diagnostic_label(&self) -> &'static str {
        match self {
            TrimEdgeConstruction::ExistingSeam => "existing seam",
            TrimEdgeConstruction::NewExactGeneratorTrim => "new exact generator",
            TrimEdgeConstruction::NewExactSupportCylinderTrim => "new exact support/cylinder trim",
        }
    }
}

#[cfg(test)]
#[path = "repeated_apex_trim_subcomplex_tests.rs"]
mod tests;
