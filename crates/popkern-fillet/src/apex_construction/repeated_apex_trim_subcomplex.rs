//! Trim-subcomplex planning diagnostics for repeated-apex partition rewrites.

use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::error::ApexBuildError;

#[derive(Clone, Copy)]
pub(super) struct TopologyCountSummary {
    vertices: i64,
    edges: i64,
    faces: i64,
    loops: i64,
    shells: i64,
}

impl TopologyCountSummary {
    pub(super) fn any_genus_lhs(self) -> i64 {
        self.vertices - self.edges + 2 * self.faces - self.loops
    }

    pub(super) fn genus_zero_rhs(self) -> i64 {
        2 * self.shells
    }

    pub(super) fn deficit_to_genus_zero(self) -> i64 {
        self.genus_zero_rhs() - self.any_genus_lhs()
    }

    pub(super) fn trace_label(self, label: &'static str) -> String {
        format!(
            "[a2-trace §14.9/topology-counts/{label}] V={} E={} F={} L={} S={} lhs={}",
            self.vertices,
            self.edges,
            self.faces,
            self.loops,
            self.shells,
            self.any_genus_lhs()
        )
    }

    pub(super) fn diagnostic_label(self) -> String {
        format!(
            "V={} E={} F={} L={} S={} lhs={} rhs={} deficit={}",
            self.vertices,
            self.edges,
            self.faces,
            self.loops,
            self.shells,
            self.any_genus_lhs(),
            self.genus_zero_rhs(),
            self.deficit_to_genus_zero()
        )
    }
}

pub(super) fn topology_count_summary(brep: &BRep) -> TopologyCountSummary {
    TopologyCountSummary {
        vertices: brep.vertices.len() as i64,
        edges: brep.edges.len() as i64,
        faces: brep.faces.len() as i64,
        loops: brep.loops.len() as i64,
        shells: brep.shells.len() as i64,
    }
}

pub(super) struct TopologyCountDelta {
    vertices: i64,
    edges: i64,
    faces: i64,
    loops: i64,
    shells: i64,
    lhs: i64,
    rhs: i64,
    deficit: i64,
}

impl TopologyCountDelta {
    pub(super) fn zero() -> Self {
        Self::from_components(0, 0, 0, 0, 0)
    }

    pub(super) fn from_components(
        vertices: i64,
        edges: i64,
        faces: i64,
        loops: i64,
        shells: i64,
    ) -> Self {
        let lhs = vertices - edges + 2 * faces - loops;
        let rhs = 2 * shells;
        Self {
            vertices,
            edges,
            faces,
            loops,
            shells,
            lhs,
            rhs,
            deficit: rhs - lhs,
        }
    }

    pub(super) fn between(before: TopologyCountSummary, after: TopologyCountSummary) -> Self {
        Self {
            vertices: after.vertices - before.vertices,
            edges: after.edges - before.edges,
            faces: after.faces - before.faces,
            loops: after.loops - before.loops,
            shells: after.shells - before.shells,
            lhs: after.any_genus_lhs() - before.any_genus_lhs(),
            rhs: after.genus_zero_rhs() - before.genus_zero_rhs(),
            deficit: after.deficit_to_genus_zero() - before.deficit_to_genus_zero(),
        }
    }

    pub(super) fn add_assign(&mut self, other: Self) {
        self.vertices += other.vertices;
        self.edges += other.edges;
        self.faces += other.faces;
        self.loops += other.loops;
        self.shells += other.shells;
        self.lhs += other.lhs;
        self.rhs += other.rhs;
        self.deficit += other.deficit;
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "dV={} dE={} dF={} dL={} dS={} dLhs={} dRhs={} dDeficit={}",
            self.vertices,
            self.edges,
            self.faces,
            self.loops,
            self.shells,
            self.lhs,
            self.rhs,
            self.deficit
        )
    }
}

pub(super) struct TopologyBalanceRequirement {
    required_lhs_increase: i64,
}

impl TopologyBalanceRequirement {
    pub(super) fn from_deficit_increase(deficit_increase: i64) -> Self {
        Self {
            required_lhs_increase: deficit_increase.max(0),
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        if self.required_lhs_increase == 0 {
            return "no additional Euler compensation required at this stage".to_string();
        }
        format!(
            "needs +{} to V-E+2F-L before commit; with unchanged V/S, that is at least {} \
             additional face-loop unit(s) or an equivalent validated Euler-operation sequence",
            self.required_lhs_increase, self.required_lhs_increase
        )
    }
}

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
mod tests {
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    use super::*;

    fn cube_handles() -> (FaceHandle, FaceHandle, VertexHandle, VertexHandle) {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let faces = brep.faces.keys().collect::<Vec<_>>();
        let vertices = brep.vertices.keys().collect::<Vec<_>>();
        (faces[0], faces[1], vertices[0], vertices[1])
    }

    #[test]
    fn subcomplex_plan_owns_support_group_derivation() {
        let (cylinder_face, support_face, start, end) = cube_handles();
        let plan = TrimRewriteSubcomplexPlan::from_edge_plans(
            vec![TrimRewriteEdgePlan::new_generator_trim(
                0,
                cylinder_face,
                support_face,
                start,
                end,
            )],
            TopologyCountDelta::zero(),
        );

        assert_eq!(plan.support_groups().len(), 1);
        assert_eq!(plan.support_groups()[0].support_face(), support_face);
        assert_eq!(plan.support_groups()[0].generated_trim_count(), 1);
        assert_eq!(plan.support_groups()[0].cylinder_faces(), &[cylinder_face]);
    }

    #[test]
    fn planned_generator_support_face_requires_unique_owner() {
        let (cylinder_face, support_face, start, end) = cube_handles();
        let edge_plan =
            TrimRewriteEdgePlan::new_generator_trim(0, cylinder_face, support_face, start, end);
        let plan =
            TrimRewriteSubcomplexPlan::from_edge_plans(vec![edge_plan], TopologyCountDelta::zero());

        assert_eq!(
            plan.planned_generator_support_face(cylinder_face, 0, start, end, "generator")
                .expect("unique owner"),
            support_face
        );

        let duplicate = TrimRewriteSubcomplexPlan::from_edge_plans(
            vec![edge_plan, edge_plan],
            TopologyCountDelta::zero(),
        );
        assert!(
            duplicate
                .planned_generator_support_face(cylinder_face, 0, start, end, "generator")
                .is_err()
        );
    }

    #[test]
    fn subcomplex_builder_owns_generator_delta_and_support_groups() {
        let (cylinder_face, support_face, start, end) = cube_handles();
        let mut builder = TrimRewriteSubcomplexBuilder::new();
        builder.push_new_generator_trim(0, cylinder_face, support_face, start, end);
        let plan = builder.finish();

        assert_eq!(plan.edge_plans().len(), 1);
        assert_eq!(plan.support_groups().len(), 1);
        assert!(
            plan.diagnostic_label()
                .contains("edge-only delta dV=0 dE=1 dF=0 dL=0")
        );
    }

    #[test]
    fn subcomplex_plan_projects_seed_faces_and_trim_intents() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let faces = brep.faces.keys().collect::<Vec<_>>();
        let vertices = brep.vertices.keys().collect::<Vec<_>>();
        let cylinder_face = faces[0];
        let cap_face = faces[1];
        let support_face = faces[2];
        let start = vertices[0];
        let end = vertices[1];

        let mut builder = TrimRewriteSubcomplexBuilder::new();
        builder.push_existing_seam(0, cylinder_face, cap_face, start, end);
        builder.push_new_generator_trim(0, cylinder_face, support_face, start, end);
        let plan = builder.finish();

        let seed_faces = plan
            .seed_faces()
            .into_iter()
            .collect::<std::collections::BTreeSet<_>>();
        assert!(seed_faces.contains(&cylinder_face));
        assert!(seed_faces.contains(&cap_face));
        assert!(seed_faces.contains(&support_face));

        let trim_intents = plan.trim_intents();
        assert_eq!(trim_intents.len(), 2);
        let generator_intent = trim_intents
            .iter()
            .find(|intent| intent.is_new_generator_trim)
            .expect("generator intent is projected");
        assert_eq!(generator_intent.cylinder_face, cylinder_face);
        assert_eq!(generator_intent.support_face, Some(support_face));
        assert_eq!(generator_intent.start, start);
        assert_eq!(generator_intent.end, end);
    }
}
