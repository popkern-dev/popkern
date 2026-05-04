//! Trim-subcomplex planning diagnostics for repeated-apex partition rewrites.

use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

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
    pub(super) edge_plans: Vec<TrimRewriteEdgePlan>,
    pub(super) expected_edge_only_delta: TopologyCountDelta,
    pub(super) support_groups: Vec<TrimRewriteSupportGroup>,
}

pub(super) struct TrimRewriteEdgePlan {
    pub(super) fillet_index: usize,
    pub(super) owner_pair: TrimOwnerPair,
    pub(super) construction: TrimEdgeConstruction,
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
}

pub(super) struct TrimRewriteSupportGroup {
    support_face: FaceHandle,
    generated_trim_count: usize,
    cylinder_faces: Vec<FaceHandle>,
}

#[derive(Clone, Copy)]
pub(super) enum TrimOwnerPair {
    CylinderCap {
        cylinder_face: FaceHandle,
        cap_face: FaceHandle,
    },
    CylinderSupport {
        cylinder_face: FaceHandle,
        support_face: FaceHandle,
    },
}

#[derive(Clone, Copy)]
pub(super) enum TrimEdgeConstruction {
    ExistingSeam,
    NewExactGeneratorTrim,
    NewExactSupportCylinderTrim,
}

pub(super) fn trim_rewrite_support_groups(
    edge_plans: &[TrimRewriteEdgePlan],
) -> Vec<TrimRewriteSupportGroup> {
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

impl TrimRewriteSubcomplexPlan {
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
    fn diagnostic_label(&self) -> String {
        format!(
            "{:?}: {} generated trim(s) from cylinder face(s) {:?}",
            self.support_face, self.generated_trim_count, self.cylinder_faces
        )
    }
}

impl TrimRewriteEdgePlan {
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
}

impl TrimOwnerPair {
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
    pub(super) fn is_generated(&self) -> bool {
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
