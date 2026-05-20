use popkern_core::topology::BRep;

#[derive(Clone, Copy)]
pub(in crate::apex_construction) struct TopologyCountSummary {
    vertices: i64,
    edges: i64,
    faces: i64,
    loops: i64,
    shells: i64,
}

impl TopologyCountSummary {
    pub(in crate::apex_construction) fn any_genus_lhs(self) -> i64 {
        self.vertices - self.edges + 2 * self.faces - self.loops
    }

    pub(in crate::apex_construction) fn genus_zero_rhs(self) -> i64 {
        2 * self.shells
    }

    pub(in crate::apex_construction) fn deficit_to_genus_zero(self) -> i64 {
        self.genus_zero_rhs() - self.any_genus_lhs()
    }

    pub(in crate::apex_construction) fn trace_label(self, label: &'static str) -> String {
        format!(
            "[direct-apex-trace §14.9/topology-counts/{label}] V={} E={} F={} L={} S={} lhs={}",
            self.vertices,
            self.edges,
            self.faces,
            self.loops,
            self.shells,
            self.any_genus_lhs()
        )
    }

    pub(in crate::apex_construction) fn diagnostic_label(self) -> String {
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

pub(in crate::apex_construction) fn topology_count_summary(brep: &BRep) -> TopologyCountSummary {
    TopologyCountSummary {
        vertices: brep.vertices.len() as i64,
        edges: brep.edges.len() as i64,
        faces: brep.faces.len() as i64,
        loops: brep.loops.len() as i64,
        shells: brep.shells.len() as i64,
    }
}

pub(in crate::apex_construction) struct TopologyCountDelta {
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
    pub(in crate::apex_construction) fn zero() -> Self {
        Self::from_components(0, 0, 0, 0, 0)
    }

    pub(in crate::apex_construction) fn from_components(
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

    pub(in crate::apex_construction) fn between(
        before: TopologyCountSummary,
        after: TopologyCountSummary,
    ) -> Self {
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

    pub(in crate::apex_construction) fn add_assign(&mut self, other: Self) {
        self.vertices += other.vertices;
        self.edges += other.edges;
        self.faces += other.faces;
        self.loops += other.loops;
        self.shells += other.shells;
        self.lhs += other.lhs;
        self.rhs += other.rhs;
        self.deficit += other.deficit;
    }

    pub(in crate::apex_construction) fn diagnostic_label(&self) -> String {
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

pub(in crate::apex_construction) struct TopologyBalanceRequirement {
    required_lhs_increase: i64,
}

impl TopologyBalanceRequirement {
    pub(in crate::apex_construction) fn from_deficit_increase(deficit_increase: i64) -> Self {
        Self {
            required_lhs_increase: deficit_increase.max(0),
        }
    }

    pub(in crate::apex_construction) fn diagnostic_label(&self) -> String {
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
