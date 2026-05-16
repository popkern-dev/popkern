//! Public API entry points for shared BRep repair knowledge.

use crate::topology::{
    BRep, BrepRepairKnowledgeGraph, BrepRepairKnowledgeOptions, BrepRepairKnowledgeScope,
    LocalToleranceSource,
};

/// Analyze a BRep with a scoped, non-mutating repair knowledge collector.
pub fn analyze_brep_repair_knowledge(
    brep: &BRep,
    scope: BrepRepairKnowledgeScope,
    options: BrepRepairKnowledgeOptions,
) -> BrepRepairKnowledgeGraph {
    let samples = options.sample_budget.max(2);
    BrepRepairKnowledgeGraph::build_with_local_tolerance(
        brep,
        scope,
        brep.tolerance.linear(),
        LocalToleranceSource::Model,
        samples,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BrepRepairScopeSeed;

    #[test]
    fn public_api_builds_non_mutating_scoped_graph() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let before = crate::topology::BrepKnowledgeFingerprint::from_brep(&brep);
        let graph = analyze_brep_repair_knowledge(
            &brep,
            BrepRepairKnowledgeScope::seeds_only(vec![BrepRepairScopeSeed::Face(face)]),
            BrepRepairKnowledgeOptions {
                sample_budget: 5,
                ..Default::default()
            },
        );

        assert_eq!(
            crate::topology::BrepKnowledgeFingerprint::from_brep(&brep),
            before
        );
        assert!(graph.validate_fresh(&brep));
        assert!(!graph.facts.is_empty());
        assert!(!graph.diagnostics.is_empty());
    }
}
