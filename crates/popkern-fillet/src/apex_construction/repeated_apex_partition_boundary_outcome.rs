use super::*;

pub(super) struct PartitionBoundaryOutcome {
    pub(super) cap_cycle: Vec<HalfEdgeHandle>,
    strip_side_overrides: std::collections::BTreeMap<usize, HalfEdgeHandle>,
    strip_loop_overrides: std::collections::BTreeMap<usize, Vec<HalfEdgeHandle>>,
}

impl PartitionBoundaryOutcome {
    pub(super) fn from_current_cap_cycle(cap_cycle: Vec<HalfEdgeHandle>) -> Self {
        Self {
            cap_cycle,
            strip_side_overrides: std::collections::BTreeMap::new(),
            strip_loop_overrides: std::collections::BTreeMap::new(),
        }
    }

    pub(super) fn from_reauthored_cap_cycle(
        cap_cycle: Vec<HalfEdgeHandle>,
        strip_side_overrides: std::collections::BTreeMap<usize, HalfEdgeHandle>,
    ) -> Self {
        Self {
            cap_cycle,
            strip_side_overrides,
            strip_loop_overrides: std::collections::BTreeMap::new(),
        }
    }

    pub(super) fn strip_loop(
        &self,
        fillet_index: usize,
        ownership: &FilletEndOwnership,
    ) -> Vec<HalfEdgeHandle> {
        if let Some(loop_half_edges) = self.strip_loop_overrides.get(&fillet_index) {
            return loop_half_edges.clone();
        }
        let shared = self
            .strip_side_overrides
            .get(&fillet_index)
            .copied()
            .unwrap_or(ownership.cylinder_candidate_half_edge);
        vec![
            shared,
            ownership.retired_apex_cell.t_a_to_apex,
            ownership.retired_apex_cell.apex_to_t_b,
        ]
    }

    pub(super) fn diagnostic_label(&self) -> String {
        let overrides = self
            .strip_side_overrides
            .iter()
            .map(|(fillet_index, half_edge)| format!("FilletEnd[{fillet_index}]={half_edge:?}"))
            .collect::<Vec<_>>()
            .join(", ");
        let strip_loops = self
            .strip_loop_overrides
            .iter()
            .map(|(fillet_index, half_edges)| format!("FilletEnd[{fillet_index}]={half_edges:?}"))
            .collect::<Vec<_>>()
            .join(", ");
        format!(
            "selected cap_cycle={:?}; strip-side overrides=[{}]; strip-loop overrides=[{}]",
            self.cap_cycle, overrides, strip_loops
        )
    }

    pub(super) fn matches_current_cap_cycle(&self, current_cap_cycle: &[HalfEdgeHandle]) -> bool {
        self.cap_cycle == current_cap_cycle
    }

    pub(super) fn plan_strip_loop_overrides(
        &mut self,
        brep: &BRep,
        fillet_ends: &[FilletEndOwnership],
    ) -> Result<(), ApexBuildError> {
        for ownership in fillet_ends {
            let candidates = self.strip_loop_candidates(brep, ownership)?;
            let mut rejected = Vec::new();
            let mut selected = None;
            for candidate in candidates {
                match validate_half_edge_cycle(
                    brep,
                    &candidate,
                    "selected FilletEnd strip boundary",
                    ownership.target_cylinder_face,
                ) {
                    Ok(()) => {
                        selected = Some(candidate);
                        break;
                    }
                    Err(err) => rejected.push(format!("{candidate:?}: {err:?}")),
                }
            }
            let Some(selected) = selected else {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition target-cylinder trim rewrite cannot yet \
                         mutate: no contiguous strip loop candidate for FilletEnd[{}]; rejected {}",
                        ownership.fillet_index,
                        rejected.join("; ")
                    ),
                });
            };
            self.strip_loop_overrides
                .insert(ownership.fillet_index, selected);
        }
        Ok(())
    }

    fn strip_loop_candidates(
        &self,
        brep: &BRep,
        ownership: &FilletEndOwnership,
    ) -> Result<Vec<Vec<HalfEdgeHandle>>, ApexBuildError> {
        let shared = self
            .strip_side_overrides
            .get(&ownership.fillet_index)
            .copied()
            .unwrap_or(ownership.cylinder_candidate_half_edge);
        let ta_to_apex = ownership.retired_apex_cell.t_a_to_apex;
        let apex_to_tb = ownership.retired_apex_cell.apex_to_t_b;
        let apex_to_ta = brep
            .half_edges
            .get(ta_to_apex)
            .map(|half_edge| half_edge.twin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] strip candidate references \
                     stale T_a->apex half-edge {ta_to_apex:?}",
                    ownership.fillet_index
                ),
            })?;
        let tb_to_apex = brep
            .half_edges
            .get(apex_to_tb)
            .map(|half_edge| half_edge.twin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] strip candidate references \
                     stale apex->T_b half-edge {apex_to_tb:?}",
                    ownership.fillet_index
                ),
            })?;
        Ok(vec![
            vec![shared, ta_to_apex, apex_to_tb],
            vec![shared, tb_to_apex, apex_to_ta],
        ])
    }
}
