use super::*;

pub(super) struct PartitionBoundaryOutcome {
    pub(super) cap_cycle: Vec<HalfEdgeHandle>,
    fillet_end_cap_sides: std::collections::BTreeMap<usize, HalfEdgeHandle>,
    strip_side_overrides: std::collections::BTreeMap<usize, HalfEdgeHandle>,
    strip_loop_overrides: std::collections::BTreeMap<usize, Vec<HalfEdgeHandle>>,
    existing_target_cylinder_strips: std::collections::BTreeSet<usize>,
}

impl PartitionBoundaryOutcome {
    pub(super) fn from_current_cap_cycle(cap_cycle: Vec<HalfEdgeHandle>) -> Self {
        Self {
            cap_cycle,
            fillet_end_cap_sides: std::collections::BTreeMap::new(),
            strip_side_overrides: std::collections::BTreeMap::new(),
            strip_loop_overrides: std::collections::BTreeMap::new(),
            existing_target_cylinder_strips: std::collections::BTreeSet::new(),
        }
    }

    pub(super) fn from_reauthored_cap_cycle(
        cap_cycle: Vec<HalfEdgeHandle>,
        fillet_end_cap_sides: std::collections::BTreeMap<usize, HalfEdgeHandle>,
        strip_side_overrides: std::collections::BTreeMap<usize, HalfEdgeHandle>,
    ) -> Self {
        Self {
            cap_cycle,
            fillet_end_cap_sides,
            strip_side_overrides,
            strip_loop_overrides: std::collections::BTreeMap::new(),
            existing_target_cylinder_strips: std::collections::BTreeSet::new(),
        }
    }

    pub(super) fn mark_existing_target_cylinder_strip(&mut self, fillet_index: usize) {
        self.existing_target_cylinder_strips.insert(fillet_index);
    }

    pub(super) fn uses_existing_target_cylinder_strip(&self, fillet_index: usize) -> bool {
        self.existing_target_cylinder_strips.contains(&fillet_index)
    }

    pub(super) fn strip_side_override(&self, fillet_index: usize) -> Option<HalfEdgeHandle> {
        self.strip_side_overrides.get(&fillet_index).copied()
    }

    pub(super) fn fillet_end_cap_side(&self, fillet_index: usize) -> Option<HalfEdgeHandle> {
        self.fillet_end_cap_sides.get(&fillet_index).copied()
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
        let cap_sides = self
            .fillet_end_cap_sides
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
        let existing_strips = self
            .existing_target_cylinder_strips
            .iter()
            .map(|fillet_index| format!("FilletEnd[{fillet_index}]"))
            .collect::<Vec<_>>()
            .join(", ");
        format!(
            "selected cap_cycle={:?}; strip-side overrides=[{}]; strip-loop overrides=[{}]; \
             fillet-end cap-sides=[{}]; existing target-cylinder strips=[{}]",
            self.cap_cycle, overrides, strip_loops, cap_sides, existing_strips
        )
    }

    pub(super) fn matches_current_cap_cycle(&self, current_cap_cycle: &[HalfEdgeHandle]) -> bool {
        self.cap_cycle == current_cap_cycle
    }

    pub(super) fn refresh_from_current_cap_loop(
        &mut self,
        brep: &BRep,
        cap_face: FaceHandle,
    ) -> Result<(), ApexBuildError> {
        let cap_loop = brep
            .faces
            .get(cap_face)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition selected cap boundary references stale cap \
                     face {cap_face:?}"
                ),
            })?
            .outer_loop;
        let actual = collect_loop_half_edges(brep, cap_loop)?;
        if actual == self.cap_cycle {
            return Ok(());
        }
        if !cap_cycles_are_geometrically_equivalent(brep, &actual, &self.cap_cycle)
            && !cap_cycles_are_single_handle_refresh(&actual, &self.cap_cycle)
        {
            return Ok(());
        }
        for (actual_side, expected_side) in
            actual.iter().copied().zip(self.cap_cycle.iter().copied())
        {
            for cap_side in self.fillet_end_cap_sides.values_mut() {
                if *cap_side == expected_side {
                    *cap_side = actual_side;
                }
            }
            let strip_replacements = self
                .strip_side_overrides
                .iter_mut()
                .filter_map(|(fillet_index, strip_side)| {
                    if *strip_side == expected_side {
                        *strip_side = actual_side;
                        Some(*fillet_index)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
            for fillet_index in strip_replacements {
                if let Some(cap_side) = self.fillet_end_cap_sides.get_mut(&fillet_index) {
                    *cap_side = actual_side;
                }
            }
            for strip_loop in self.strip_loop_overrides.values_mut() {
                for strip_side in strip_loop {
                    if *strip_side == expected_side {
                        *strip_side = actual_side;
                    }
                }
            }
        }
        self.cap_cycle = actual;
        Ok(())
    }

    pub(super) fn plan_strip_loop_overrides(
        &mut self,
        brep: &BRep,
        fillet_ends: &[FilletEndOwnership],
    ) -> Result<(), ApexBuildError> {
        for ownership in fillet_ends {
            if self.uses_existing_target_cylinder_strip(ownership.fillet_index) {
                continue;
            }
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

fn cap_cycles_are_geometrically_equivalent(
    brep: &BRep,
    actual: &[HalfEdgeHandle],
    expected: &[HalfEdgeHandle],
) -> bool {
    actual.len() == expected.len()
        && actual.iter().copied().zip(expected.iter().copied()).all(
            |(actual_side, expected_side)| {
                cap_boundary_side_is_geometrically_equivalent(brep, actual_side, expected_side)
            },
        )
}

fn cap_cycles_are_single_handle_refresh(
    actual: &[HalfEdgeHandle],
    expected: &[HalfEdgeHandle],
) -> bool {
    actual.len() == expected.len()
        && !actual.is_empty()
        && actual
            .iter()
            .zip(expected)
            .filter(|(actual_side, expected_side)| actual_side == expected_side)
            .count()
            + 1
            >= expected.len()
}

fn cap_boundary_side_is_geometrically_equivalent(
    brep: &BRep,
    actual: HalfEdgeHandle,
    expected: HalfEdgeHandle,
) -> bool {
    let Some(actual_start) = brep
        .half_edges
        .get(actual)
        .map(|half_edge| half_edge.origin)
    else {
        return false;
    };
    let Some(expected_start) = brep
        .half_edges
        .get(expected)
        .map(|half_edge| half_edge.origin)
    else {
        return false;
    };
    let Some(actual_end) = half_edge_end_vertex(brep, actual) else {
        return false;
    };
    let Some(expected_end) = half_edge_end_vertex(brep, expected) else {
        return false;
    };
    vertices_are_coincident(brep, actual_start, expected_start)
        && vertices_are_coincident(brep, actual_end, expected_end)
}
