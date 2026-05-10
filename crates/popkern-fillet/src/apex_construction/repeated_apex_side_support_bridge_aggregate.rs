//! Aggregate cycle selection for face-family bridge primitives.

use std::collections::BTreeSet;

use popkern_core::topology::{BRep, EdgeHandle};

use super::{FaceFamilyBridgeBoundaryCycle, FaceFamilyBridgeRegionPrimitivePlan};

pub(super) fn select_aggregate_region_cycles(
    brep: &BRep,
    regions: Vec<FaceFamilyBridgeRegionPrimitivePlan>,
) -> Vec<FaceFamilyBridgeRegionPrimitivePlan> {
    if regions.len() <= 1 || !regions_share_one_outside_loop(&regions) {
        return regions;
    }
    let Some(selected) = disjoint_aggregate_cycles(brep, &regions) else {
        return regions;
    };
    regions
        .into_iter()
        .zip(selected)
        .map(|(region, cycle)| region.with_boundary_cycles(brep, vec![cycle]))
        .collect()
}

fn regions_share_one_outside_loop(regions: &[FaceFamilyBridgeRegionPrimitivePlan]) -> bool {
    let Some(first) = regions.first() else {
        return false;
    };
    regions
        .iter()
        .all(|region| region.outside_loop() == first.outside_loop())
}

fn disjoint_aggregate_cycles(
    brep: &BRep,
    regions: &[FaceFamilyBridgeRegionPrimitivePlan],
) -> Option<Vec<FaceFamilyBridgeBoundaryCycle>> {
    let candidates = regions
        .iter()
        .map(|region| aggregate_cycle_candidates(brep, region))
        .collect::<Option<Vec<_>>>()?;
    if candidates.iter().any(Vec::is_empty) {
        return None;
    }
    let mut selected = Vec::with_capacity(regions.len());
    let mut consumed = BTreeSet::<EdgeHandle>::new();
    select_disjoint_cycle(&candidates, 0, &mut consumed, &mut selected).then_some(selected)
}

fn aggregate_cycle_candidates(
    brep: &BRep,
    region: &FaceFamilyBridgeRegionPrimitivePlan,
) -> Option<Vec<(FaceFamilyBridgeBoundaryCycle, BTreeSet<EdgeHandle>)>> {
    let support_intervals = [
        region.support_paths().first_to_second().to_vec(),
        region.support_paths().second_to_first().to_vec(),
    ];
    let outside_intervals = [
        region.outside_paths().first_to_second().to_vec(),
        region.outside_paths().second_to_first().to_vec(),
    ];
    let mut candidates = support_intervals
        .into_iter()
        .flat_map(|support_interval| {
            outside_intervals
                .iter()
                .cloned()
                .map(move |outside_interval| {
                    FaceFamilyBridgeBoundaryCycle::new(support_interval.clone(), outside_interval)
                })
        })
        .filter_map(|cycle| cycle_edge_set(brep, &cycle).map(|edges| (cycle, edges)))
        .collect::<Vec<_>>();
    candidates.sort_by_key(|(cycle, edges)| {
        (
            cycle.support_interval().len() + cycle.outside_interval().len(),
            cycle.support_interval().len(),
            cycle.outside_interval().len(),
            edges.len(),
        )
    });
    Some(candidates)
}

fn select_disjoint_cycle(
    candidates: &[Vec<(FaceFamilyBridgeBoundaryCycle, BTreeSet<EdgeHandle>)>],
    region_index: usize,
    consumed: &mut BTreeSet<EdgeHandle>,
    selected: &mut Vec<FaceFamilyBridgeBoundaryCycle>,
) -> bool {
    if region_index == candidates.len() {
        return true;
    }
    for (cycle, edges) in &candidates[region_index] {
        if edges.iter().any(|edge| consumed.contains(edge)) {
            continue;
        }
        for edge in edges {
            consumed.insert(*edge);
        }
        selected.push(cycle.clone());
        if select_disjoint_cycle(candidates, region_index + 1, consumed, selected) {
            return true;
        }
        selected.pop();
        for edge in edges {
            consumed.remove(edge);
        }
    }
    false
}

fn cycle_edge_set(
    brep: &BRep,
    cycle: &FaceFamilyBridgeBoundaryCycle,
) -> Option<BTreeSet<EdgeHandle>> {
    let mut edges = BTreeSet::new();
    for half_edge in cycle
        .support_interval()
        .iter()
        .take(cycle.support_interval().len().checked_sub(1)?)
        .chain(
            cycle
                .outside_interval()
                .iter()
                .take(cycle.outside_interval().len().checked_sub(1)?),
        )
    {
        edges.insert(brep.half_edges.get(*half_edge)?.edge);
    }
    Some(edges)
}
