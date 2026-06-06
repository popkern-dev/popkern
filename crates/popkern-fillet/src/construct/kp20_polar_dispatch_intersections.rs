//! Self-intersection screening helpers for KP20 polar dispatch.

use popkern_core::geometry::gregory::{GregoryPatch, GregorySelfIntersectionReport};

use super::super::kp_polar_bi6_promote::SELF_INT_TOL;

type IndexedIntersection = (usize, GregorySelfIntersectionReport);

pub(super) fn fast_self_intersecting_patches(
    candidate: &[GregoryPatch],
) -> Vec<IndexedIntersection> {
    candidate
        .iter()
        .enumerate()
        .filter_map(|(idx, gp)| {
            let report = gp.jacobian_self_intersection_report();
            report.intersects().then_some((idx, report))
        })
        .collect()
}

pub(super) fn exact_self_intersecting_patches(
    candidate: &[GregoryPatch],
) -> Vec<IndexedIntersection> {
    candidate
        .iter()
        .enumerate()
        .filter_map(|(idx, gp)| {
            let report = gp.self_intersection_report(SELF_INT_TOL);
            report.intersects().then_some((idx, report))
        })
        .collect()
}

pub(super) fn acceptance_screen(candidate: &[GregoryPatch]) -> Vec<IndexedIntersection> {
    let fast = fast_self_intersecting_patches(candidate);
    if fast.is_empty() {
        exact_self_intersecting_patches(candidate)
    } else {
        fast
    }
}
