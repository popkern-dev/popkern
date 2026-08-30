//! Honest attribution of apply-phase failures on ridge-completed
//! strips (PK-380).
//!
//! Two non-adjacent fillets whose removal strips cross on a shared
//! support face (rails crossed: `r_a + r_b > gap`) are completed by
//! the parallel-cylinder seam construction (PK-159) instead of being
//! declined. Every per-edge radius gate has already passed for such
//! an edge, so when the apply phase then fails on it the only
//! radius-related cause left is the interference itself. Reporting
//! the trim stage's generic `InvalidTopology` ("heal the input") for a
//! valid body was measured on the oversized cube fillet; the honest
//! payload is [`BlendFailure::RadiusTooLarge`] with the proportional
//! bound `r_a · gap / (r_a + r_b)` — the radius at which both strips,
//! shrunk together, exactly share the face.

use crate::construct::EdgeGeometry;
use crate::error::BlendFailure;
use crate::overlap_strip_interference::StripInterference;
use crate::{FilletRadius, FilletSpec};
use popkern_core::topology::{BRep, EdgeHandle};

/// Peak radius of a spec.
fn requested_radius(radius: &FilletRadius) -> f64 {
    crate::geometry_compute_feasibility::largest_requested_radius(radius)
}

/// The tightest proportional-shrink bound over every ridge-completed
/// pair the geometry at `index` belongs to, with the request itself.
fn ridge_bound(
    index: usize,
    ridge_completed: &[StripInterference],
    specs: &[FilletSpec],
) -> Option<(f64, f64)> {
    let r_self = requested_radius(&specs.get(index)?.radius);
    let mut best: Option<f64> = None;
    for hit in ridge_completed {
        let other = if hit.fillet_a == index {
            hit.fillet_b
        } else if hit.fillet_b == index {
            hit.fillet_a
        } else {
            continue;
        };
        let r_other = requested_radius(&specs.get(other)?.radius);
        let sum = r_self + r_other;
        if !(sum > 0.0 && hit.gap.is_finite() && hit.gap > 0.0) {
            continue;
        }
        let bound = r_self * hit.gap / sum;
        best = Some(best.map_or(bound, |b: f64| b.min(bound)));
    }
    best.map(|b| (r_self, b))
}

/// Replace apply-phase `InvalidTopology` reports on ridge-completed
/// edges with `RadiusTooLarge { requested, max_radius }` located at the
/// input edge's midpoint; every other failure passes through
/// untouched. Order is preserved.
pub(crate) fn attribute_ridge_completed_failures(
    input: &BRep,
    failures: Vec<(EdgeHandle, BlendFailure)>,
    ridge_completed: &[StripInterference],
    specs: &[FilletSpec],
    edge_geometries: &[EdgeGeometry],
) -> Vec<(EdgeHandle, BlendFailure)> {
    if ridge_completed.is_empty() {
        return failures;
    }
    failures
        .into_iter()
        .map(|(edge, failure)| {
            let BlendFailure::InvalidTopology { location, .. } = &failure else {
                return (edge, failure);
            };
            let Some(index) = edge_geometries.iter().position(|g| g.analysis.edge == edge) else {
                return (edge, failure);
            };
            match ridge_bound(index, ridge_completed, specs) {
                Some((requested, max_radius)) => (
                    edge,
                    BlendFailure::RadiusTooLarge {
                        requested,
                        max_radius: Some(max_radius),
                        // The input edge's midpoint: the working copy
                        // has already been split by the trim stage.
                        location: crate::geometry_compute::edge_midpoint_by_handle(input, edge)
                            .or(*location),
                    },
                ),
                None => (edge, failure),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(edge: EdgeHandle, r: f64) -> FilletSpec {
        FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        }
    }

    #[test]
    fn proportional_bound_shares_the_gap() {
        let e = EdgeHandle::default();
        let specs = vec![spec(e, 0.9), spec(e, 0.9)];
        let hits = vec![StripInterference {
            fillet_a: 0,
            fillet_b: 1,
            gap: 1.0,
        }];
        let (requested, bound) = ridge_bound(0, &hits, &specs).expect("member");
        assert_eq!(requested, 0.9);
        assert!((bound - 0.5).abs() < 1e-15, "bound {bound}");
        // Unequal radii shrink proportionally: 0.3 and 0.9 over 1.0.
        let specs = vec![spec(e, 0.3), spec(e, 0.9)];
        let (_, b0) = ridge_bound(0, &hits, &specs).expect("member");
        let (_, b1) = ridge_bound(1, &hits, &specs).expect("member");
        assert!((b0 - 0.25).abs() < 1e-15 && (b1 - 0.75).abs() < 1e-15);
        assert!(
            (b0 + b1 - 1.0).abs() < 1e-15,
            "shrunk strips exactly share the gap"
        );
        // Non-members and degenerate gaps attribute nothing.
        assert!(ridge_bound(2, &hits, &specs).is_none());
        let zero = vec![StripInterference {
            fillet_a: 0,
            fillet_b: 1,
            gap: 0.0,
        }];
        assert!(ridge_bound(0, &zero, &specs).is_none());
    }
}
