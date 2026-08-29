//! PK-192 — does the PnS cap beat the Coons fan on the low valences?
//!
//! The ticket's rule is explicit: change dispatch **only** if the
//! measured numbers win, and record a falsification if they do not.
//! These gates hold the candidate's construction contract; the
//! end-to-end A/B against the shipped fan lives in
//! `tests/pns_low_valence_gates.rs`, which can drive the whole
//! pipeline.

use super::*;
use popkern_core::geometry::fill_n_sided_boundary;

/// A synthetic ring standing in for a corner's sectors: `n` sides,
/// each spanning an arc-then-bridge with a plausible surround.
fn synthetic_ring(n: usize) -> Vec<FillBoundaryEdge> {
    let mut ring = Vec::with_capacity(n);
    for k in 0..n {
        let a = |i: f64| 2.0 * std::f64::consts::PI * i / (n as f64);
        let p = |t: f64| Point3::new(f64::cos(t), f64::sin(t), 0.0);
        let (start, mid, end) = (p(a(k as f64)), p(a(k as f64 + 0.5)), p(a(k as f64 + 1.0)));
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let tangent = (end - start)
            .normalized()
            .expect("synthetic side has length");
        ring.push(FillBoundaryEdge {
            start,
            mid,
            end,
            normal_start: normal,
            normal_mid: normal,
            normal_end: normal,
            tangent_start: tangent,
            tangent_mid: tangent,
            tangent_end: tangent,
            extra_samples: Vec::new(),
        });
    }
    ring
}

/// The candidate's whole point: one patch per sector, so it maps
/// one-to-one onto the cap faces the exact-boundary route already
/// commits. A mismatch here would mean the surfaces could not be
/// committed without changing topology.
#[test]
fn one_patch_per_sector_at_the_low_valences() {
    for n in [3_usize, 4, 5] {
        let fill = fill_n_sided_boundary(&synthetic_ring(n), 1e-9)
            .unwrap_or_else(|e| panic!("n={n} ring must fill, got {e}"));
        assert_eq!(
            fill.surfaces.len(),
            n,
            "n={n}: the cap must produce exactly one patch per sector, \
             or it cannot replace the fan without a topology change"
        );
    }
}

/// The candidate interpolates the corner boundary it was handed.
///
/// A cap that does not meet the existing rim cannot be substituted
/// for the fan at all, whatever its internal continuity.
#[test]
fn candidate_interpolates_the_corner_boundary() {
    for n in [3_usize, 4, 5, 6] {
        let ring = synthetic_ring(n);
        let fill = fill_n_sided_boundary(&ring, 1e-9).expect("ring fills");
        assert!(
            fill.worst_boundary_gap < 1e-12,
            "n={n}: cap deviates from the rim by {:.3e}",
            fill.worst_boundary_gap
        );
    }
}

/// A sector pair whose arc and bridge do not meet must decline with
/// the measured gap rather than silently building a cap over a rim
/// that is not closed.
#[test]
fn a_broken_sector_join_is_reported_with_its_gap() {
    let mut ring = synthetic_ring(4);
    ring[2].end = Point3::new(9.0, 9.0, 9.0);
    match fill_n_sided_boundary(&ring, 1e-6) {
        Err(NSidedFillError::BoundaryNotClosed { side, gap, .. }) => {
            assert_eq!(side, 2);
            assert!(gap > 1.0, "the reported gap must be the real one");
        }
        other => panic!("a broken ring must decline, got {other:?}"),
    }
}
