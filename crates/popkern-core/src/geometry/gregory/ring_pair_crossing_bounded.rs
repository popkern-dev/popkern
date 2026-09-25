//! The bounded ring-crossing verdict for optimizers holding a fallback
//! candidate (PK-547): no surface-surface solve on spoke-neighbour
//! pairs.

use super::super::RectangularGregoryBiNonic;
use super::ring_pair_crossing_certificate::pair_certified_clean;
use super::ring_pair_crossing_setup::{pair_setup, patch_aabbs, spoke_polylines};
use super::{BinonicRingCrossing, EXEMPTION_BAND_FACTOR, RingPairRelation, pair_crossing};
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

/// A ring-crossing verdict that never runs the surface-surface solver on
/// a spoke-neighbour pair (PK-547). See [`binonic_ring_crossing_bounded`].
#[derive(Debug, Clone)]
pub enum BoundedRingCrossing {
    /// Every pair is certified clean, or, for an apex-only pair, settled
    /// clean by the solver.
    Clean,
    /// A proven crossing.
    Crossing(BinonicRingCrossing),
    /// A spoke-neighbour pair the separation certificate could not clear:
    /// a crossing, or a near-crossing past the exemption band.
    Undecided {
        /// The pair `(a, b)`, `a < b`.
        pair: (usize, usize),
    },
}

/// [`binonic_ring_crossing`] for optimizers that hold a fallback
/// candidate. Spoke-neighbour pairs are judged by the separation
/// certificate alone. Legitimate spoke contact lies inside its macro
/// exemption, and the served rings certify in microseconds. A pair it
/// cannot clear is reported [`BoundedRingCrossing::Undecided`] rather
/// than settled by the solver, which took 223 s on a pyramid4_uniform
/// candidate whose crossing sat 4% past the band. Apex-only pairs, whose
/// point contact the strict band cannot always clear, fall through to
/// the solver as in the full gate.
pub fn binonic_ring_crossing_bounded(
    ring: &[RectangularGregoryBiNonic],
    tolerance: f64,
) -> BoundedRingCrossing {
    let n = ring.len();
    if n < 2 {
        return BoundedRingCrossing::Clean;
    }
    let band = EXEMPTION_BAND_FACTOR * tolerance;
    let aabbs = patch_aabbs(ring);
    let spokes = spoke_polylines(ring);
    let ssi_tol = Tolerance::default().with_linear(tolerance.max(1e-6));
    let mut nurbs: Vec<Option<SurfaceGeometry>> = vec![None; n];
    for a in 0..n {
        for b in (a + 1)..n {
            if let (Some(bb_a), Some(bb_b)) = (&aabbs[a], &aabbs[b])
                && !bb_a.overlaps(bb_b, band)
            {
                continue;
            }
            let setup = pair_setup(ring, &aabbs, &spokes, (a, b), band, false);
            if pair_certified_clean(
                &ring[a].grid,
                &ring[b].grid,
                &setup.shared_spokes,
                setup.apex,
                band,
                setup.exempt_depth,
            ) {
                continue;
            }
            if setup.relation == RingPairRelation::SpokeNeighbors {
                return BoundedRingCrossing::Undecided { pair: (a, b) };
            }
            for k in [a, b] {
                if nurbs[k].is_none() {
                    nurbs[k] = Some(SurfaceGeometry::Nurbs(ring[k].to_nurbs_exact()));
                }
            }
            let (Some(surf_a), Some(surf_b)) = (&nurbs[a], &nurbs[b]) else {
                continue;
            };
            if let Some(crossing) = pair_crossing(
                surf_a,
                surf_b,
                (a, b),
                setup.relation,
                &setup.shared_spokes,
                setup.apex,
                band,
                setup.macro_depth,
                false,
                &ssi_tol,
            ) {
                return BoundedRingCrossing::Crossing(crossing);
            }
        }
    }
    BoundedRingCrossing::Clean
}
