//! Per-pair exemption geometry shared by the ring crossing gate and its
//! cheap crossing witness: patch AABBs, spoke polylines, and each
//! pair's relation, shared spokes, apex and exemption depths.

use super::super::RectangularGregoryBiNonic;
use super::{MACRO_DEPTH_FACTOR, RingPairRelation, SPOKE_POLYLINE_SAMPLES, SharedSpoke};
use crate::geometry::Point3;
use crate::geometry::intersection::surface_subdivision::Aabb;

/// The exemption geometry of one patch pair (see the module docs).
pub(super) struct PairSetup<'r> {
    pub(super) relation: RingPairRelation,
    pub(super) shared_spokes: Vec<SharedSpoke<'r>>,
    pub(super) apex: Point3,
    pub(super) macro_depth: f64,
    pub(super) exempt_depth: f64,
}

/// Grid AABBs of every ring patch: tensor-product Bézier control grids
/// have the convex-hull property, so the grid AABB bounds the surface.
pub(super) fn patch_aabbs(ring: &[RectangularGregoryBiNonic]) -> Vec<Option<Aabb>> {
    ring.iter()
        .map(|patch| {
            let points: Vec<Point3> = patch.grid.iter().flatten().copied().collect();
            Aabb::from_points(&points)
        })
        .collect()
}

/// Spoke `k` polyline: patch `k`'s `u = 1` boundary iso, sampled once.
pub(super) fn spoke_polylines(ring: &[RectangularGregoryBiNonic]) -> Vec<Vec<Point3>> {
    ring.iter()
        .map(|patch| {
            (0..SPOKE_POLYLINE_SAMPLES)
                .map(|s| {
                    let t = s as f64 / (SPOKE_POLYLINE_SAMPLES - 1) as f64;
                    patch.point_at(1.0, t)
                })
                .collect()
        })
        .collect()
}

/// Relation, shared spokes, apex and exemption depths of the pair
/// `(a, b)`: spoke `k` glues patch `k`'s UOne side to patch
/// `(k + 1) % n`'s VZero side, and both ring directions are checked
/// mod `n`.
pub(super) fn pair_setup<'r>(
    ring: &'r [RectangularGregoryBiNonic],
    aabbs: &[Option<Aabb>],
    spokes: &'r [Vec<Point3>],
    (a, b): (usize, usize),
    band: f64,
    macro_only: bool,
) -> PairSetup<'r> {
    let n = ring.len();
    let mut shared_spokes: Vec<SharedSpoke> = Vec::with_capacity(2);
    if (a + 1) % n == b {
        shared_spokes.push(SharedSpoke {
            polyline: &spokes[a],
            patch: &ring[a],
        });
    }
    if (b + 1) % n == a {
        shared_spokes.push(SharedSpoke {
            polyline: &spokes[b],
            patch: &ring[b],
        });
    }
    let relation = if shared_spokes.is_empty() {
        RingPairRelation::ApexOnly
    } else {
        RingPairRelation::SpokeNeighbors
    };
    let apex = ring[a].point_at(1.0, 0.0);
    // Rail-hover depth cap, scale-relative to the pair's joint extent
    // (see MACRO_DEPTH_FACTOR).
    let macro_depth = match (&aabbs[a], &aabbs[b]) {
        (Some(bb_a), Some(bb_b)) => {
            let joint = Aabb {
                min: Point3::new(
                    bb_a.min.x.min(bb_b.min.x),
                    bb_a.min.y.min(bb_b.min.y),
                    bb_a.min.z.min(bb_b.min.z),
                ),
                max: Point3::new(
                    bb_a.max.x.max(bb_b.max.x),
                    bb_a.max.y.max(bb_b.max.y),
                    bb_a.max.z.max(bb_b.max.z),
                ),
            };
            MACRO_DEPTH_FACTOR * joint.diagonal()
        }
        _ => band,
    };
    let exempt_depth = if macro_only || relation == RingPairRelation::SpokeNeighbors {
        macro_depth.max(band)
    } else {
        band
    };
    PairSetup {
        relation,
        shared_spokes,
        apex,
        macro_depth,
        exempt_depth,
    }
}
