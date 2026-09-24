//! Shared-boundary-CP re-unification for the setback cap fan.
//!
//! `reconcile_cap_rim_to_sweeps` lands each rim patch's boundary
//! rows on the live sweeps by moving its depth-0..2 rows. Those
//! moves also displace the END control points of the patch's
//! PERPENDICULAR (spoke) boundary rows — which adjacent patches
//! share. A neighbour that is not itself rim-moved (or is moved by
//! a different per-patch Bernstein solve) keeps stale copies, and
//! the fan tears: measured 2.5e-2 G0 on hexagon6_variable spoke
//! seams near rim corners while whole-side gaps elsewhere sit at
//! 1e-15..3.5e-4 (PK-087 stage 14b).
//!
//! The patches are assembled from a shared net, so boundary CPs
//! that must remain identical are exactly equal BEFORE the rim
//! moves. This pass records those equivalence groups first; after
//! the moves every group is re-unified to the mean of its MOVED
//! members (movers carry the on-sweep positions; stale members
//! adopt them).

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::GregoryPatch;

/// Boundary-CP slots of a bi-9 grid: i or j on the 0/9 border.
fn boundary_slots() -> impl Iterator<Item = (usize, usize)> {
    (0..10usize).flat_map(|i| {
        (0..10usize).filter_map(move |j| (i == 0 || i == 9 || j == 0 || j == 9).then_some((i, j)))
    })
}

/// Equivalence groups of exactly-shared boundary CPs, recorded from
/// the pre-move grids. Key: quantized position.
pub(super) struct SharedBoundaryGroups {
    groups: Vec<Vec<(usize, usize, usize)>>,
    original: Vec<Vec<Point3>>,
}

const QUANTUM: f64 = 1.0e-9;

fn key_of(p: &Point3) -> (i64, i64, i64) {
    (
        (p.x / QUANTUM).round() as i64,
        (p.y / QUANTUM).round() as i64,
        (p.z / QUANTUM).round() as i64,
    )
}

/// `(patch index, grid i, grid j)` of one shared-CP group member.
type GroupMember = (usize, usize, usize);

pub(super) fn record_shared_boundary_groups(patches: &[GregoryPatch]) -> SharedBoundaryGroups {
    let mut map: std::collections::BTreeMap<(i64, i64, i64), Vec<GroupMember>> =
        std::collections::BTreeMap::new();
    for (pi, patch) in patches.iter().enumerate() {
        let GregoryPatch::RectangularBiNonic(bp) = patch else {
            continue;
        };
        for (i, j) in boundary_slots() {
            map.entry(key_of(&bp.grid[i][j]))
                .or_default()
                .push((pi, i, j));
        }
    }
    let groups: Vec<Vec<(usize, usize, usize)>> = map
        .into_values()
        .filter(|members| {
            // Shared across at least two PATCHES (a patch's own
            // corner appears once per patch here).
            members
                .iter()
                .map(|(pi, _, _)| pi)
                .collect::<std::collections::BTreeSet<_>>()
                .len()
                >= 2
        })
        .collect();
    let original = groups
        .iter()
        .map(|members| {
            members
                .iter()
                .map(|&(pi, i, j)| match &patches[pi] {
                    GregoryPatch::RectangularBiNonic(bp) => bp.grid[i][j],
                    _ => Point3::new(0.0, 0.0, 0.0),
                })
                .collect()
        })
        .collect();
    SharedBoundaryGroups { groups, original }
}

/// Re-unify every recorded group: members that moved carry the new
/// (on-sweep) positions; all members adopt the movers' mean. Returns
/// the worst pre-unification disagreement within any group — the
/// measured size of the tear this pass closed.
pub(super) fn reunify_shared_boundary_cps(
    patches: &mut [GregoryPatch],
    recorded: &SharedBoundaryGroups,
) -> f64 {
    let mut worst_tear = 0.0_f64;
    for (members, originals) in recorded.groups.iter().zip(&recorded.original) {
        let mut moved_sum = [0.0_f64; 3];
        let mut moved_count = 0usize;
        let mut spread = 0.0_f64;
        let mut positions = Vec::with_capacity(members.len());
        for (&(pi, i, j), original) in members.iter().zip(originals) {
            let Some(grid) = patches[pi].binonic_grid() else {
                continue;
            };
            let p = grid[i][j];
            if p.distance_to(original) > 1e-12 {
                moved_sum[0] += p.x;
                moved_sum[1] += p.y;
                moved_sum[2] += p.z;
                moved_count += 1;
            }
            positions.push(p);
        }
        for a in &positions {
            for b in &positions {
                spread = spread.max(a.distance_to(b));
            }
        }
        if moved_count == 0 || spread <= 1e-12 {
            continue;
        }
        worst_tear = worst_tear.max(spread);
        let target = Point3::new(
            moved_sum[0] / moved_count as f64,
            moved_sum[1] / moved_count as f64,
            moved_sum[2] / moved_count as f64,
        );
        for &(pi, i, j) in members {
            if let Some(grid) = patches[pi].binonic_grid_mut() {
                grid[i][j] = target;
            }
        }
    }
    worst_tear
}
