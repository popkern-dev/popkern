//! Subdivision separation certificate for the ring pairwise crossing
//! gate — proves a patch pair crossing-free without the surface-surface
//! solver (see the parent module's "Cost structure" docs).

use super::{RectangularGregoryBiNonic, SharedSpoke};
use crate::geometry::Point3;
use crate::geometry::intersection::surface_subdivision::Aabb;

/// Subdivision limit for the separation certificate. Both grids halve
/// along one axis per level, so the box diameter halves every two
/// levels: 40 levels reach ~1e-6 of the original extent — below the
/// strict band from any practical ring scale.
const CERTIFICATE_MAX_DEPTH: usize = 40;
/// Node budget for one pair's certificate. Clean rings certify in a few
/// hundred nodes (contact is a curve or a point); exhausting the budget
/// means genuinely entangled geometry, which falls through to the full
/// solver for an honest verdict.
const CERTIFICATE_NODE_BUDGET: usize = 50_000;

/// A degree-9 tensor-product Bézier control grid (the binonic's `grid`
/// field). The convex-hull property makes its AABB a surface bound, and
/// de Casteljau subdivision at `t = 0.5` is exact.
type BezierGrid = [[Point3; 10]; 10];

/// Try to prove a patch pair crossing-free without the surface-surface
/// solver: recursively subdivide both control grids and certify each
/// sub-box pair as either disjoint (padded by the strict `band`) or
/// close only within `exempt_depth` of the shared contact locus
/// (`center distance + half diagonal` interval bound). Returns `false`
/// — *not certified*, never "crossing" — when the subdivision budget or
/// depth runs out.
pub(super) fn pair_certified_clean(
    grid_a: &BezierGrid,
    grid_b: &BezierGrid,
    shared_spokes: &[SharedSpoke<'_>],
    apex: Point3,
    band: f64,
    exempt_depth: f64,
) -> bool {
    let mut budget = CERTIFICATE_NODE_BUDGET;
    certify_boxes(
        grid_a,
        grid_b,
        shared_spokes,
        apex,
        band,
        exempt_depth,
        0,
        &mut budget,
    )
}

#[allow(clippy::too_many_arguments)]
fn certify_boxes(
    grid_a: &BezierGrid,
    grid_b: &BezierGrid,
    shared_spokes: &[SharedSpoke<'_>],
    apex: Point3,
    band: f64,
    exempt_depth: f64,
    depth: usize,
    budget: &mut usize,
) -> bool {
    let Some(bb_a) = grid_box(grid_a) else {
        return true;
    };
    let Some(bb_b) = grid_box(grid_b) else {
        return true;
    };
    if !bb_a.overlaps(&bb_b, band) {
        return true;
    }
    // Interval bound over the overlap of the padded boxes: any close
    // approach between the two sub-surfaces lies inside it, so if its
    // farthest point from the contact locus is within the exemption
    // depth, this pair region is legitimate contact.
    let overlap = Aabb {
        min: Point3::new(
            (bb_a.min.x - band).max(bb_b.min.x - band),
            (bb_a.min.y - band).max(bb_b.min.y - band),
            (bb_a.min.z - band).max(bb_b.min.z - band),
        ),
        max: Point3::new(
            (bb_a.max.x + band).min(bb_b.max.x + band),
            (bb_a.max.y + band).min(bb_b.max.y + band),
            (bb_a.max.z + band).min(bb_b.max.z + band),
        ),
    };
    let center = Point3::new(
        0.5 * (overlap.min.x + overlap.max.x),
        0.5 * (overlap.min.y + overlap.max.y),
        0.5 * (overlap.min.z + overlap.max.z),
    );
    let half_diagonal = 0.5 * overlap.diagonal();
    let mut locus_distance = center.distance_to(&apex);
    for spoke in shared_spokes {
        locus_distance = locus_distance.min(polyline_distance(center, spoke.polyline));
    }
    if locus_distance + half_diagonal <= exempt_depth {
        return true;
    }
    if depth >= CERTIFICATE_MAX_DEPTH || *budget == 0 {
        return false;
    }
    *budget = budget.saturating_sub(1);
    // Halve both grids along alternating axes and require every child
    // pair to certify.
    let ((a0, a1), (b0, b1)) = if depth.is_multiple_of(2) {
        (split_grid_u(grid_a), split_grid_u(grid_b))
    } else {
        (split_grid_v(grid_a), split_grid_v(grid_b))
    };
    for child_a in [&a0, &a1] {
        for child_b in [&b0, &b1] {
            if !certify_boxes(
                child_a,
                child_b,
                shared_spokes,
                apex,
                band,
                exempt_depth,
                depth + 1,
                budget,
            ) {
                return false;
            }
        }
    }
    true
}

/// AABB of a control grid (a surface bound by the convex-hull property).
fn grid_box(grid: &BezierGrid) -> Option<Aabb> {
    let points: Vec<Point3> = grid.iter().flatten().copied().collect();
    Aabb::from_points(&points)
}

/// Exact de Casteljau split of a degree-9 Bézier grid at `u = 0.5`.
fn split_grid_u(grid: &BezierGrid) -> (BezierGrid, BezierGrid) {
    let mut left = *grid;
    let mut right = *grid;
    for j in 0..10 {
        let mut work: [Point3; 10] = std::array::from_fn(|i| grid[i][j]);
        left[0][j] = work[0];
        right[9][j] = work[9];
        for level in 1..10 {
            for i in 0..(10 - level) {
                work[i] = Point3::new(
                    0.5 * (work[i].x + work[i + 1].x),
                    0.5 * (work[i].y + work[i + 1].y),
                    0.5 * (work[i].z + work[i + 1].z),
                );
            }
            left[level][j] = work[0];
            right[9 - level][j] = work[9 - level];
        }
    }
    (left, right)
}

/// Exact de Casteljau split of a degree-9 Bézier grid at `v = 0.5`.
fn split_grid_v(grid: &BezierGrid) -> (BezierGrid, BezierGrid) {
    let mut left = *grid;
    let mut right = *grid;
    for (i, row) in grid.iter().enumerate() {
        let mut work: [Point3; 10] = *row;
        left[i][0] = work[0];
        right[i][9] = work[9];
        for level in 1..10 {
            for k in 0..(10 - level) {
                work[k] = Point3::new(
                    0.5 * (work[k].x + work[k + 1].x),
                    0.5 * (work[k].y + work[k + 1].y),
                    0.5 * (work[k].z + work[k + 1].z),
                );
            }
            left[i][level] = work[0];
            right[i][9 - level] = work[9 - level];
        }
    }
    (left, right)
}

/// Point-to-segment distance over a polyline (no refinement — the
/// certificate's interval bound already absorbs the chord sagitta at
/// the scales it certifies).
fn polyline_distance(point: Point3, polyline: &[Point3]) -> f64 {
    let mut best_d2 = f64::INFINITY;
    for i in 0..(polyline.len() - 1) {
        let a = polyline[i];
        let ab = polyline[i + 1] - a;
        let ab_len2 = ab.length_squared();
        let s = if ab_len2 > 0.0 {
            ((point - a).dot(&ab) / ab_len2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let foot = Point3::new(a.x + s * ab.x, a.y + s * ab.y, a.z + s * ab.z);
        best_d2 = best_d2.min((point - foot).length_squared());
    }
    best_d2.sqrt()
}

/// Distance from `point` to patch `spoke_patch`'s `u = 1` boundary iso.
///
/// Point-to-segment distance over the precomputed polyline (exact for
/// straight spokes; off by at most the chord sagitta for curved ones —
/// orders below the exemption band at ring scales), then tightened by a
/// fixed-step true-curve refinement around the best segment's parameter
/// (no golden-section, whose descending-range form has bitten before).
/// The polyline term may undershoot the true curve distance by the same
/// sagitta; for an exemption test that direction is conservative.
pub(super) fn distance_to_spoke(
    point: Point3,
    polyline: &[Point3],
    spoke_patch: &RectangularGregoryBiNonic,
) -> f64 {
    let n = polyline.len();
    let step = 1.0 / (n - 1) as f64;
    let mut best_d2 = f64::INFINITY;
    let mut best_t = 0.0_f64;
    for i in 0..(n - 1) {
        let a = polyline[i];
        let ab = polyline[i + 1] - a;
        let ab_len2 = ab.length_squared();
        let s = if ab_len2 > 0.0 {
            ((point - a).dot(&ab) / ab_len2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let foot = Point3::new(a.x + s * ab.x, a.y + s * ab.y, a.z + s * ab.z);
        let d2 = (point - foot).length_squared();
        if d2 < best_d2 {
            best_d2 = d2;
            best_t = (i as f64 + s) * step;
        }
    }
    let mut best = best_d2.sqrt();
    const REFINE: usize = 20;
    let lo = (best_t - step).max(0.0);
    let hi = (best_t + step).min(1.0);
    for s in 0..=REFINE {
        let t = lo + (hi - lo) * (s as f64) / REFINE as f64;
        best = best.min(point.distance_to(&spoke_patch.point_at(1.0, t)));
    }
    best
}
