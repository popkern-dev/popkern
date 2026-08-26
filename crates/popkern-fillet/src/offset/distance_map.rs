//! Distance-map self-fold scan for free-form offsets (PK-327 tier 3,
//! after Seong/Elber/Kim 2006, "Trimming local and global
//! self-intersections in offset curves/surfaces using distance maps").
//!
//! An offset point `O(r,t)` belongs to a self-intersecting region iff
//! some source point `S(u,v)` lies closer than `d − ρ` to it —
//! `F(u,v,r,t) = ‖S(u,v) − O(r,t)‖² − (d − ρ)² ≤ 0`. Matched points
//! `(u,v) = (r,t)` sit at distance exactly `d` and are excluded by the
//! trimming margin `ρ` without any explicit diagonal carve-out; the
//! paper recommends ρ between 1% and 5% of `d`, and popkern fixes
//! **ρ = 2% of d** (`RHO_FRACTION`). Tiny self-intersection loops
//! whose depth stays inside the ρ band escape this scan by
//! construction — the paper's own recommendation routes them to the
//! LOCAL curvature test, which popkern runs as tier 1
//! ([`super::self_intersect::check_offset_self_intersection`]); that
//! division of labor is the tiny-loop guard.
//!
//! The zero-set machinery here is interval subdivision with AABB
//! distance certificates (a lower bound on the true pair distance
//! prunes a box pair whenever it exceeds `d − ρ`), not the paper's
//! rational-spline zero-set solver — same certificates, the
//! subdivision substrate popkern already trusts for SSI seeding.
//! Surviving box pairs at the resolution floor project onto the
//! offset's `rt`-space as fold regions.

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Fixed trimming margin: ρ = 2% of the offset distance.
pub const RHO_FRACTION: f64 = 0.02;

/// Subdivision budget (box pairs examined) — the truncation flag
/// records when the scan gave up rather than silently under-reporting.
const MAX_BOXES: usize = 200_000;

/// Resolution floor as a fraction of the initial parameter span.
const MIN_SPAN_FRACTION: f64 = 1.0 / 32.0;

/// One detected fold region in the OFFSET surface's parameter space.
#[derive(Debug, Clone, Copy)]
pub struct FoldRegion {
    /// Offset `r` (≙ u) window.
    pub r_range: (f64, f64),
    /// Offset `t` (≙ v) window.
    pub t_range: (f64, f64),
}

/// Result of [`scan_offset_self_fold`].
#[derive(Debug, Clone, Default)]
pub struct FoldScanReport {
    /// Fold regions (rt-projections of the surviving box pairs,
    /// merged where overlapping).
    pub regions: Vec<FoldRegion>,
    /// Box pairs examined (observability).
    pub boxes_examined: usize,
    /// True when the subdivision budget was exhausted — the report is
    /// then a lower bound, never a clean bill.
    pub truncated: bool,
}

/// Scan a free-form offset for global self-folding against its source.
///
/// `source` and `offset` share the same parameterization (popkern's
/// offsets are same-parameterized); `distance` is the |offset|.
#[must_use]
pub fn scan_offset_self_fold(
    source: &SurfaceGeometry,
    offset: &SurfaceGeometry,
    distance: f64,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> FoldScanReport {
    let mut report = FoldScanReport::default();
    let threshold = distance.abs() * (1.0 - RHO_FRACTION);
    if threshold <= 0.0 {
        return report;
    }
    let min_span_u = (u_range.1 - u_range.0) * MIN_SPAN_FRACTION;
    let min_span_v = (v_range.1 - v_range.0) * MIN_SPAN_FRACTION;
    if min_span_u <= 0.0 || min_span_v <= 0.0 {
        return report;
    }

    type Box2 = ((f64, f64), (f64, f64));
    let mut stack: Vec<(Box2, Box2)> = vec![((u_range, v_range), (u_range, v_range))];
    let mut hits: Vec<FoldRegion> = Vec::new();

    while let Some((sb, ob)) = stack.pop() {
        if report.boxes_examined >= MAX_BOXES {
            report.truncated = true;
            break;
        }
        report.boxes_examined += 1;
        let source_aabb = sample_aabb(source, sb);
        let offset_aabb = sample_aabb(offset, ob);
        let (Some(sa), Some(oa)) = (source_aabb, offset_aabb) else {
            continue;
        };
        // Certificate: AABB gap is a lower bound on the true pair
        // distance (sampled AABBs are inflated by the box's sampled
        // diagonal to cover the un-sampled interior).
        if aabb_gap(&sa, &oa) > threshold {
            continue;
        }
        let (su, sv) = sb;
        let (ou, ov) = ob;
        // Near-diagonal exemption at fine scale: parameter-adjacent
        // patches sit at pairwise distance ≥ d (offsets are PARALLEL
        // surfaces), and anything below that is a LOCAL fold — tier
        // 1's territory (the tiny-loop guard). Without this the
        // diagonal band subdivides to the floor everywhere.
        let fine = (su.1 - su.0) <= 4.0 * min_span_u
            && (sv.1 - sv.0) <= 4.0 * min_span_v
            && (ou.1 - ou.0) <= 4.0 * min_span_u
            && (ov.1 - ov.0) <= 4.0 * min_span_v;
        if fine
            && ranges_overlap(su, ou, 2.0 * min_span_u)
            && ranges_overlap(sv, ov, 2.0 * min_span_v)
        {
            continue;
        }
        // Early accept at fine granularity: a sampled pair at or below
        // the threshold is a WITNESS of F ≤ 0 (sampled distance upper-
        // bounds the pair minimum). Coarse witnesses keep subdividing
        // so the reported region localizes.
        if fine && sampled_pair_min(source, sb, offset, ob) <= threshold {
            hits.push(FoldRegion {
                r_range: ou,
                t_range: ov,
            });
            continue;
        }
        let at_floor = (su.1 - su.0) <= min_span_u
            && (sv.1 - sv.0) <= min_span_v
            && (ou.1 - ou.0) <= min_span_u
            && (ov.1 - ov.0) <= min_span_v;
        if at_floor {
            // Unproven at the floor: the AABB gap admitted the pair
            // but no sampled witness confirmed — discard (tier-1 owns
            // anything shallower than ρ).
            continue;
        }
        // Split the wider box of the pair (largest parametric span).
        let s_span = (su.1 - su.0).max(sv.1 - sv.0);
        let o_span = (ou.1 - ou.0).max(ov.1 - ov.0);
        if s_span >= o_span {
            for half in split_box(sb) {
                stack.push((half, ob));
            }
        } else {
            for half in split_box(ob) {
                stack.push((sb, half));
            }
        }
    }

    report.regions = merge_regions(hits);
    report
}

fn ranges_overlap((a0, a1): (f64, f64), (b0, b1): (f64, f64), margin: f64) -> bool {
    a0 <= b1 + margin && b0 <= a1 + margin
}

fn split_box(((u0, u1), (v0, v1)): ((f64, f64), (f64, f64))) -> [((f64, f64), (f64, f64)); 2] {
    if (u1 - u0) >= (v1 - v0) {
        let um = (u0 + u1) * 0.5;
        [((u0, um), (v0, v1)), ((um, u1), (v0, v1))]
    } else {
        let vm = (v0 + v1) * 0.5;
        [((u0, u1), (v0, vm)), ((u0, u1), (vm, v1))]
    }
}

/// Sampled AABB over a parameter box, inflated by the largest adjacent
/// sample gap so the un-sampled interior is covered.
fn sample_aabb(
    surface: &SurfaceGeometry,
    ((u0, u1), (v0, v1)): ((f64, f64), (f64, f64)),
) -> Option<Aabb> {
    const N: usize = 4;
    let mut points = Vec::with_capacity((N + 1) * (N + 1));
    for i in 0..=N {
        for j in 0..=N {
            let u = u0 + (u1 - u0) * (i as f64) / N as f64;
            let v = v0 + (v1 - v0) * (j as f64) / N as f64;
            points.push(surface.point_at(u, v));
        }
    }
    let mut max_gap: f64 = 0.0;
    for i in 0..=N {
        for j in 0..=N {
            let p = points[i * (N + 1) + j];
            if i < N {
                max_gap = max_gap.max(p.distance_to(&points[(i + 1) * (N + 1) + j]));
            }
            if j < N {
                max_gap = max_gap.max(p.distance_to(&points[i * (N + 1) + j + 1]));
            }
        }
    }
    let aabb = Aabb::from_points(&points)?;
    let pad = max_gap * 0.5;
    Some(Aabb {
        min: Point3::new(aabb.min.x - pad, aabb.min.y - pad, aabb.min.z - pad),
        max: Point3::new(aabb.max.x + pad, aabb.max.y + pad, aabb.max.z + pad),
    })
}

/// Distance between two AABBs (0 when overlapping).
fn aabb_gap(a: &Aabb, b: &Aabb) -> f64 {
    let axis_gap = |a0: f64, a1: f64, b0: f64, b1: f64| (b0 - a1).max(a0 - b1).max(0.0);
    let dx = axis_gap(a.min.x, a.max.x, b.min.x, b.max.x);
    let dy = axis_gap(a.min.y, a.max.y, b.min.y, b.max.y);
    let dz = axis_gap(a.min.z, a.max.z, b.min.z, b.max.z);
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Minimum sampled distance between two boxed patches (5×5 each).
fn sampled_pair_min(
    source: &SurfaceGeometry,
    ((su0, su1), (sv0, sv1)): ((f64, f64), (f64, f64)),
    offset: &SurfaceGeometry,
    ((ou0, ou1), (ov0, ov1)): ((f64, f64), (f64, f64)),
) -> f64 {
    const N: usize = 2;
    let mut best = f64::INFINITY;
    for i in 0..=N {
        for j in 0..=N {
            let s = source.point_at(
                su0 + (su1 - su0) * (i as f64) / N as f64,
                sv0 + (sv1 - sv0) * (j as f64) / N as f64,
            );
            for k in 0..=N {
                for l in 0..=N {
                    let o = offset.point_at(
                        ou0 + (ou1 - ou0) * (k as f64) / N as f64,
                        ov0 + (ov1 - ov0) * (l as f64) / N as f64,
                    );
                    best = best.min(s.distance_to(&o));
                }
            }
        }
    }
    best
}

/// Merge overlapping rt-boxes (loose union — reports are regions, not
/// exact zero-sets).
fn merge_regions(mut boxes: Vec<FoldRegion>) -> Vec<FoldRegion> {
    let overlaps = |a: &FoldRegion, b: &FoldRegion| {
        a.r_range.0 <= b.r_range.1
            && b.r_range.0 <= a.r_range.1
            && a.t_range.0 <= b.t_range.1
            && b.t_range.0 <= a.t_range.1
    };
    let mut merged: Vec<FoldRegion> = Vec::new();
    while let Some(mut current) = boxes.pop() {
        let mut changed = true;
        while changed {
            changed = false;
            let mut i = 0;
            while i < boxes.len() {
                if overlaps(&current, &boxes[i]) {
                    let other = boxes.swap_remove(i);
                    current.r_range.0 = current.r_range.0.min(other.r_range.0);
                    current.r_range.1 = current.r_range.1.max(other.r_range.1);
                    current.t_range.0 = current.t_range.0.min(other.t_range.0);
                    current.t_range.1 = current.t_range.1.max(other.t_range.1);
                    changed = true;
                } else {
                    i += 1;
                }
            }
        }
        merged.push(current);
    }
    merged.reverse();
    merged
}
