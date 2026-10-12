//! General surface-curve intersection fallback.
//!
//! Orchestrates: overlap detection -> curve Bezier decomposition -> surface
//! bounding-box grid (8x8) -> AABB culling -> Newton 3x3 polish -> dedup -> sort.

use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::intersection::SurfaceCurveIntersection;
use crate::geometry::intersection::surface_curve_newton::newton_surface_curve;
use crate::geometry::intersection::surface_curve_overlap::detect_surface_curve_overlaps;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

#[path = "surface_curve_grid.rs"]
mod surface_curve_grid;
use surface_curve_grid::{
    GRID_N, SurfaceCell, exact_conic_segment_box, interior_knots_in, padded_cell_memo,
    surface_grid_cached,
};
#[path = "surface_curve_bezier.rs"]
mod surface_curve_bezier;
pub(super) use surface_curve_bezier::decompose_to_bezier_segments;
#[path = "surface_curve_pair_roots.rs"]
mod surface_curve_pair_roots;
use surface_curve_pair_roots::{
    PairCell, Root, curve_period, isolate_further_roots, lift_into_window,
};

// ---------------------------------------------------------------------------
// Curve AABB
// ---------------------------------------------------------------------------

struct CurveSegment {
    t_range: (f64, f64),
    min: Point3,
    max: Point3,
}

/// Compute AABB for a NURBS curve by sampling a set of points.
fn curve_aabb(nurbs: &NurbsCurve, t_range: (f64, f64)) -> (Point3, Point3) {
    let n_samples = 16;
    let (t0, t1) = t_range;
    let mut p = nurbs.point_at(t0);
    let mut bmin = p;
    let mut bmax = p;
    for i in 1..=n_samples {
        let frac = i as f64 / n_samples as f64;
        p = nurbs.point_at(t0 + frac * (t1 - t0));
        if p.x < bmin.x {
            bmin.x = p.x;
        }
        if p.y < bmin.y {
            bmin.y = p.y;
        }
        if p.z < bmin.z {
            bmin.z = p.z;
        }
        if p.x > bmax.x {
            bmax.x = p.x;
        }
        if p.y > bmax.y {
            bmax.y = p.y;
        }
        if p.z > bmax.z {
            bmax.z = p.z;
        }
    }
    (bmin, bmax)
}

/// Check AABB overlap with a tolerance margin.
/// Whether the straight segment `a`→`b` meets the box `[lo, hi]`
/// grown by `margin` (slab test). Exact for a segment, where the
/// segment's own AABB is a gross over-approximation: a diagonal
/// classification ray spans a box that overlaps EVERY cell of a
/// thread face, and each overlap paid span refinement and Newton
/// seeding (PK-513).
fn segment_meets_box(a: &Point3, b: &Point3, lo: &Point3, hi: &Point3, margin: f64) -> bool {
    let (mut t0, mut t1) = (0.0_f64, 1.0_f64);
    for (p, d, l, h) in [
        (a.x, b.x - a.x, lo.x - margin, hi.x + margin),
        (a.y, b.y - a.y, lo.y - margin, hi.y + margin),
        (a.z, b.z - a.z, lo.z - margin, hi.z + margin),
    ] {
        if d.abs() < 1e-300 {
            if p < l || p > h {
                return false;
            }
            continue;
        }
        let (mut ta, mut tb) = ((l - p) / d, (h - p) / d);
        if ta > tb {
            std::mem::swap(&mut ta, &mut tb);
        }
        t0 = t0.max(ta);
        t1 = t1.min(tb);
        if t0 > t1 {
            return false;
        }
    }
    true
}

fn aabb_overlap(
    a_min: &Point3,
    a_max: &Point3,
    b_min: &Point3,
    b_max: &Point3,
    margin: f64,
) -> bool {
    a_min.x - margin <= b_max.x
        && a_max.x + margin >= b_min.x
        && a_min.y - margin <= b_max.y
        && a_max.y + margin >= b_min.y
        && a_min.z - margin <= b_max.z
        && a_max.z + margin >= b_min.z
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// General surface-curve intersection solver.
///
/// Orchestrates: overlap detection, curve Bezier decomposition, surface
/// bounding-box grid, AABB culling, Newton 3x3 polish, dedup, sort.
#[allow(clippy::too_many_arguments)]
pub(crate) fn intersect_surface_curve_general(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    curve: &CurveGeometry,
    t_range: (f64, f64),
    tolerance: &Tolerance,
    dense_grid: bool,
) -> Vec<SurfaceCurveIntersection> {
    // --- Step 1: overlap detection ---
    let mut results =
        detect_surface_curve_overlaps(surface, u_range, v_range, curve, t_range, tolerance);

    // --- Step 2: promote curve to NURBS and decompose into Bezier segments ---
    let nurbs = curve.to_nurbs_in_range(t_range);
    let bezier_segments = decompose_to_bezier_segments(&nurbs);

    // Compute AABB for each Bezier segment.
    let curve_segments: Vec<CurveSegment> = bezier_segments
        .into_iter()
        .map(|(seg, tr)| {
            let (bmin, bmax) =
                exact_conic_segment_box(curve, &seg, tr).unwrap_or_else(|| curve_aabb(&seg, tr));
            CurveSegment {
                t_range: tr,
                min: bmin,
                max: bmax,
            }
        })
        .collect();

    // --- Step 3: build surface bbox grid ---
    let surface_cells = surface_grid_cached(surface, u_range, v_range, GRID_N, dense_grid);

    // --- Step 4: AABB culling + Newton polish ---
    let margin = tolerance.linear() * 10.0; // small margin for AABB overlap test
    let mut point_results: Vec<(f64, f64, f64, Point3)> = Vec::new();
    let dedup_dist = tolerance.linear() * 10.0;
    let t_period = curve_period(curve);
    let (u_period, v_period) = surface.periodicity();

    // Per-span refinement for hit cells: the density clamp keeps the
    // coarse grid bounded (a 198-span thread wall grids at 48), so a
    // hit cell may span several knot spans — and a Newton seed at its
    // midpoint sits spans away from the crossing on a wrapped chart
    // and diverges (the PK-441 class, re-measured on the PK-460
    // thread coil at 4-span cells). Hit cells refine to single-span
    // sub-cells (padded boxes, re-culled) before seeding; only the
    // ray's corridor pays.
    let fingerprint = match surface {
        SurfaceGeometry::Nurbs(s) => Some(s.content_fingerprint()),
        _ => None,
    };
    let (knots_u, knots_v): (&[f64], &[f64]) = match surface {
        SurfaceGeometry::Nurbs(s) => (&s.knots_u, &s.knots_v),
        _ => (&[], &[]),
    };
    // Straight curves (classification rays) cull by the exact segment
    // instead of its AABB.
    let straight = matches!(curve, CurveGeometry::Line { .. });
    let straight_meets = |cseg: &CurveSegment, lo: &Point3, hi: &Point3| -> bool {
        !straight || {
            let a = curve.point_at(cseg.t_range.0);
            let b = curve.point_at(cseg.t_range.1);
            segment_meets_box(&a, &b, lo, hi, margin)
        }
    };
    for cseg in &curve_segments {
        for cell in surface_cells.iter() {
            if !aabb_overlap(&cseg.min, &cseg.max, &cell.min, &cell.max, margin)
                || !straight_meets(cseg, &cell.min, &cell.max)
            {
                if crate::trace::trace_enabled("scg-cull") {
                    eprintln!(
                        "  [scg-cull] u({:.4},{:.4}) v({:.4},{:.4}) box ({:.2},{:.2},{:.2})..({:.2},{:.2},{:.2}) vs seg ({:.2},{:.2},{:.2})..({:.2},{:.2},{:.2})",
                        cell.u_range.0,
                        cell.u_range.1,
                        cell.v_range.0,
                        cell.v_range.1,
                        cell.min.x,
                        cell.min.y,
                        cell.min.z,
                        cell.max.x,
                        cell.max.y,
                        cell.max.z,
                        cseg.min.x,
                        cseg.min.y,
                        cseg.min.z,
                        cseg.max.x,
                        cseg.max.y,
                        cseg.max.z
                    );
                }
                continue;
            }
            let mut u_edges = vec![cell.u_range.0];
            u_edges.extend(interior_knots_in(knots_u, cell.u_range.0, cell.u_range.1));
            u_edges.push(cell.u_range.1);
            let mut v_edges = vec![cell.v_range.0];
            v_edges.extend(interior_knots_in(knots_v, cell.v_range.0, cell.v_range.1));
            v_edges.push(cell.v_range.1);
            let refined = u_edges.len() > 2 || v_edges.len() > 2;
            for ui in 0..u_edges.len() - 1 {
                for vi in 0..v_edges.len() - 1 {
                    let sub = if refined {
                        padded_cell_memo(
                            surface,
                            fingerprint,
                            (u_edges[ui], u_edges[ui + 1]),
                            (v_edges[vi], v_edges[vi + 1]),
                        )
                    } else {
                        SurfaceCell {
                            u_range: cell.u_range,
                            v_range: cell.v_range,
                            min: cell.min,
                            max: cell.max,
                        }
                    };
                    if refined
                        && (!aabb_overlap(&cseg.min, &cseg.max, &sub.min, &sub.max, margin)
                            || !straight_meets(cseg, &sub.min, &sub.max))
                    {
                        continue;
                    }
                    let (sub_u, sub_v) = (sub.u_range, sub.v_range);

                    // Seed Newton at the curve sample NEAREST the
                    // sub-cell, not the segment's t-midpoint: a
                    // classification RAY arrives as one Bezier
                    // segment spanning t ∈ (0, 1e6), so the midpoint
                    // seed sat at t = 5e5 — half a million units from
                    // every real crossing — and Newton either
                    // diverged or converged onto the EXTENDED surface
                    // basis (u ≈ −304) and was range-rejected: rays
                    // through thread coils counted zero wall
                    // crossings and interior points classified
                    // Outside (PK-460). Nine samples per segment keep
                    // the seed within an eighth of the segment.
                    let u_mid = 0.5 * (sub_u.0 + sub_u.1);
                    let v_mid = 0.5 * (sub_v.0 + sub_v.1);
                    let cell_centre = Point3::new(
                        0.5 * (sub.min.x + sub.max.x),
                        0.5 * (sub.min.y + sub.max.y),
                        0.5 * (sub.min.z + sub.max.z),
                    );
                    // Seed selection is scoped: proper curve segments
                    // keep the tuned midpoint seed (nearest-sample
                    // seeding shifts Newton basins and flipped the
                    // revolved-sphere cap-rim fixture's crossings —
                    // measured). Only a midpoint LANDING FAR from the
                    // cell — the classification-RAY case, one Bezier
                    // segment spanning t ∈ (0, 1e6) whose midpoint
                    // sits half a million units out — switches to the
                    // curve sample nearest the sub-cell.
                    let t_plain = 0.5 * (cseg.t_range.0 + cseg.t_range.1);
                    let cell_diag = sub.min.distance_to(&sub.max).max(tolerance.linear());
                    let t_nearest = (0..=8)
                        .map(|k| {
                            cseg.t_range.0 + (cseg.t_range.1 - cseg.t_range.0) * k as f64 / 8.0
                        })
                        .min_by(|&a, &b| {
                            let da = curve.point_at(a).distance_to(&cell_centre);
                            let db = curve.point_at(b).distance_to(&cell_centre);
                            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                        })
                        .unwrap_or(t_plain);
                    let t_mid =
                        if curve.point_at(t_plain).distance_to(&cell_centre) <= 4.0 * cell_diag {
                            t_plain
                        } else {
                            t_nearest
                        };

                    let mut solved =
                        newton_surface_curve(curve, t_mid, surface, u_mid, v_mid, tolerance);
                    // Seed rescue for a FAILED cell only: the tuned seed
                    // keeps its basin (changing it flips fixtures —
                    // PK-449), but a cell whose AABB overlaps the curve
                    // and whose Newton diverged from the segment
                    // t-midpoint deserves a second attempt from the
                    // curve sample nearest the cell — the root, if any,
                    // is IN the cell, so the near sample is the
                    // better-informed seed (measured: the angled-tap
                    // mouth rim crossing its thread flank at t 0.537 in
                    // the t∈[0.5,0.75] Bezier segment — the midpoint
                    // seed 2 arc-units away diverged on the shallow
                    // crossing and the rim×flank junction vanished from
                    // the 1D pierce inventory, stranding the groove-mouth
                    // record chain, PK-465).
                    if solved.is_none() && t_nearest != t_mid {
                        solved = newton_surface_curve(
                            curve, t_nearest, surface, u_mid, v_mid, tolerance,
                        );
                    }
                    if crate::trace::trace_enabled("scg-newton") {
                        eprintln!(
                            "  [scg-newton] seed u={u_mid:.4} v={v_mid:.4} t={t_mid:.3} -> {:?}",
                            solved.as_ref().map(|(t, u, v, _)| (*t, *u, *v))
                        );
                    }
                    // A pair with a root stays ambiguous while a piece
                    // of the segment no root explains still meets the
                    // cell (PK-781): isolate its further roots.
                    let mut pair_roots: Vec<Root> = solved.into_iter().collect();
                    let pair = PairCell {
                        curve,
                        surface,
                        lo: sub.min,
                        hi: sub.max,
                        u_mid,
                        v_mid,
                        margin,
                        tolerance,
                    };
                    let first = solved.map(|r| r.0);
                    if first.is_some() {
                        isolate_further_roots(
                            &pair,
                            cseg.t_range,
                            first,
                            dedup_dist,
                            &mut pair_roots,
                        );
                    }
                    for (t, u, v, pt) in pair_roots {
                        // Within the requested ranges; periodic
                        // parameters lift by whole periods into them.
                        let pb = tolerance.parametric_base();
                        let lifted = lift_into_window(t, t_range, pb, t_period).and_then(|t| {
                            let u = lift_into_window(u, u_range, pb, u_period)?;
                            Some((t, u, lift_into_window(v, v_range, pb, v_period)?))
                        });
                        if let Some((t, u, v)) = lifted {
                            point_results.push((t, u, v, pt));
                        }
                    }
                }
            }
        }
    }

    // --- Step 5: deduplicate ---
    let mut unique: Vec<(f64, f64, f64, Point3)> = Vec::new();
    'outer: for (t, u, v, pt) in point_results {
        for &(et, _eu, _ev, ref ep) in &unique {
            if (t - et).abs() < dedup_dist || pt.distance_to(ep) < dedup_dist {
                continue 'outer;
            }
        }
        unique.push((t, u, v, pt));
    }

    // A tangential touch arrives as roots spread over the tolerance
    // band around it: one contact, one root.
    let unique =
        super::surface_curve_touch::consolidate_touches(curve, surface, unique, t_range, tolerance);

    // Also skip points that fall inside any overlap range.
    for (t, u, v, pt) in unique {
        let in_overlap = results.iter().any(|r| {
            if let SurfaceCurveIntersection::Overlap { t_range: (ta, tb) } = r {
                t >= *ta && t <= *tb
            } else {
                false
            }
        });
        if !in_overlap {
            results.push(SurfaceCurveIntersection::Point { point: pt, t, u, v });
        }
    }

    // --- Step 6: sort by t ---
    results.sort_by(|a, b| {
        let ta = match a {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            SurfaceCurveIntersection::Overlap { t_range } => t_range.0,
        };
        let tb = match b {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            SurfaceCurveIntersection::Overlap { t_range } => t_range.0,
        };
        ta.partial_cmp(&tb).unwrap_or(std::cmp::Ordering::Equal)
    });

    results
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "surface_curve_general_tests.rs"]
mod tests;
