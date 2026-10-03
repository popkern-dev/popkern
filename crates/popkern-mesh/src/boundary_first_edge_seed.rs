//! Knot-seeded initial fractions with alias-proof coarsening for
//! canonical edge sampling (PK-494) — split from
//! `boundary_first_edge_samples.rs` under the 500-line budget. The
//! measured rationale lives on the code below, verbatim from the
//! parent.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveGeometry;

/// Turn measure across ordered interval samples `(fractions, points)`
/// (first and last are the interval ends), in radians.
pub(super) type TurnOver<'a> = dyn Fn(&[f64], &[Point3]) -> f64 + 'a;

/// Seed `fractions` at every distinct interior knot (budgeted), then
/// coarsen alias-proof: merge seeded intervals back while every
/// DISCARDED seed point stays within `chord_height` of the merged
/// chord and (with an angular bound `(max, turn_over)`) the turn
/// measured over the merged interval's seeds stays within `max`.
/// `turn_over(fractions, points)` measures the turn across its
/// ordered samples (first and last are the interval ends).
pub(super) fn knot_seed_and_coarsen(
    geometry: &CurveGeometry,
    t_start: f64,
    t_end: f64,
    chord_height: f64,
    max_turn: Option<(f64, &TurnOver<'_>)>,
    at: &dyn Fn(f64) -> Point3,
    fractions: &mut Vec<f64>,
) {
    // KNOT-SEEDED initial fractions (PK-494): starting from the
    // bare endpoints, the breadth-first sagitta refinement can
    // ALIAS on periodic curves — an interval spanning a whole
    // number of coil turns has its parametric midpoint back in
    // phase with the chord, the deviation test reads ~0, and the
    // interval never splits (measured on a 41-turn crease edge:
    // one surviving segment spanned exactly 20.0 turns — 322
    // metric units of a 776-long chart — and its diametral
    // protecting disk evacuated the ENTIRE interior lattice of
    // the root-band face). Seeding a fraction at every distinct
    // interior knot bounds every initial interval by one knot
    // span (a quarter turn on the canonical helix), where the
    // sagitta test cannot alias. Non-NURBS curves keep the bare
    // endpoints (their sagitta is unimodal per interval).
    if let CurveGeometry::Nurbs(nc) = geometry {
        let (t0, t1) = (t_start, t_end);
        let span = (t1 - t0).abs().max(1e-300);
        let mut seeds: Vec<f64> = nc
            .knots
            .iter()
            .map(|&k| (k - t0) / (t1 - t0))
            .filter(|&f| f > 1e-9 && f < 1.0 - 1e-9)
            .collect();
        seeds.dedup_by(|a, b| (*a - *b).abs() * span < 1e-12 * span.max(1.0));
        // Respect the sample budget: subsample evenly if a very
        // fine knot vector alone would exceed it.
        let budget = 4096usize;
        if seeds.len() > budget {
            let step = seeds.len() as f64 / budget as f64;
            seeds = (0..budget)
                .map(|i| seeds[(i as f64 * step) as usize])
                .collect();
        }
        fractions.extend(seeds);
        fractions.sort_by(f64::total_cmp);
        fractions.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        // Alias-proof coarsening: merge seeded intervals back
        // while every DISCARDED seed point stays within the chord
        // tolerance of the merged chord. Judged on the real seed
        // points — not a parametric midpoint — so a merged
        // interval can never hide interior excursions the way an
        // in-phase midpoint does, and a loose tolerance coarsens
        // back to the historical minimal sample sets (the
        // 2-triangle loose-sphere floor).
        if fractions.len() > 2 {
            let pts: Vec<Point3> = fractions.iter().map(|&f| at(f)).collect();
            let dev = |a: usize, b: usize, k: usize| -> f64 {
                let (pa, pb, pk) = (pts[a], pts[b], pts[k]);
                let (dx, dy, dz) = (pb.x - pa.x, pb.y - pa.y, pb.z - pa.z);
                let l2 = dx * dx + dy * dy + dz * dz;
                if l2 < 1e-300 {
                    return pk.distance_to(&pa);
                }
                let t = (((pk.x - pa.x) * dx + (pk.y - pa.y) * dy + (pk.z - pa.z) * dz) / l2)
                    .clamp(0.0, 1.0);
                let q = Point3::new(pa.x + dx * t, pa.y + dy * t, pa.z + dz * t);
                pk.distance_to(&q)
            };
            let mut kept: Vec<usize> = vec![0];
            let mut s_idx = 0usize;
            let mut j = 1usize;
            while j < fractions.len() {
                let mut end = j;
                while end + 1 < fractions.len() {
                    let candidate = end + 1;
                    let ok = (s_idx + 1..candidate)
                        .all(|k| dev(s_idx, candidate, k) <= chord_height)
                        && max_turn.is_none_or(|(m, turn_over)| {
                            turn_over(&fractions[s_idx..=candidate], &pts[s_idx..=candidate]) <= m
                        });
                    if ok {
                        end = candidate;
                    } else {
                        break;
                    }
                }
                kept.push(end);
                s_idx = end;
                j = end + 1;
            }
            *fractions = kept.iter().map(|&i| fractions[i]).collect();
        }
    }
}

/// Seed a CLOSED edge (a full circle bounding a face on its own) with
/// quarter turns. Its coincident ends start as one zero-chord interval
/// that neither refinement test can split when the curve is smaller
/// than the chord height: the midpoint sags by the diameter, which is
/// below the chord height, and the end normals are the same normal.
/// Left at two samples it is no ring at all, and the face it bounds
/// falls off the boundary-first path (measured: clevis21.stp's
/// 4e-5-radius drill-point circles — six cone faces declined and
/// meshed open). Quarter turns are what bisection reaches on any
/// larger closed curve, so those sample identically.
pub(super) fn seed_closed_curve(at: &dyn Fn(f64) -> Point3, fractions: &mut Vec<f64>) {
    if fractions.len() != 2 {
        return;
    }
    let (start, end) = (at(0.0), at(1.0));
    if start.distance_to(&end) <= 1.0e-12 * (1.0 + at(0.5).distance_to(&start)) {
        *fractions = vec![0.0, 0.25, 0.5, 0.75, 1.0];
    }
}
