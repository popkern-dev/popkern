//! Knot-seeded initial fractions with alias-proof coarsening for
//! canonical edge sampling (PK-494) — split from
//! `boundary_first_edge_samples.rs` under the 500-line budget. The
//! measured rationale lives on the code below, verbatim from the
//! parent.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveGeometry;

/// Seed `fractions` at every distinct interior knot (budgeted), then
/// coarsen alias-proof: merge seeded intervals back while every
/// DISCARDED seed point stays within `chord_height` of the merged
/// chord.
pub(super) fn knot_seed_and_coarsen(
    geometry: &CurveGeometry,
    t_start: f64,
    t_end: f64,
    chord_height: f64,
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
                    let ok =
                        (s_idx + 1..candidate).all(|k| dev(s_idx, candidate, k) <= chord_height);
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
