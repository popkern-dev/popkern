//! Exact even-odd containment for loops on a chart that is periodic
//! in `u` and open in `v` (cylinders, cones, closed-profile
//! extrusions).
//!
//! A region on such a chart is the projection of its loops' lifts.
//! A ray cast from a point along +v never wraps, and it crosses the
//! region's boundary an odd number of times exactly when the point
//! lies inside. The count includes crossings with EVERY period lift
//! of every loop; outer and inner loops count alike, so holes
//! subtract by construction. Consecutive samples are joined by their
//! SHORT representative (`|Δu| ≤ period/2`). A run authored on
//! another period branch then contributes the same cylinder segment
//! as its aligned copy, and a branch jump contributes nothing. The
//! aliased-seam convention is handled the same way: a helical band's
//! turn pieces each end at u = 0 and restart at u = period, and those
//! jumps also contribute nothing.
//!
//! PK-566 measured two consumers whose local tests failed where this
//! is exact:
//! - The cellular periodic-fold interior sampler read "material" from
//!   the side of the nearest boundary sample. That vote ties at a seam
//!   edge pair (one 3D segment walked both ways). On a nut-mouth
//!   chamfer cone both threaded-band fragments sampled each other's
//!   region, and their verdicts swapped.
//! - The ray classifier's multi-turn band gate ran a horizontal
//!   even-odd test on the FOLDED trim polygon. That polygon has no
//!   seam segments, so a crest strip that wraps the seam at the query
//!   height reads as its complement. A bore-void ray took a phantom
//!   crest-band hit, and a void point on the M6 nut classified Inside.

/// Whether `(u, v)` lies inside the region bounded by `loops` (outer
/// plus inner loops, each a closed UV polyline in walk order) on a
/// chart periodic in `u` with `period` and open in `v`.
///
/// Crossings use the half-open convention in `u`, so a ray through a
/// sample vertex counts once. A ray lying exactly along a vertical
/// segment (a seam edge) contributes no crossing from it; callers
/// sample at arc midpoints, away from loop crossings.
pub(crate) fn contains_on_u_periodic_chart(
    loops: &[&[(f64, f64)]],
    u: f64,
    v: f64,
    period: f64,
) -> bool {
    contains_on_u_periodic_chart_with_poles(loops, u, v, period, &|_| false)
}

/// [`contains_on_u_periodic_chart`] on a chart with collapsed rows:
/// `pole_row(v)` names the rows (a sphere's poles, a cone's apex) that
/// are one 3D point. A step of whole periods along such a row is the
/// loop walking it, a real boundary in the chart, and counts at its
/// full length: a sphere face's outer loop runs both seam copies and
/// joins them along its pole rows, and with those rows read as branch
/// jumps the face's parity inverted and its holes read as material
/// (PK-937).
pub(crate) fn contains_on_u_periodic_chart_with_poles(
    loops: &[&[(f64, f64)]],
    u: f64,
    v: f64,
    period: f64,
    pole_row: &dyn Fn(f64) -> bool,
) -> bool {
    let mut crossings = 0usize;
    for poly in loops {
        let n = poly.len();
        for i in 0..n {
            let (u0, v0) = poly[i];
            let (u1_raw, v1) = poly[(i + 1) % n];
            let raw = u1_raw - u0;
            let wrapped = raw - period * (raw / period).round();
            let along_pole = wrapped.abs() <= 1e-9 * period
                && raw.abs() > 0.5 * period
                && (v1 - v0).abs() <= 1e-9 * (1.0 + v0.abs())
                && pole_row(v0);
            let du = if along_pole { raw } else { wrapped };
            let u1 = u0 + du;
            if du == 0.0 {
                continue;
            }
            // The lifts of the ray's u within reach of the segment: a
            // short step (|du| ≤ period/2) straddles at most one among
            // the start's nearest and its neighbours; a pole row spans
            // its whole length.
            let base = u + period * ((u0 - u) / period).round();
            let reach = (du.abs() / period).ceil().max(1.0) as i64;
            for k in -reach..=reach {
                let x = base + k as f64 * period;
                if (u0 > x) != (u1 > x) {
                    let t = (x - u0) / du;
                    if v0 + t * (v1 - v0) > v {
                        crossings += 1;
                    }
                }
            }
        }
    }
    crossings % 2 == 1
}

/// Sort folded crossings and cancel coincident PAIRS: two crossings
/// at one folded position toggle containment twice (the seam copies
/// u = 0 and u = period of one seam edge, or a tangential touch).
/// Keeping one of the pair would leave an odd count and a spurious
/// toggle.
pub(crate) fn cancel_coincident_pairs(mut xs: Vec<f64>, eps: f64) -> Vec<f64> {
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut out: Vec<f64> = Vec::with_capacity(xs.len());
    for x in xs {
        match out.last() {
            Some(&last) if (x - last).abs() < eps => {
                out.pop();
            }
            _ => out.push(x),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    /// The PK-566 cone band shape: a void band between an apex ring
    /// (v = 1) and an S-shaped cut curve, walked with seam edges on
    /// both sides of the strip.
    fn s_band() -> Vec<(f64, f64)> {
        let mut poly = vec![(0.0, 1.0), (TAU, 1.0), (TAU, 2.0)];
        // S-curve from u = 2π back to 0, doubling back in u.
        for &(u, v) in &[(4.0, 1.5), (3.6, 1.5), (5.0, 3.0), (2.7, 3.0), (0.0, 2.0)] {
            poly.push((u, v));
        }
        poly
    }

    #[test]
    fn s_band_containment_follows_the_band_not_the_seam() {
        let band = s_band();
        let loops: [&[(f64, f64)]; 1] = [&band];
        // Just above the apex ring, on both sides of the seam: inside.
        assert!(contains_on_u_periodic_chart(&loops, 0.2, 1.1, TAU));
        assert!(contains_on_u_periodic_chart(&loops, TAU - 0.2, 1.1, TAU));
        // At u = 4.5 the vertical line meets three strands: the
        // descending strand (v ≈ 1.61), the doubled-back strand
        // (v ≈ 2.46), and the upper strand (v = 3.0). The region
        // alternates: in, out, in.
        assert!(contains_on_u_periodic_chart(&loops, 4.5, 1.2, TAU));
        assert!(!contains_on_u_periodic_chart(&loops, 4.5, 2.0, TAU));
        assert!(contains_on_u_periodic_chart(&loops, 4.5, 2.9, TAU));
        // Lifting every sample by a period changes nothing.
        let lifted: Vec<(f64, f64)> = band.iter().map(|&(u, v)| (u + TAU, v)).collect();
        let lifted_loops: [&[(f64, f64)]; 1] = [&lifted];
        for (u, v) in [(0.2, 1.1), (4.5, 2.9), (4.5, 2.0), (4.5, 1.2), (1.0, 2.5)] {
            assert_eq!(
                contains_on_u_periodic_chart(&loops, u, v, TAU),
                contains_on_u_periodic_chart(&lifted_loops, u, v, TAU)
            );
        }
    }

    #[test]
    fn ring_band_with_a_branch_jump_counts_rings() {
        // A band between two full rings, the upper ring authored on
        // the next period branch halfway (a branch-glue jump).
        let lower: Vec<(f64, f64)> = (0..=16).map(|k| (TAU * k as f64 / 16.0, 1.0)).collect();
        let upper: Vec<(f64, f64)> = (0..=16)
            .rev()
            .map(|k| {
                let u = TAU * k as f64 / 16.0;
                (if k < 8 { u + TAU } else { u }, 2.0)
            })
            .collect();
        let loops: [&[(f64, f64)]; 2] = [&lower, &upper];
        assert!(contains_on_u_periodic_chart(&loops, 1.3, 1.5, TAU));
        assert!(!contains_on_u_periodic_chart(&loops, 1.3, 2.5, TAU));
        assert!(!contains_on_u_periodic_chart(&loops, 1.3, 0.5, TAU));
    }

    /// PK-937: a sphere face's outer loop as the chart walk embeds it,
    /// up the u = 0 seam copy, along the north pole row, down the
    /// u = 2π copy and back along the south row, holed at u = 1. Read
    /// as branch jumps the pole rows bound nothing and the hole reads
    /// as the material; read as the rows they are, the region is the
    /// chart less the hole, on every branch.
    #[test]
    fn pole_rows_bound_a_sphere_face() {
        use std::f64::consts::PI;
        let outer = vec![(0.0, PI), (0.0, 0.0), (TAU, 0.0), (TAU, PI)];
        let hole: Vec<(f64, f64)> = (0..16)
            .map(|k| {
                let t = -TAU * k as f64 / 16.0;
                (1.0 + 0.3 * t.cos(), 1.5 + 0.3 * t.sin())
            })
            .collect();
        let loops: [&[(f64, f64)]; 2] = [&outer, &hole];
        let poles = |v: f64| v == 0.0 || v == PI;
        for lift in [-TAU, 0.0, TAU] {
            let inside = |u: f64, v: f64| {
                contains_on_u_periodic_chart_with_poles(&loops, u + lift, v, TAU, &poles)
            };
            assert!(inside(3.0, 1.5) && inside(1.0, 0.5) && inside(5.9, 2.9));
            assert!(!inside(1.0, 1.5), "the hole is not material");
        }
        // Without the pole rows the parity inverts: the reading this
        // replaces.
        assert!(contains_on_u_periodic_chart(&loops, 1.0, 1.5, TAU));
    }

    #[test]
    fn coincident_crossings_cancel_in_pairs() {
        assert_eq!(
            cancel_coincident_pairs(vec![0.0, 4.0, -1.0, 0.0], 1e-9),
            vec![-1.0, 4.0]
        );
        assert_eq!(
            cancel_coincident_pairs(vec![0.0, 0.0, 0.0], 1e-9),
            vec![0.0]
        );
    }
}
