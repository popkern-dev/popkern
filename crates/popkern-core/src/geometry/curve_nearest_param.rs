//! Nearest-parameter curve queries robust to multi-sheet geometry
//! (extracted from `curve_inversion.rs` under the structural line
//! rule). See [`nearest_param_on_curve`].

use crate::geometry::Point3;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::curve::CurveOnEdge;

/// Nearest parameter on a curve-on-edge window to `point`, returned
/// with its distance. Robust to MULTI-SHEET geometry: a coiled sweep
/// seam passes near a query point once per turn, one distance basin
/// each, and a fixed-count coarse scan can rank an adjacent turn's
/// basin above the true one (PK-449, measured: a 16-sample scan on a
/// 6.5-turn helical seam put the plate-crossing 0.209 from the wrong
/// turn's sample vs 0.219 from the right turn's — the ternary bracket
/// then never contained the crossing and a legitimate edge split was
/// classified off-curve).
///
/// The structure is deliberately the classic scan + single-bracket
/// ternary (first strict minimum wins, best EVALUATED station
/// returned): density, not machinery, was the defect — the coarse
/// count adapts to the NURBS knot-span count so a many-turn curve is
/// sampled several times per span. A refine-the-best-k-brackets
/// variant with distance-tie policies was built and FALSIFIED three
/// times over (filleted-torus adoption gates read the bracket-midpoint
/// distance; a cap-cut sphere's rim windows flipped ends): the
/// legacy sites' first-minimum single-bracket semantics are
/// load-bearing across the boolean corpus.
pub fn nearest_param_on_curve(curve: &CurveOnEdge, point: &Point3) -> (f64, f64) {
    use crate::geometry::curve::CurveEval as _;
    let (lo, hi) = (
        curve.t_start.min(curve.t_end),
        curve.t_start.max(curve.t_end),
    );
    let d_at = |t: f64| curve.geometry.point_at(t).distance_to(point);
    if !(hi - lo).is_finite() || hi - lo <= 0.0 {
        return (curve.t_start, d_at(curve.t_start));
    }
    // Coarse density: 64 base, at least 2 samples per knot span in
    // the window for NURBS (capped).
    let n: usize = match &curve.geometry {
        CurveGeometry::Nurbs(nurbs) => {
            let mut spans = 0usize;
            let mut last = f64::NEG_INFINITY;
            for &k in &nurbs.knots {
                if k >= lo && k <= hi && k > last {
                    spans += 1;
                    last = k;
                }
            }
            (2 * spans).clamp(64, 4096)
        }
        _ => 64,
    };
    let mut best = (0usize, f64::MAX);
    for i in 0..=n {
        let t = lo + (hi - lo) * i as f64 / n as f64;
        let d = d_at(t);
        if d < best.1 {
            best = (i, d);
        }
    }
    let (mut a, mut b) = (
        lo + (hi - lo) * (best.0.saturating_sub(1)) as f64 / n as f64,
        lo + (hi - lo) * (best.0 + 1).min(n) as f64 / n as f64,
    );
    let mut best_t = lo + (hi - lo) * best.0 as f64 / n as f64;
    let mut best_d = best.1;
    for _ in 0..48 {
        let m1 = a + (b - a) / 3.0;
        let m2 = b - (b - a) / 3.0;
        let d1 = d_at(m1);
        let d2 = d_at(m2);
        if d1 < d2 {
            b = m2;
            if d1 < best_d {
                best_d = d1;
                best_t = m1;
            }
        } else {
            a = m1;
            if d2 < best_d {
                best_d = d2;
                best_t = m2;
            }
        }
    }
    (best_t, best_d)
}
