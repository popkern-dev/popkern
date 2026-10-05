//! Boundary-ring UV conditioning for the Delaunay incircle test
//! (extracted from `boundary_first.rs` under the structural line
//! budget): collinear-run inward nudge and deterministic cocircular
//! micro-jitter. Both perturb UV only — 3D positions stay canonical,
//! so shared-edge closure is unaffected.

use super::BoundaryRing;
use popkern_math::PortableMath;

/// Nudge strictly-collinear boundary points INWARD by an invisible
/// epsilon in UV. Collinear boundary runs (straight edges sampled at
/// several stations) sit exactly on one line in UV; Bowyer–Watson
/// then connects the run's endpoints directly and the interior
/// stations never join the hull — the coverage check would reject
/// the whole face.
pub(super) fn nudge_collinear_runs(rings: &mut [BoundaryRing]) {
    for ring in rings.iter_mut() {
        let n = ring.uv.len();
        if n < 3 {
            continue;
        }
        let mut scale = 0.0_f64;
        for i in 0..n {
            let a = ring.uv[i];
            let b = ring.uv[(i + 1) % n];
            scale = scale.max((a.0 - b.0).abs().max((a.1 - b.1).abs()));
        }
        let eps = scale.max(1e-9) * 1.0e-7;
        let mut area2 = 0.0;
        for i in 0..n {
            let (x1, y1) = ring.uv[i];
            let (x2, y2) = ring.uv[(i + 1) % n];
            area2 += x1 * y2 - x2 * y1;
        }
        let inward = if area2 > 0.0 { 1.0 } else { -1.0 };
        let orig: Vec<(f64, f64)> = ring.uv.clone();
        for i in 0..n {
            let p = orig[i];
            let q = orig[(i + n - 1) % n];
            let r = orig[(i + 1) % n];
            let cross = (p.0 - q.0) * (r.1 - q.1) - (p.1 - q.1) * (r.0 - q.0);
            let span = ((r.0 - q.0).powi(2) + (r.1 - q.1).powi(2))
                .sqrt()
                .max(1e-300);
            if (cross / span).abs() < eps {
                // Left normal of q→r is inward for CCW rings.
                let (dx, dy) = ((r.0 - q.0) / span, (r.1 - q.1) / span);
                ring.uv[i] = (p.0 - -(dy * eps * inward), p.1 + -(dx * eps * inward));
            }
        }
    }
}

/// Deterministic micro-jitter on every boundary point's UV (≈1e-7 of
/// the ring scale, direction hashed from the index): exact COCIRCULAR
/// point groups — font-glyph arcs flatten into them by construction —
/// read as degenerate to the incircle test the same way
/// exactly-collinear runs do, and cascade coverage holes through the
/// Delaunay. The jitter is orders of magnitude below any chord
/// tolerance.
pub(super) fn jitter_cocircular(rings: &mut [BoundaryRing]) {
    let mut scale = 0.0_f64;
    for ring in rings.iter() {
        let n = ring.uv.len();
        for i in 0..n {
            let a = ring.uv[i];
            let b = ring.uv[(i + 1) % n];
            scale = scale.max((a.0 - b.0).abs().max((a.1 - b.1).abs()));
        }
    }
    let eps = scale.max(1e-9) * 1.0e-7;
    let mut k = 0u64;
    for ring in rings.iter_mut() {
        for uv in &mut ring.uv {
            // Cheap deterministic hash → angle.
            k = k
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let angle = (k >> 11) as f64 / (1u64 << 53) as f64 * std::f64::consts::TAU;
            uv.0 += eps * angle.pcos();
            uv.1 += eps * angle.psin();
        }
    }
}
