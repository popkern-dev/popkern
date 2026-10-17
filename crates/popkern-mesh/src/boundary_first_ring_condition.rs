//! Boundary-ring UV conditioning for the Delaunay incircle test
//! (extracted from `boundary_first.rs` under the structural line
//! budget): collinear-run inward nudge and deterministic cocircular
//! micro-jitter. Both perturb UV only — 3D positions stay canonical,
//! so shared-edge closure is unaffected — and only for the
//! triangulation: the served vertices take back their exact UV
//! (`restore_exact_uv`).

use super::BoundaryRing;
use popkern_math::PortableMath;

/// Nudge strictly-collinear boundary points INWARD by an invisible
/// epsilon in UV. Collinear boundary runs (straight edges sampled at
/// several stations) sit exactly on one line in UV; Bowyer–Watson
/// then connects the run's endpoints directly and the interior
/// stations never join the hull — the coverage check would reject
/// the whole face.
fn nudge_collinear_runs(rings: &mut [BoundaryRing]) {
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
fn jitter_cocircular(rings: &mut [BoundaryRing]) {
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

/// Condition the rings for the incircle test (collinear nudge,
/// cocircular jitter, pinch copies re-snapped together), returning
/// their UV from before, in point order (the order the face's point
/// set lists them) — for [`restore_exact_uv`].
pub(super) fn condition_for_incircle(rings: &mut [BoundaryRing]) -> Vec<(f64, f64)> {
    let shared = shared_points(rings);
    let exact = rings.iter().flat_map(|r| r.uv.iter().copied()).collect();
    nudge_collinear_runs(rings);
    jitter_cocircular(rings);
    resnap_shared_points(rings, &shared);
    exact
}

/// Give the boundary points (the point set's first `exact.len()`) back
/// their unconditioned UV once the triangulation no longer needs the
/// conditioning (PK-877). The nudge and the jitter move a boundary UV
/// off its trim by up to 1e-7 of the ring scale while its 3-D position
/// stays on the trim, and every consumer of the served mesh's UV read
/// that offset as geometry: the exact-surface volume quadrature
/// integrated the (u, v) triangles of a drilled dome's sphere zone
/// 1.3e-8 rad off its exact latitude rim, 1.1e-8 of volume (the face's
/// slivers measured chord-to-trim, not the conditioning).
pub(super) fn restore_exact_uv(uv_pts: &mut [(f64, f64)], exact: &[(f64, f64)]) {
    for (uv, &e) in uv_pts.iter_mut().zip(exact) {
        *uv = e;
    }
}

/// Ring points that coincide exactly — a pinch, where a lobe cut from a
/// ring (`dissolve_pinches`) shares its point with the ring it left —
/// as `(ring, index)` groups, recorded before the conditioning passes.
fn shared_points(rings: &[BoundaryRing]) -> Vec<Vec<(usize, usize)>> {
    let mut by_point: std::collections::BTreeMap<(u64, u64), Vec<(usize, usize)>> =
        std::collections::BTreeMap::new();
    for (r, ring) in rings.iter().enumerate() {
        for (i, uv) in ring.uv.iter().enumerate() {
            by_point
                .entry((uv.0.to_bits(), uv.1.to_bits()))
                .or_default()
                .push((r, i));
        }
    }
    by_point.into_values().filter(|g| g.len() > 1).collect()
}

/// Give every group of [`shared_points`] its first member's conditioned
/// position: the nudge and the jitter move each ring point on its own,
/// and two copies of a pinch moved apart let the lobes' first segments
/// cross, which no constrained triangulation can recover (measured,
/// PK-732: the equal-radius tee's main wall, the copies 3.6e-6 rad
/// apart at the singular point; the wall fell back to a mesh off the
/// canonical edge samples, 893 open mesh edges).
fn resnap_shared_points(rings: &mut [BoundaryRing], groups: &[Vec<(usize, usize)>]) {
    for group in groups {
        let (r0, i0) = group[0];
        let uv = rings[r0].uv[i0];
        for &(r, i) in &group[1..] {
            rings[r].uv[i] = uv;
        }
    }
}
