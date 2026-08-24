//! Periodic-branch unwrapping for boundary-first rings (split out
//! of `boundary_first.rs` for the 500-line ratchet).
//!
//! Pcurves on periodic charts can be authored in different period
//! branches; the unwrap walks the ring and shifts each sample to the
//! branch nearest its predecessor (positions are canonical 3-D points
//! and unaffected). Chart-degenerate samples (poles/apexes, arbitrary
//! azimuth) are skipped as anchors and rewritten from their flanking
//! regular samples. Returns `false` on a genuine net full wrap — the
//! caller routes those to the structured-revolution or legacy paths.

use super::BoundaryRing;

pub(super) fn unwrap_dir(
    ring: &mut BoundaryRing,
    period: f64,
    pick: fn(&mut (f64, f64)) -> &mut f64,
    degenerate: &[bool],
) -> bool {
    let is_degen = |i: usize| degenerate.get(i).copied().unwrap_or(false);
    let mut prev: Option<f64> = None;
    for i in 0..ring.uv.len() {
        if is_degen(i) {
            continue;
        }
        let cur = pick(&mut ring.uv[i]);
        if let Some(p) = prev {
            while *cur - p > 0.5 * period {
                *cur -= period;
            }
            while *cur - p < -0.5 * period {
                *cur += period;
            }
        }
        prev = Some(*cur);
    }
    // Rewrite each degenerate run's coordinate to the
    // mean of its flanking regular samples (circular):
    // any value evaluates to the same 3D point, and the
    // mean keeps both adjacent UV segments non-zero for
    // the Delaunay coverage walk.
    let n = ring.uv.len();
    for i in 0..n {
        if !is_degen(i) {
            continue;
        }
        let before = (1..n)
            .map(|k| (i + n - k) % n)
            .find(|&j| !is_degen(j))
            .map(|j| *pick(&mut ring.uv[j].clone()));
        let after = (1..n)
            .map(|k| (i + k) % n)
            .find(|&j| !is_degen(j))
            .map(|j| *pick(&mut ring.uv[j].clone()));
        if let (Some(b), Some(a)) = (before, after) {
            *pick(&mut ring.uv[i]) = 0.5 * (b + a);
        }
    }
    // Net winding across the regular samples = a genuine
    // full wrap; boundary-first cannot chart it.
    let regular: Vec<usize> = (0..n).filter(|&i| !is_degen(i)).collect();
    match (regular.first(), regular.last()) {
        (Some(&f), Some(&l)) => {
            let first = *pick(&mut ring.uv[f].clone());
            let last = *pick(&mut ring.uv[l].clone());
            (last - first).abs() <= 0.5 * period
        }
        _ => true,
    }
}
// On decline, dump the unwrapped ring samples: the jump
// pattern distinguishes genuine winding (monotone drift
// over the ring) from branch-incoherent pcurves (one
// edge's run offset by a full period) — they need
// different fixes.
pub(super) fn dump_ring(ring: &BoundaryRing, dir: &str) {
    if std::env::var_os("PK_BF_TRACE").is_some() {
        eprintln!("    [bf] full {dir}-wrap ring — cannot chart");
        for (i, ((u, v), p)) in ring.uv.iter().zip(ring.pos.iter()).enumerate() {
            eprintln!(
                "      [bf] ring[{i}] uv=({u:.4},{v:.4}) pos=({:.4},{:.4},{:.4})",
                p.x, p.y, p.z
            );
        }
    }
}
