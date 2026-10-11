//! Slit dissolution for boundary-first rings (PK-335).
//!
//! A slit face (the shell rim wall over a closed removed-face
//! boundary is the canonical case: an annulus whose outer loop runs
//! outer circle → radial slit → inner circle → radial slit) traverses
//! its slit edges TWICE within one loop, and the two traversals are
//! COINCIDENT in the chart: the ring polygon has zero-width arms.
//! Delaunay coverage can never cover a zero-width segment pair, and
//! the encroachment escalation that follows splits the slit spans
//! forever — worse, those futile splits enrich the canonical samples
//! of edges SHARED with neighbouring faces, un-equalizing their
//! ruled-strip chains (measured on the shelled cylinder: 7 rounds of
//! splits, then uniform-grid fallback for the rim and Delaunay fallback for
//! the walls).
//!
//! The slit is purely a B-rep bookkeeping artifact for the MESH: both
//! half-edges of a slit edge lie on this same face, so no other face
//! borders it and the tessellation needs no boundary there at all.
//! Dissolution cuts the ring at the coincident runs, yielding the
//! genuine boundary components (outer rim + hole), and Delaunay
//! covers them directly.
//!
//! A periodic SEAM also traverses one edge twice in a loop. Its uv
//! images differ by the chart period AFTER branch unwrapping — but
//! dissolution runs on the raw ring, where branch-incoherent seam
//! pcurves can coincide exactly like a slit (measured: the shelled
//! cylinder's wall seams). Dissolution is therefore scoped to charts
//! with NO periodic axis (the shell-rim class is planar); a slit on
//! a periodic chart keeps the pre-existing escalation path — a
//! recorded fallback boundary, not silent incompleteness.

use super::BoundaryRing;
use popkern_math::PortableMath;

/// Split `ring` at every coincident double-traversal (slit) into its
/// genuine boundary components. Rings without slits — and rings on
/// periodic charts, where raw-branch seam coincidence is
/// indistinguishable from a slit — pass through unchanged.
pub(super) fn dissolve_slits(
    ring: BoundaryRing,
    periodicity: (Option<f64>, Option<f64>),
) -> Vec<BoundaryRing> {
    if periodicity.0.is_some() || periodicity.1.is_some() {
        return vec![ring];
    }
    let mut out = Vec::new();
    let mut work = vec![ring];
    // Bounded: each dissolution strictly shrinks total point count.
    while let Some(ring) = work.pop() {
        match find_and_cut_slit(&ring) {
            Some(pieces) => work.extend(pieces),
            None => out.push(ring),
        }
    }
    // Deterministic output order: largest component first, then by
    // first uv (the work stack permutes pieces).
    out.sort_by(|a, b| {
        b.uv.len().cmp(&a.uv.len()).then_with(|| {
            a.uv.first()
                .map(|p| (p.0, p.1))
                .partial_cmp(&b.uv.first().map(|p| (p.0, p.1)))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    });
    out
}

fn find_and_cut_slit(ring: &BoundaryRing) -> Option<Vec<BoundaryRing>> {
    let n = ring.uv.len();
    if n < 6 || ring.sources.len() != n {
        return None;
    }
    // Chart-scale epsilons: uv coincidence must hold in the chart (a
    // seam's runs differ by a full period there); positions are
    // canonical samples of the same edge, so they match near-exactly.
    let (mut u0, mut u1, mut v0, mut v1) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &(u, v) in &ring.uv {
        u0 = u0.min(u);
        u1 = u1.max(u);
        v0 = v0.min(v);
        v1 = v1.max(v);
    }
    let eps_uv = ((u1 - u0).phypot(v1 - v0)).max(1e-12) * 1e-6;
    let scale = ring
        .pos
        .iter()
        .zip(ring.pos.iter().cycle().skip(1))
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0_f64, f64::max);
    let eps_pos = (1.0 + scale) * 1e-9;

    // Maximal same-edge segment runs around the ring: (start, len).
    let mut runs: Vec<(usize, usize)> = Vec::new();
    let mut i = 0;
    while i < n {
        let eh = ring.sources[i].0;
        let mut len = 1;
        while len < n && ring.sources[(i + len) % n].0 == eh {
            len += 1;
        }
        runs.push((i, len));
        i += len;
    }
    // Merge a wrap-around run (last and first share the edge).
    if runs.len() >= 2 {
        let (fs, fl) = runs[0];
        let (ls, ll) = *runs.last()?;
        if ring.sources[fs].0 == ring.sources[ls].0 && ls + ll == n && fs == 0 {
            runs[0] = (ls, ll + fl);
            runs.pop();
        }
    }

    for a in 0..runs.len() {
        for b in (a + 1)..runs.len() {
            let (s1, k) = runs[a];
            let (s2, k2) = runs[b];
            if ring.sources[s1].0 != ring.sources[s2].0 || k != k2 {
                continue;
            }
            // Reversed coincidence over the k+1 run points.
            let coincident = (0..=k).all(|t| {
                let p = (s1 + t) % n;
                let q = (s2 + k - t) % n;
                let (pu, pv) = ring.uv[p];
                let (qu, qv) = ring.uv[q];
                (pu - qu).phypot(pv - qv) <= eps_uv
                    && ring.pos[p].distance_to(&ring.pos[q]) <= eps_pos
            });
            if !coincident {
                continue;
            }
            let piece = |from: usize, to: usize| -> Option<BoundaryRing> {
                let mut uv = Vec::new();
                let mut pos = Vec::new();
                let mut sources = Vec::new();
                let mut i = from;
                loop {
                    uv.push(ring.uv[i]);
                    pos.push(ring.pos[i]);
                    sources.push(ring.sources[i]);
                    if i == to {
                        break;
                    }
                    i = (i + 1) % n;
                }
                // A closed component repeats its junction point.
                if uv.len() >= 2
                    && pos
                        .first()
                        .zip(pos.last())
                        .map(|(x, y)| x.distance_to(y) <= eps_pos)
                        .unwrap_or(false)
                {
                    uv.pop();
                    pos.pop();
                    sources.pop();
                }
                (uv.len() >= 3).then_some(BoundaryRing { uv, pos, sources })
            };
            let piece_a = piece((s1 + k) % n, s2);
            let piece_b = piece((s2 + k) % n, s1);
            let pieces: Vec<BoundaryRing> = [piece_a, piece_b].into_iter().flatten().collect();
            if pieces.is_empty() {
                return None; // nothing meshable would remain — keep the ring
            }
            if popkern_core::trace::trace_enabled("boundary-first") {
                eprintln!(
                    "    [bf] slit dissolved: edge {:?} double-traversal cut into {} component(s)",
                    ring.sources[s1].0,
                    pieces.len()
                );
            }
            return Some(pieces);
        }
    }
    None
}

/// Split `ring` at every PINCH — a chart point the ring visits twice
/// (same UV within the chart epsilon and the same 3-D position) —
/// into the simple polygons that touch there.
///
/// A band around a doubly-periodic chart whose loop traverses the
/// v-seam twice is such a ring: a torus cut by a plane keeps the
/// bottom band, and its loop runs one half of the band on the base
/// u-sheet and the other half one period over, the two halves
/// meeting only at the equator/tube-seam corner (2π, 2π). The chain
/// is continuous, so the unwrap cannot fold it; as one polygon it is
/// a bow-tie (measured: Delaunay coverage failed at the pinch, the
/// encroachment rounds split the equator forever, and the fallback
/// grid left 16 open edges). Each lobe is a plain rectangle on its
/// own sheet; both consume the canonical equator samples, so they
/// weld. Runs after the unwrap (the pinch is a property of the
/// unwrapped chart) and on every chart.
pub(super) fn dissolve_pinches(ring: BoundaryRing) -> Vec<BoundaryRing> {
    dissolve_pinches_among(ring, pinch_candidates)
}

/// The pairs of ring points [`find_and_cut_pinch`] tests, given the
/// ring and its uv coincidence epsilon.
type PinchCandidates = fn(&BoundaryRing, f64) -> Vec<(usize, usize)>;

fn dissolve_pinches_among(ring: BoundaryRing, candidates: PinchCandidates) -> Vec<BoundaryRing> {
    let mut out = Vec::new();
    let mut work = vec![ring];
    while let Some(ring) = work.pop() {
        match find_and_cut_pinch(&ring, candidates) {
            Some(pieces) => work.extend(pieces),
            None => out.push(ring),
        }
    }
    out.sort_by(|a, b| {
        b.uv.len().cmp(&a.uv.len()).then_with(|| {
            a.uv.first()
                .map(|p| (p.0, p.1))
                .partial_cmp(&b.uv.first().map(|p| (p.0, p.1)))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    });
    out
}

fn find_and_cut_pinch(
    ring: &BoundaryRing,
    candidates: PinchCandidates,
) -> Option<Vec<BoundaryRing>> {
    let n = ring.uv.len();
    if n < 6 || ring.sources.len() != n {
        return None;
    }
    let (mut u0, mut u1, mut v0, mut v1) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &(u, v) in &ring.uv {
        u0 = u0.min(u);
        u1 = u1.max(u);
        v0 = v0.min(v);
        v1 = v1.max(v);
    }
    let eps_uv = ((u1 - u0).phypot(v1 - v0)).max(1e-12) * 1e-6;
    let scale = ring
        .pos
        .iter()
        .zip(ring.pos.iter().cycle().skip(1))
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0_f64, f64::max);
    let eps_pos = (1.0 + scale) * 1e-9;
    for (i, j) in candidates(ring, eps_uv) {
        if i == 0 && j + 3 > n {
            continue; // too close around the wrap to bound a lobe
        }
        let (pu, pv) = ring.uv[i];
        let (qu, qv) = ring.uv[j];
        if (pu - qu).phypot(pv - qv) > eps_uv || ring.pos[i].distance_to(&ring.pos[j]) > eps_pos {
            continue;
        }
        let take = |indices: Vec<usize>| BoundaryRing {
            uv: indices.iter().map(|&k| ring.uv[k]).collect(),
            pos: indices.iter().map(|&k| ring.pos[k]).collect(),
            sources: indices.iter().map(|&k| ring.sources[k]).collect(),
        };
        let lobe_a = take((i..j).collect());
        let lobe_b = take((j..n).chain(0..i).collect());
        if lobe_a.uv.len() < 3 || lobe_b.uv.len() < 3 {
            continue;
        }
        if popkern_core::trace::trace_enabled("boundary-first") {
            eprintln!(
                "    [bf] pinch dissolved at uv ({pu:.4},{pv:.4}): ring cut into {} + {} points",
                lobe_a.uv.len(),
                lobe_b.uv.len()
            );
        }
        return Some(vec![lobe_a, lobe_b]);
    }
    None
}

/// Every pair `(i, j)`, `j >= i + 3`, that could be a pinch, in
/// lexicographic order: the pairs whose uv points share or neighbour a
/// grid cell of side `eps_uv`, which covers every pair within `eps_uv`.
/// The caller applies the exact pinch test to each, so the first pair
/// it accepts is the one an all-pairs scan would accept. That scan was
/// O(n^2) per ring and took 0.9 s on one 80 mm thread wall's
/// 15.6k-sample ring (PK-794). A ring with a non-finite coordinate
/// keeps the all-pairs order: NaN passes the scan's distance test.
fn pinch_candidates(ring: &BoundaryRing, eps_uv: f64) -> Vec<(usize, usize)> {
    let finite = ring.uv.iter().all(|&(u, v)| u.is_finite() && v.is_finite())
        && ring
            .pos
            .iter()
            .all(|p| p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
    if !finite {
        return all_pinch_pairs(ring, eps_uv);
    }
    let cell = |(u, v): (f64, f64)| ((u / eps_uv).floor() as i64, (v / eps_uv).floor() as i64);
    let mut grid: std::collections::HashMap<(i64, i64), Vec<usize>> =
        std::collections::HashMap::new();
    for (k, &uv) in ring.uv.iter().enumerate() {
        grid.entry(cell(uv)).or_default().push(k);
    }
    let mut pairs = Vec::new();
    for (i, &uv) in ring.uv.iter().enumerate() {
        let (cu, cv) = cell(uv);
        for du in -1..=1 {
            for dv in -1..=1 {
                if let Some(members) = grid.get(&(cu.saturating_add(du), cv.saturating_add(dv))) {
                    pairs.extend(members.iter().filter(|&&j| j >= i + 3).map(|&j| (i, j)));
                }
            }
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

/// Every pair `(i, j)`, `j >= i + 3`, in lexicographic order.
fn all_pinch_pairs(ring: &BoundaryRing, _eps_uv: f64) -> Vec<(usize, usize)> {
    let n = ring.uv.len();
    (0..n)
        .flat_map(|i| ((i + 3)..n).map(move |j| (i, j)))
        .collect()
}

#[cfg(test)]
#[path = "boundary_first_slit_tests.rs"]
mod tests;
