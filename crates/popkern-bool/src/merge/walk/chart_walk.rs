//! Permutation walk over a face's chart subdivision, and the dual-walk
//! comparator (PK-104 stages 2–3 — see
//! docs/architecture/fragment-walk-chart-subdivision.md §3.2–3.3).
//!
//! `next()` is the classical DCEL rotational successor for
//! material-on-left traversal: among half-edges departing the branch
//! point where `h` arrives, take the one with the smallest strictly
//! positive CLOCKWISE angle from `h`'s backward direction in oriented
//! UV, twin ranked last. It is a pure function of the half-edge, so
//! loops are the orbits of a permutation — the partition cannot depend
//! on seed order.

use std::collections::BTreeSet;

use super::chart_subdivision::FaceChartSubdivision;
use crate::merge::graph::IntersectionGraph;

/// The rotational successor of `h` on this face's chart.
///
/// Returns `None` only on substrate defects (missing branch point or
/// no departures) — impossible on a subdivision that passed I1–I5,
/// and treated as a face-level fallback trigger by callers.
pub(crate) fn next_in_chart(
    sub: &FaceChartSubdivision,
    graph: &IntersectionGraph,
    h: usize,
) -> Option<usize> {
    let bp = sub.end_bp.get(&h)?;
    let deps = sub.departures.get(bp)?;
    let twin = graph.half_edges[h].twin;

    // POLE branch points rotate by AZIMUTH: every incident half-edge
    // lands at the pole carrying its meridian's u, and the geometric
    // angle around the pole IS that u — local UV directions there all
    // point along ±v and cannot order the fan. The material-left
    // successor is the departure with the smallest positive
    // (u_arrive − u_depart)·dir mod period, dir folding σ with the
    // pole side (rotation sense reverses between poles); the twin
    // sits at the arrival azimuth (angle → period) and ranks last
    // naturally.
    if let Some(&dir) = sub.poles.get(&bp.node)
        && let Some(period) = sub.u_period
    {
        let emb_h = sub.embeddings.get(&h)?;
        let u_a = emb_h.uv[emb_h.uv.len() - 1].0;
        let mut best: Option<(f64, usize)> = None;
        for &c in deps {
            let Some(emb_c) = sub.embeddings.get(&c) else {
                continue;
            };
            let u_c = emb_c.uv[0].0;
            let mut a = ((u_a - u_c) * dir).rem_euclid(period);
            if a <= 1e-12 || c == twin {
                a += period;
            }
            let better = match best {
                None => true,
                Some((ba, bc)) => a < ba - 1e-15 || ((a - ba).abs() <= 1e-15 && c < bc),
            };
            if better {
                best = Some((a, c));
            }
        }
        return best.map(|(_, c)| c);
    }

    let arrive = sub.embeddings.get(&h)?.end_dir()?;
    let back = (-arrive.0, -arrive.1);

    let mut best: Option<(f64, usize)> = None;
    for &c in deps {
        let angle = if c == twin {
            // U-turn along the twin: strictly last. Still eligible so
            // a genuine spur walks back rather than dead-ending.
            std::f64::consts::TAU + 1.0
        } else {
            let Some(dep) = sub.embeddings.get(&c).and_then(|e| e.start_dir()) else {
                continue;
            };
            // Clockwise angle from `back` to `dep` in oriented UV:
            // positive cross = CCW in raw UV; sigma flips handedness
            // for orientation-false faces; negation converts to CW.
            let cross = back.0 * dep.1 - back.1 * dep.0;
            let dot = back.0 * dep.0 + back.1 * dep.1;
            let a = (-sub.sigma * cross).atan2(dot);
            if a <= 1e-12 {
                a + std::f64::consts::TAU
            } else {
                a
            }
        };
        let better = match best {
            None => true,
            Some((ba, bc)) => angle < ba - 1e-15 || ((angle - ba).abs() <= 1e-15 && c < bc),
        };
        if better {
            best = Some((angle, c));
        }
    }
    best.map(|(_, c)| c)
}

/// Enumerate the permutation's orbits, ascending from the lowest
/// unconsumed half-edge (deterministic, seed-order independent by
/// construction). Returns `None` when any `next` fails or an orbit
/// fails to close within the face's half-edge budget.
pub(crate) fn walk_chart_orbits(
    sub: &FaceChartSubdivision,
    graph: &IntersectionGraph,
) -> Option<Vec<Vec<usize>>> {
    let mut consumed: BTreeSet<usize> = BTreeSet::new();
    let mut orbits = Vec::new();
    let budget = sub.embeddings.len() + 1;
    for &seed in sub.embeddings.keys() {
        if consumed.contains(&seed) {
            continue;
        }
        let mut orbit = vec![seed];
        consumed.insert(seed);
        let mut h = seed;
        loop {
            let n = next_in_chart(sub, graph, h)?;
            if n == seed {
                break;
            }
            // A permutation orbit cannot revisit a non-seed member; if
            // it does, the substrate lied — bail to fallback.
            if consumed.contains(&n) || orbit.len() > budget {
                return None;
            }
            orbit.push(n);
            consumed.insert(n);
            h = n;
        }
        orbits.push(orbit);
    }
    Some(orbits)
}

/// Canonical form for loop-set comparison: each loop rotated so its
/// smallest half-edge index leads (loops are cyclic; direction is
/// preserved because direction is part of the contract). Test-only
/// since the dual-walk probe retired (PK-108); the convention pins
/// use it to compare partitions.
#[cfg(test)]
pub(crate) fn canonical_loops(loops: &[Vec<usize>]) -> BTreeSet<Vec<usize>> {
    loops
        .iter()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let min_pos = l
                .iter()
                .enumerate()
                .min_by_key(|&(_, &v)| v)
                .map(|(i, _)| i)
                .unwrap_or(0);
            let mut rotated = Vec::with_capacity(l.len());
            rotated.extend_from_slice(&l[min_pos..]);
            rotated.extend_from_slice(&l[..min_pos]);
            rotated
        })
        .collect()
}

/// Divergence between two walk partitions of one face (test-only
/// since the dual-walk probe retired — PK-108).
#[cfg(test)]
#[derive(Debug)]
pub(crate) struct DualWalkDiff {
    /// Loops present in both partitions (canonical form).
    pub matched: usize,
    /// Legacy loops the orbit partition does not contain.
    pub legacy_only: Vec<Vec<usize>>,
    /// Orbits the legacy walk never emitted. Orbits consuming only
    /// half-edges legacy left unwalked are benign extras (legacy
    /// skips some cross-registered seeds); orbits sharing half-edges
    /// with differing legacy loops are true routing divergences.
    pub orbit_only: Vec<Vec<usize>>,
}

#[cfg(test)]
impl DualWalkDiff {
    pub fn is_match(&self) -> bool {
        self.legacy_only.is_empty() && self.orbit_only.is_empty()
    }

    /// Orbit-only loops that share half-edges with legacy-only loops —
    /// the genuinely divergent routing (as opposed to benign extras).
    pub fn contested(&self) -> usize {
        let legacy_hes: BTreeSet<usize> = self.legacy_only.iter().flatten().copied().collect();
        self.orbit_only
            .iter()
            .filter(|o| o.iter().any(|h| legacy_hes.contains(h)))
            .count()
    }
}

/// Compare two walk partitions (test-only since PK-108).
#[cfg(test)]
pub(crate) fn compare_walks(legacy: &[Vec<usize>], orbits: &[Vec<usize>]) -> DualWalkDiff {
    let l = canonical_loops(legacy);
    let o = canonical_loops(orbits);
    let matched = l.intersection(&o).count();
    DualWalkDiff {
        matched,
        legacy_only: l.difference(&o).cloned().collect(),
        orbit_only: o.difference(&l).cloned().collect(),
    }
}
