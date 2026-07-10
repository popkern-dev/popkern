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
    let arrive = sub.embeddings.get(&h)?.end_dir()?;
    let back = (-arrive.0, -arrive.1);
    let twin = graph.half_edges[h].twin;

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
/// preserved because direction is part of the contract).
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

/// Divergence between the legacy walk's loops and the chart orbits
/// for one face.
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

/// Compare the legacy walk's loop set against the chart orbits.
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

/// Verbose substrate dump: one line per embedded half-edge with its
/// chart endpoints, branch points, directions, and source — the data
/// needed to adjudicate an incomplete or collision case by hand.
fn dump_substrate(sub: &FaceChartSubdivision, graph: &IntersectionGraph) {
    eprintln!("    [dual-sub] sigma={}", sub.sigma);
    for (&h, emb) in &sub.embeddings {
        let (u0, v0) = emb.uv[0];
        let (u1, v1) = emb.uv[emb.uv.len() - 1];
        let sd = emb.start_dir().unwrap_or((0.0, 0.0));
        let ed = emb.end_dir().unwrap_or((0.0, 0.0));
        let us: Vec<String> = emb.uv.iter().map(|&(u, _)| format!("{u:.2}")).collect();
        eprintln!(
            "    [dual-sub] he {h} tw {} {:?} ({u0:.4},{v0:.4})->({u1:.4},{v1:.4}) sdir=({:.3},{:.3}) edir=({:.3},{:.3}) sbp={:?} ebp={:?} seam={} u=[{}]",
            graph.half_edges[h].twin,
            emb.source,
            sd.0,
            sd.1,
            ed.0,
            ed.1,
            sub.start_bp.get(&h),
            sub.end_bp.get(&h),
            emb.seam_like,
            us.join(",")
        );
    }
}

/// Stage-2 dual-walk instrument (env-gated by the caller): build the
/// subdivision, walk the orbits, and report agreement with the legacy
/// walk for this face. Output is one line per face so corpus sweeps
/// aggregate with grep; PK_WALK_DUAL_VERBOSE adds the diverging loops.
pub(crate) fn dual_walk_probe(
    graph: &IntersectionGraph,
    face: &crate::merge::graph::FaceRef,
    he_indices: &[usize],
    legacy_loops: &[Vec<usize>],
    brep_a: &popkern_core::topology::BRep,
    brep_b: &popkern_core::topology::BRep,
) {
    use super::chart_subdivision::{build_face_subdivision, face_surface};
    let Some((surface, orientation)) = face_surface(face, brep_a, brep_b) else {
        eprintln!("[dual] {face:?}: no-surface");
        return;
    };
    let sub = build_face_subdivision(graph, face, he_indices, surface, orientation);
    let verbose = std::env::var_os("PK_WALK_DUAL_VERBOSE").is_some();
    if !sub.complete() {
        let issues: Vec<String> = sub.report.issues.iter().map(|i| i.to_string()).collect();
        eprintln!("[dual] {face:?}: incomplete: {}", issues.join("; "));
        if verbose {
            dump_substrate(&sub, graph);
        }
        return;
    }
    match walk_chart_orbits(&sub, graph) {
        None => {
            eprintln!("[dual] {face:?}: orbit-walk-failed");
            if verbose {
                dump_substrate(&sub, graph);
            }
        }
        Some(orbits) => {
            let diff = compare_walks(legacy_loops, &orbits);
            if diff.is_match() {
                eprintln!("[dual] {face:?}: MATCH {}", diff.matched);
                if std::env::var_os("PK_WALK_DUAL_DUMP_ALL").is_some() {
                    dump_substrate(&sub, graph);
                }
            } else {
                eprintln!(
                    "[dual] {face:?}: DIVERGE matched={} legacy_only={} orbit_only={} contested={}",
                    diff.matched,
                    diff.legacy_only.len(),
                    diff.orbit_only.len(),
                    diff.contested()
                );
                if std::env::var_os("PK_WALK_DUAL_VERBOSE").is_some() {
                    for l in &diff.legacy_only {
                        eprintln!("    [dual-legacy] {l:?}");
                    }
                    for o in &diff.orbit_only {
                        eprintln!("    [dual-orbit]  {o:?}");
                    }
                }
            }
        }
    }
}
