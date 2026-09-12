//! The rotational-successor step of the chart walk (split from
//! `chart_walk.rs` under the 500-line file budget): DCEL
//! material-on-left successor selection with endpoint-refined wedge
//! directions and the near-tie twin-match repair.

use super::super::chart_subdivision::FaceChartSubdivision;
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

    // GLOBAL WEDGE PAIRING (permutation by construction): the former
    // per-arrival min-CW rule picked each arrival's successor
    // independently, and at a slop-cluster junction two arrivals
    // could elect the SAME departure — the "orbits" then stop being a
    // permutation, the walk collides, and the face falls back to the
    // 3D-Newell walker mis-partitioned (measured: rotated torus × torus,
    // a node where the correct matching requires one arrival to
    // U-turn; the mis-partition cost 5.5% of the union volume before
    // the strand gate declined it). The DCEL-correct assignment is
    // the cyclic wedge matching over ALL slots at the branch point:
    // arrivals (at their back directions) and departures (at their
    // start directions) sort into one σ-oriented circular order and
    // match parenthesis-style sweeping clockwise — each departure
    // pairs the nearest pending arrival, injective by construction.
    let bp_key = *bp;
    let psi = |d: (f64, f64)| -> f64 {
        (sub.sigma * d.1)
            .atan2(d.0)
            .rem_euclid(std::f64::consts::TAU)
    };
    // Slots: (angle, is_departure, he). Tie order in the clockwise
    // sweep puts the DEPARTURE first so an exactly-parallel departure
    // is NOT taken by its parallel arrival (the historic strictly-
    // positive-angle rule; the twin's departure is exactly parallel
    // to its arrival's back and ranks last this way).
    let mut slots: Vec<(f64, bool, usize)> = Vec::new();
    for (hh, e) in sub.end_bp.iter() {
        if *e == bp_key
            && let Some(d) = sub.embeddings.get(hh).and_then(|em| em.end_dir())
        {
            slots.push((psi((-d.0, -d.1)), false, *hh));
        }
    }
    for &c in deps {
        if let Some(d) = sub.embeddings.get(&c).and_then(|em| em.start_dir()) {
            slots.push((psi(d), true, c));
        }
    }
    let _ = twin;
    // TWIN-SLOT ANGLE SNAP: an arrival's back direction and its own
    // twin's departure direction are the same geometric ray, but they
    // are measured from two independently built embeddings (each side
    // of a single-face residue pair builds its own), and the
    // endpoint-refined δ-samples can disagree below the acceptance
    // cone — the ψ tie then breaks EITHER way. When the arrival
    // sorts a hair before its twin's departure, the clockwise sweep
    // pairs the arrival with its own twin and the chain U-turns
    // (measured, PK-465: both residue pairs at one junction node of
    // the angled tap's coil flank hairpinned, the plate-top ring
    // never closed, and the flank's above-plate part classified with
    // its in-material part — 55 masked ghost strands). Snap the
    // arrival's angle onto the twin departure's angle whenever they
    // are within the near-tie cone; the departure-first tie rule
    // then guarantees the departure sweeps first and is only taken
    // by its own arrival when the parenthesis structure demands a
    // genuine U-turn.
    let near = |a: f64, b: f64| -> bool {
        let d = (a - b).abs();
        d.min(std::f64::consts::TAU - d) < 0.02
    };
    for i in 0..slots.len() {
        let (a_ang, is_dep, hh) = slots[i];
        if is_dep {
            continue;
        }
        let t = graph.half_edges[hh].twin;
        if let Some(&(d_ang, _, _)) = slots.iter().find(|(_, d, h2)| *d && *h2 == t)
            && near(a_ang, d_ang)
        {
            slots[i].0 = d_ang;
        }
    }
    // Clockwise = descending ψ. Deterministic ties: departures first,
    // then half-edge index.
    slots.sort_by(|a, b| {
        b.0.partial_cmp(&a.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.1.cmp(&a.1))
            .then_with(|| a.2.cmp(&b.2))
    });
    // Two clockwise passes with a pending-arrival stack: a departure
    // matches the most recent pending arrival (wedge nesting); the
    // second pass closes pairs straddling the angular origin.
    if crate::trace_enabled_wedge() {
        let names: Vec<String> = slots
            .iter()
            .map(|(a, d, hh)| format!("{}{}@{a:.3}", if *d { "dep" } else { "arr" }, hh))
            .collect();
        eprintln!("      [wedge] bp {bp_key:?} slots {}", names.join(" "));
        for &(_, is_dep, hh) in &slots {
            if let Some(em) = sub.embeddings.get(&hh) {
                let n = em.uv.len();
                eprintln!(
                    "        [wedge-uv] {}{hh}: {} pts ({:.4},{:.4})..({:.4},{:.4})",
                    if is_dep { "dep" } else { "arr" },
                    n,
                    em.uv[0].0,
                    em.uv[0].1,
                    em.uv[n - 1].0,
                    em.uv[n - 1].1
                );
            }
        }
    }
    let mut pending: Vec<usize> = Vec::new();
    let mut matched: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut used_deps: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for _pass in 0..2 {
        for &(_, is_dep, hh) in &slots {
            if is_dep {
                if used_deps.contains(&hh) {
                    continue;
                }
                if let Some(a) = pending.pop() {
                    matched.insert(a, hh);
                    used_deps.insert(hh);
                }
            } else if !matched.contains_key(&hh) && !pending.contains(&hh) {
                pending.push(hh);
            }
        }
    }
    // Exact-tie degeneracy fallback: pair leftovers deterministically.
    let mut leftover_deps: Vec<usize> = deps
        .iter()
        .copied()
        .filter(|d| !used_deps.contains(d))
        .collect();
    leftover_deps.sort_unstable();
    pending.sort_unstable();
    for a in pending {
        if matched.contains_key(&a) {
            continue;
        }
        if let Some(d) = leftover_deps.pop() {
            matched.insert(a, d);
        }
    }
    // NEAR-TIE TWIN-MATCH REPAIR: when a record chain meets a face
    // boundary almost tangentially (a thread groove exiting through
    // a rim at an x.75 phase meets it 0.013 rad off the tangent),
    // the slot angles of the record pair and the boundary departure
    // tie within the pcurve-fit slop, the sweep order flips, and the
    // parenthesis matching pairs the record arrival with its OWN
    // TWIN — the walk then traverses the chain out-and-back as a
    // zero-area ring and the region beside it never gets a fragment
    // (measured: the thread cap's kept crescent vanished and the
    // whole cap classified away, stranding the mouth records —
    // PK-460). A twin pairing inside a near-tie cluster is the
    // degenerate signature; swapping with the tied pair restores
    // both rings. Far-from-tie twin pairings are LEGITIMATE U-turns
    // (the rotated torus × torus junction this matching was built
    // for) and are left alone.
    const NEAR_TIE: f64 = 0.02;
    let angle_of = |he: usize, is_dep: bool| -> Option<f64> {
        slots
            .iter()
            .find(|(_, d, hh)| *d == is_dep && *hh == he)
            .map(|(a, _, _)| *a)
    };
    let circ_close = |a: f64, b: f64| -> bool {
        let d = (a - b).abs();
        d.min(std::f64::consts::TAU - d) < NEAR_TIE
    };
    let mut twin_matched: Vec<(usize, usize)> = matched
        .iter()
        .filter(|(a, d)| graph.half_edges[**a].twin == **d)
        .map(|(a, d)| (*a, *d))
        .collect();
    // HashMap iteration order is not deterministic; the repair order
    // decides which tied pair swaps first (measured: repeat runs of
    // the determinism fixture differed in an edge window once the
    // residue records enlarged the tie population — PK-465).
    twin_matched.sort_unstable();
    for (a, d) in twin_matched {
        let (Some(a_ang), Some(d_ang)) = (angle_of(a, false), angle_of(d, true)) else {
            continue;
        };
        let swap = matched
            .iter()
            .filter(|(a2, _)| **a2 != a)
            .filter(|(a2, d2)| {
                let close = angle_of(**d2, true).is_some_and(|x| circ_close(x, d_ang))
                    || angle_of(**a2, false).is_some_and(|x| circ_close(x, a_ang));
                close && graph.half_edges[a].twin != **d2 && graph.half_edges[**a2].twin != d
            })
            .map(|(a2, d2)| (*a2, *d2))
            .min_by_key(|&(a2, _)| a2);
        if let Some((a2, d2)) = swap {
            matched.insert(a, d2);
            matched.insert(a2, d);
            if crate::trace_enabled_wedge() {
                eprintln!("      [wedge] near-tie twin repair at {bp_key:?}: {a}->{d2}, {a2}->{d}");
            }
        }
    }
    if crate::trace_enabled_wedge() {
        eprintln!("      [wedge] h {h} -> {:?}", matched.get(&h));
    }
    matched.get(&h).copied()
}
