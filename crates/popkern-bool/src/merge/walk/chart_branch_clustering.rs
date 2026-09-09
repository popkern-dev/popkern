//! Per-node branch clustering, degree-balance checks, duplicate
//! resolution, and the boundary shoelace σ calibration for chart
//! subdivisions.
//!
//! Split from `chart_branches.rs` to honour the 500-line budget.

use std::collections::BTreeMap;

use super::super::chart_subdivision::{BranchPoint, FaceChartSubdivision, SubdivisionIssue};
use crate::merge::graph::IntersectionGraph;

/// Cluster endpoint u values per node; branch index ascends with u.
///
/// `wrap_nodes`: period-lattice junction nodes with NO seam edge —
/// the overwrap retry's barrier set minus the true seam nodes. Their
/// u = 0 and u = period endpoint copies are the SAME azimuth (a
/// boundary helix run split per cell arrives at one border copy and
/// departs at the other), so they cluster as ONE branch, like a
/// pole: keeping them split left every per-turn junction of the
/// nut's bore-crest band 1-in/0-out + 0-in/1-out and both caps of
/// the walk declined (PK-464). True seam nodes keep their split —
/// the seam-edge pair owns the border crossing there.
#[allow(clippy::type_complexity)]
pub(crate) fn cluster_branches(
    hes: &[usize],
    ends: &dyn Fn(usize) -> (usize, usize),
    u_start: &BTreeMap<usize, f64>,
    u_end: &BTreeMap<usize, f64>,
    period: f64,
    poles: &BTreeMap<usize, f64>,
    wrap_nodes: &std::collections::BTreeSet<usize>,
) -> (BTreeMap<usize, BranchPoint>, BTreeMap<usize, BranchPoint>) {
    // node -> [(u, he, is_start)]
    let mut per_node: BTreeMap<usize, Vec<(f64, usize, bool)>> = BTreeMap::new();
    for &h in hes {
        let (s, e) = ends(h);
        per_node.entry(s).or_default().push((u_start[&h], h, true));
        per_node.entry(e).or_default().push((u_end[&h], h, false));
    }
    let mut start_bp = BTreeMap::new();
    let mut end_bp = BTreeMap::new();
    let trace = popkern_core::trace::trace_enabled("chart-branches");
    for (node, mut entries) in per_node {
        if trace {
            for &(u, h, is_start) in &entries {
                eprintln!(
                    "  [branch] node {node}: he {h} u={u:.4} {}",
                    if is_start { "START" } else { "END" }
                );
            }
        }
        // A POLE node is one geometric point whose incident endpoints
        // carry their meridians' azimuths — never split it by u. A
        // WRAP node's copies are one azimuth mod the period — same.
        let is_pole = poles.contains_key(&node) || wrap_nodes.contains(&node);
        entries.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut branch = 0u32;
        let mut prev_u: Option<f64> = None;
        for (u, h, is_start) in entries {
            if !is_pole
                && let Some(p) = prev_u
                && u - p > period * 0.5
            {
                branch += 1;
            }
            prev_u = Some(u);
            let bp = BranchPoint { node, branch };
            if is_start {
                start_bp.insert(h, bp);
            } else {
                end_bp.insert(h, bp);
            }
        }
    }
    (start_bp, end_bp)
}

/// Per-node rebalance of border-endpoint entries at TRUE seam nodes.
///
/// A non-seam run terminating exactly ON the chart border (endpoint u
/// on the period lattice) at a seam node can legitimately register at
/// either border copy — u = 0 and u = period are the same azimuth.
/// The floating-component search cannot express this freedom when the
/// run's component spans a full period (a split rim chain is anchored,
/// yet its endpoint AT the seam node still has two-copy freedom).
/// Measured on the angled tap (PK-465): the tilted wall's top-rim
/// pieces both clustered by body-side sample (arrival low, departure
/// high) while the seam pair also splits arrival-low/departure-high,
/// leaving node 9 with branches 2-in/0-out and 0-in/2-out; the only
/// balanced assignment is the cross-pairing.
///
/// The repair is node-local — an entry belongs to exactly one node, so
/// rebalancing each imbalanced seam node independently composes to
/// global balance. Per node: exactly two branches, candidates are
/// non-seam entries whose EXACT endpoint u lies on the lattice
/// (within `period * 5e-4`, the anchoring lattice eps), bounded
/// enumeration (≤ 8 candidates, fewest flips first, deterministic
/// tie order), first assignment balancing the node wins. A wrong but
/// balanced flip is the same ambiguity class the seam config search
/// already navigates by first-balanced-wins; downstream walk
/// validation still gates the result.
pub(crate) fn rebalance_seam_border_entries(
    start_bp: &mut BTreeMap<usize, BranchPoint>,
    end_bp: &mut BTreeMap<usize, BranchPoint>,
    true_seam_nodes: &std::collections::BTreeSet<usize>,
    seam_hes: &std::collections::BTreeSet<usize>,
    endpoint_u: &dyn Fn(usize, bool) -> f64,
    period: f64,
) {
    let lattice_eps = period * 5.0e-4;
    let trace = popkern_core::trace::trace_enabled("chart-branches");
    // node -> [(he, is_start, branch)]
    let mut per_node: BTreeMap<usize, Vec<(usize, bool, u32)>> = BTreeMap::new();
    for (&h, bp) in start_bp.iter() {
        per_node
            .entry(bp.node)
            .or_default()
            .push((h, true, bp.branch));
    }
    for (&h, bp) in end_bp.iter() {
        per_node
            .entry(bp.node)
            .or_default()
            .push((h, false, bp.branch));
    }
    let node_balanced = |ents: &[(usize, bool, u32)]| {
        let mut deg: BTreeMap<u32, (i64, i64)> = BTreeMap::new();
        for &(_, is_start, b) in ents {
            let e = deg.entry(b).or_default();
            if is_start {
                e.1 += 1;
            } else {
                e.0 += 1;
            }
        }
        deg.values().all(|&(arr, dep)| arr == dep)
    };
    for (node, entries) in per_node {
        if !true_seam_nodes.contains(&node) || node_balanced(&entries) {
            continue;
        }
        let branches: std::collections::BTreeSet<u32> = entries.iter().map(|e| e.2).collect();
        if branches.len() != 2 {
            continue;
        }
        let bs: Vec<u32> = branches.into_iter().collect();
        let cand: Vec<usize> = entries
            .iter()
            .enumerate()
            .filter(|&(_, &(h, is_start, _))| {
                if seam_hes.contains(&h) {
                    return false;
                }
                let u = endpoint_u(h, is_start);
                let k = (u / period).round();
                (u - k * period).abs() < lattice_eps
            })
            .map(|(i, _)| i)
            .collect();
        if cand.is_empty() || cand.len() > 8 {
            continue;
        }
        let mut masks: Vec<u32> = (1..(1u32 << cand.len())).collect();
        masks.sort_by_key(|m| (m.count_ones(), *m));
        for m in masks {
            let mut trial = entries.clone();
            for (bit, &ei) in cand.iter().enumerate() {
                if (m >> bit) & 1 == 1 {
                    let e = &mut trial[ei];
                    e.2 = if e.2 == bs[0] { bs[1] } else { bs[0] };
                }
            }
            if !node_balanced(&trial) {
                continue;
            }
            for &(h, is_start, b) in &trial {
                let bp = BranchPoint { node, branch: b };
                if is_start {
                    start_bp.insert(h, bp);
                } else {
                    end_bp.insert(h, bp);
                }
            }
            if trace {
                eprintln!(
                    "  [branch] seam-border rebalance at node {node}: mask {m:#b} over {} candidate entr{}",
                    cand.len(),
                    if cand.len() == 1 { "y" } else { "ies" }
                );
            }
            break;
        }
    }
}

/// I4 predicate: every branch point has in-degree == out-degree.
pub(crate) fn is_balanced(
    start_bp: &BTreeMap<usize, BranchPoint>,
    end_bp: &BTreeMap<usize, BranchPoint>,
) -> bool {
    let mut deg: BTreeMap<BranchPoint, (usize, usize)> = BTreeMap::new();
    for bp in start_bp.values() {
        deg.entry(*bp).or_default().1 += 1;
    }
    for bp in end_bp.values() {
        deg.entry(*bp).or_default().0 += 1;
    }
    deg.values().all(|&(a, d)| a == d)
}

/// Final I4/I3/I5 pass over the built maps, recording issues.
pub(crate) fn balance_and_duplicates(sub: &mut FaceChartSubdivision, _u_period: Option<f64>) {
    let mut deg: BTreeMap<BranchPoint, (usize, usize)> = BTreeMap::new();
    for bp in sub.start_bp.values() {
        deg.entry(*bp).or_default().1 += 1;
    }
    for bp in sub.end_bp.values() {
        deg.entry(*bp).or_default().0 += 1;
    }
    for (bp, &(arrivals, departures)) in &deg {
        if arrivals != departures {
            sub.report.issues.push(SubdivisionIssue::DegreeImbalance {
                node: bp.node,
                branch: bp.branch,
                arrivals,
                departures,
            });
            // I3 refinement: an arrival-only point sitting on another
            // run's interior is an uncompleted T-junction.
            if departures == 0
                && arrivals > 0
                && let Some((&h_arr, _)) = sub.end_bp.iter().find(|(_, b)| **b == *bp)
            {
                let emb = &sub.embeddings[&h_arr];
                let p = emb.uv[emb.uv.len() - 1];
                if node_on_some_run(sub, h_arr, p) {
                    sub.report
                        .issues
                        .push(SubdivisionIssue::UnsplitTJunction { node: bp.node });
                }
            }
        }
    }
    // I5: duplicate departures (same direction, same far point).
    let departures = sub.departures.clone();
    for (_, deps) in departures {
        for i in 0..deps.len() {
            for j in (i + 1)..deps.len() {
                let (a, b) = (deps[i], deps[j]);
                if sub.end_bp.get(&a) != sub.end_bp.get(&b) {
                    continue;
                }
                if let (Some(da), Some(db)) = (
                    sub.embeddings[&a].start_dir(),
                    sub.embeddings[&b].start_dir(),
                ) {
                    let cross = da.0 * db.1 - da.1 * db.0;
                    let dot = da.0 * db.0 + da.1 * db.1;
                    if cross.abs() < 1e-9 && dot > 0.0 {
                        sub.report
                            .issues
                            .push(SubdivisionIssue::DuplicateChartSegment(a, b));
                    }
                }
            }
        }
    }
}

/// Chart handedness from the stored boundary: the closed union of the
/// face's non-Shared (original/split boundary) runs, traversed in
/// registered direction, has shoelace sign +1 exactly when chart-CCW
/// equals CCW-about-outward (stored loops wind CCW about outward —
/// winding census E1). Returns `None` when the face has no boundary
/// runs or the area is degenerate; caller keeps the orientation-flag
/// fallback. The sum is order-independent because the boundary is a
/// closed union of loops and each directed segment contributes its
/// own wedge term.
pub(crate) fn boundary_shoelace_sigma(
    graph: &IntersectionGraph,
    sub: &FaceChartSubdivision,
) -> Option<f64> {
    use crate::merge::graph::EdgeSource;
    let mut area2 = 0.0;
    let mut any = false;
    let mut scale = 0.0_f64;
    for (&h, emb) in &sub.embeddings {
        if matches!(graph.half_edges[h].source, EdgeSource::Shared { .. }) {
            continue;
        }
        for w in emb.uv.windows(2) {
            area2 += w[0].0 * w[1].1 - w[1].0 * w[0].1;
            scale = scale.max(w[0].0.abs() + w[0].1.abs());
            any = true;
        }
    }
    if !any || area2.abs() < 1e-12 * (1.0 + scale * scale) {
        None
    } else {
        Some(area2.signum())
    }
}

/// Non-periodic charts: every node is a single branch point.
pub(crate) fn cluster_single_branch(
    sub: &mut FaceChartSubdivision,
    ends: &dyn Fn(usize) -> (usize, usize),
) {
    let hes: Vec<usize> = sub.embeddings.keys().copied().collect();
    for h in hes {
        let (s, e) = ends(h);
        let s_bp = BranchPoint { node: s, branch: 0 };
        let e_bp = BranchPoint { node: e, branch: 0 };
        sub.start_bp.insert(h, s_bp);
        sub.end_bp.insert(h, e_bp);
        sub.departures.entry(s_bp).or_default().push(h);
    }
}

/// True when chart point `p` lies on the interior of some OTHER
/// registered half-edge's run (coarse, diagnostics-only).
pub(crate) fn node_on_some_run(sub: &FaceChartSubdivision, exclude: usize, p: (f64, f64)) -> bool {
    let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
    for emb in sub.embeddings.values() {
        for &(u, v) in &emb.uv {
            lo = (lo.0.min(u), lo.1.min(v));
            hi = (hi.0.max(u), hi.1.max(v));
        }
    }
    let eps = 1e-6 * (1.0 + (hi.0 - lo.0).abs() + (hi.1 - lo.1).abs());
    let ex_twin = exclude; // run of the arriving he itself is fine to skip only
    for (&h, emb) in &sub.embeddings {
        if h == ex_twin {
            continue;
        }
        for w in emb.uv.windows(2) {
            let (a, b) = (w[0], w[1]);
            let (dx, dy) = (b.0 - a.0, b.1 - a.1);
            let l2 = (dx * dx + dy * dy).max(1e-300);
            let t = (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / l2).clamp(0.0, 1.0);
            let d = ((p.0 - a.0 - dx * t).powi(2) + (p.1 - a.1 - dy * t).powi(2)).sqrt();
            if d < eps && t > 1e-6 && t < 1.0 - 1e-6 {
                return true;
            }
        }
    }
    false
}
