//! Per-node branch clustering, degree-balance checks, duplicate
//! resolution, and the boundary shoelace σ calibration for chart
//! subdivisions.
//!
//! Split from `chart_branches.rs` to honour the 500-line budget.

use std::collections::BTreeMap;

use super::super::chart_subdivision::{BranchPoint, FaceChartSubdivision, SubdivisionIssue};
use crate::merge::graph::IntersectionGraph;

/// Cluster endpoint u values per node; branch index ascends with u.
#[allow(clippy::type_complexity)]
pub(crate) fn cluster_branches(
    hes: &[usize],
    ends: &dyn Fn(usize) -> (usize, usize),
    u_start: &BTreeMap<usize, f64>,
    u_end: &BTreeMap<usize, f64>,
    period: f64,
    poles: &BTreeMap<usize, f64>,
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
        // carry their meridians' azimuths — never split it by u.
        let is_pole = poles.contains_key(&node);
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
