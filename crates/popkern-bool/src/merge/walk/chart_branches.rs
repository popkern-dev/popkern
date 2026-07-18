//! Period-branch resolution for chart subdivisions (PK-104).
//!
//! Split from `chart_subdivision.rs` to honour the 500-line budget:
//! consensus anchoring of embedding components, deterministic seam
//! border assignment (material-side rule), the bounded floating-
//! component search, per-node branch clustering, and the boundary
//! shoelace σ calibration.

use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};

use super::chart_subdivision::{
    BranchPoint, FaceChartSubdivision, SubdivisionIssue, node_on_some_run,
};
use crate::merge::graph::IntersectionGraph;

/// Periodic charts: consensus anchoring, seam border assignment with a
/// bounded balance search, then per-node u-clustering into branches.
pub(super) fn anchor_periodic(
    graph: &IntersectionGraph,
    sub: &mut FaceChartSubdivision,
    period: f64,
    ends: &dyn Fn(usize) -> (usize, usize),
) {
    let hes: Vec<usize> = sub.embeddings.keys().copied().collect();
    let seam_hes: BTreeSet<usize> = hes
        .iter()
        .copied()
        .filter(|&h| sub.embeddings[&h].seam_like)
        .collect();
    let seam_nodes: BTreeSet<usize> = seam_hes
        .iter()
        .flat_map(|&h| {
            let (s, e) = ends(h);
            [s, e]
        })
        .collect();

    // Consensus anchoring: BFS over interior (non-seam) hes through
    // shared non-seam nodes, translating each newcomer by whole
    // periods to match the frame already established.
    let interior: Vec<usize> = hes
        .iter()
        .copied()
        .filter(|h| !seam_hes.contains(h))
        .collect();
    let mut node_to_hes: HashMap<usize, Vec<usize>> = HashMap::new();
    for &h in &interior {
        let (s, e) = ends(h);
        for n in [s, e] {
            // Pole nodes are excluded like seam nodes: their endpoint
            // u is the incident meridian's azimuth, not a shared
            // chart coordinate — translating a neighbour to match it
            // would smear the frame (the boundary-first unwrap
            // learned the same lesson with degenerate anchors).
            if !seam_nodes.contains(&n) && !sub.poles.contains_key(&n) {
                node_to_hes.entry(n).or_default().push(h);
            }
        }
    }
    let u_at = |sub: &FaceChartSubdivision, h: usize, n: usize| -> f64 {
        let (s, _) = ends(h);
        let emb = &sub.embeddings[&h];
        if s == n {
            emb.uv[0].0
        } else {
            emb.uv[emb.uv.len() - 1].0
        }
    };
    let mut comp_of: HashMap<usize, usize> = HashMap::new();
    let mut components: Vec<Vec<usize>> = Vec::new();
    for &seed in &interior {
        if comp_of.contains_key(&seed) {
            continue;
        }
        let comp_id = components.len();
        let mut members = vec![seed];
        comp_of.insert(seed, comp_id);
        let mut queue = VecDeque::from([seed]);
        while let Some(h) = queue.pop_front() {
            let (s, e) = ends(h);
            for n in [s, e] {
                if seam_nodes.contains(&n) {
                    continue;
                }
                for &g in node_to_hes.get(&n).into_iter().flatten() {
                    if comp_of.contains_key(&g) {
                        continue;
                    }
                    let delta = u_at(sub, h, n) - u_at(sub, g, n);
                    let k = (delta / period).round();
                    if k != 0.0
                        && let Some(emb) = sub.embeddings.get_mut(&g)
                    {
                        emb.shift_u(k * period);
                    }
                    comp_of.insert(g, comp_id);
                    members.push(g);
                    queue.push_back(g);
                }
            }
        }
        components.push(members);
    }

    // Reference frame: the largest component (ties: lowest first he).
    let ref_comp = components
        .iter()
        .enumerate()
        .max_by_key(|(i, m)| (m.len(), usize::MAX - i))
        .map(|(i, _)| i);
    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    if let Some(rc) = ref_comp {
        for &h in &components[rc] {
            for &(u, _) in &sub.embeddings[&h].uv {
                lo = lo.min(u);
                hi = hi.max(u);
            }
        }
    }
    if hi > lo && hi - lo > period * (1.0 + 1e-6) {
        sub.report.issues.push(SubdivisionIssue::ChartOverwrap {
            span: hi - lo,
            period,
        });
        return;
    }

    // Seam border assignment is DETERMINISTIC, not searched: the two
    // chart copies of the seam differ only by which travel direction
    // keeps the chart interior on the walk's material side. For
    // material-left in oriented UV (σ = +1), the interior (u > b_lo)
    // is on the left of v-DECREASING travel at the low border and of
    // v-INCREASING travel at the high border; σ = −1 mirrors this.
    // (Measured: a balance-only search accepted the lockstep-swapped
    // assignment on the breached-bore wall — also balanced, but it
    // re-created the legacy 16-half-edge mega-loop orbit.)
    let mut seam_pairs: Vec<(usize, usize)> = Vec::new();
    for &h in &seam_hes {
        let tw = graph.half_edges[h].twin;
        if h < tw && seam_hes.contains(&tw) {
            seam_pairs.push((h, tw));
        }
    }
    // Window the low border must land in: [hi - period, lo]. Without
    // an interior reference, anchor on the first seam he's own u.
    let (window_lo_min, window_lo_max) = if hi > lo {
        (hi - period - 1e-9, lo + 1e-9)
    } else {
        let u0 = seam_pairs
            .first()
            .map(|&(h, _)| sub.embeddings[&h].uv[0].0)
            .unwrap_or(0.0);
        (u0 - 1e-9, u0 + 1e-9)
    };
    // Snap a u value to the LOW border: choose k with u + k·period in
    // the low window; None when no k fits.
    let low_shift = |u: f64| -> Option<f64> {
        for k in -3i32..=3 {
            let shifted = u + k as f64 * period;
            if shifted >= window_lo_min && shifted <= window_lo_max {
                return Some(k as f64 * period);
            }
        }
        None
    };
    let mut sigma = sub.sigma;
    let seam_trace = std::env::var_os("PK_SEAM_TRACE").is_some();
    let mut assigned_high: BTreeMap<usize, bool> = BTreeMap::new();
    for &(h, tw) in &seam_pairs {
        for he in [h, tw] {
            let emb = &sub.embeddings[&he];
            let dv = emb.uv[emb.uv.len() - 1].1 - emb.uv[0].1;
            let u = emb.uv[0].0;
            let Some(to_low) = low_shift(u) else {
                sub.report
                    .issues
                    .push(SubdivisionIssue::SeamBorderMisfit { he });
                return;
            };
            let wants_high = (dv > 0.0) == (sigma > 0.0);
            if seam_trace {
                eprintln!(
                    "  [seam] he {he} dv {dv:+.3} u0 {u:.3} sigma {sigma:+.0} -> {}",
                    if wants_high { "HIGH" } else { "low" }
                );
            }
            assigned_high.insert(he, wants_high);
            let shift = if wants_high { to_low + period } else { to_low };
            if shift != 0.0
                && let Some(e) = sub.embeddings.get_mut(&he)
            {
                e.shift_u(shift);
            }
        }
    }
    // σ self-calibration (see the non-periodic branch): with seams
    // provisionally assigned, the full stored boundary is closed in
    // the chart and its shoelace sign is the true handedness. A
    // contradiction flips σ AND the seam borders (the dv-rule
    // mirrors with σ).
    if let Some(s) = boundary_shoelace_sigma(graph, sub)
        && s != sigma
    {
        sigma = s;
        sub.sigma = s;
        for (&he, &was_high) in &assigned_high {
            if let Some(e) = sub.embeddings.get_mut(&he) {
                e.shift_u(if was_high { -period } else { period });
            }
        }
    }
    let _ = sigma;
    // Non-reference components were, by construction, unreachable from
    // the reference through interior nodes — they interact with the
    // walk only at seam nodes (rings crossing the seam are the
    // canonical case) or not at all. Base-shift each so its min-u
    // lands on the low border window. A FULL-PERIOD component (rim or
    // trim ring) then has no remaining freedom: its endpoints span
    // both borders exactly. Only sub-period components hanging at a
    // seam node keep a genuine two-way branch choice for the config
    // search below (a drill wall carries four ring components, and
    // searching all of them blew the config budget for no reason —
    // none had any freedom).
    let mut floating: Vec<usize> = Vec::new();
    for ci in (0..components.len()).filter(|&i| Some(i) != ref_comp) {
        let (mut cmin, mut cmax) = (f64::MAX, f64::MIN);
        for h in &components[ci] {
            for &(u, _) in &sub.embeddings[h].uv {
                cmin = cmin.min(u);
                cmax = cmax.max(u);
            }
        }
        if let Some(shift) = low_shift(cmin)
            && shift != 0.0
        {
            for &h in &components[ci] {
                if let Some(emb) = sub.embeddings.get_mut(&h) {
                    emb.shift_u(shift);
                }
            }
        }
        // Branch freedom exists only where the component can MEET the
        // rest of the walk: at a seam node. An isolated interior
        // island (a trim ring, or a record loop imprinted mid-chart)
        // walks identically at either bump — searching its branches
        // doubles the config space for nothing (measured: the cube
        // imprint on a sphere carries EIGHT corner-triangle islands
        // and the 2^8 search tripped the budget; none had freedom).
        let touches_seam = components[ci].iter().any(|&h| {
            let (a, b) = ends(h);
            seam_nodes.contains(&a) || seam_nodes.contains(&b)
        });
        if touches_seam && cmax - cmin < period * (1.0 - 1e-9) {
            floating.push(ci);
        }
    }

    // The only remaining degree of freedom: floating components whose
    // span is below a full period can sit at either border of a seam
    // node. Bounded search over their ±period bumps, first balanced
    // config wins (a full-period rim has no real freedom — the wrong
    // bump is never balanced).
    if floating.len() > 4 {
        sub.report
            .issues
            .push(SubdivisionIssue::SeamSearchExhausted);
        return;
    }
    let n_configs = 1usize << floating.len();
    sub.report.seam_search_used = n_configs > 1;

    let mut chosen: Option<(BTreeMap<usize, BranchPoint>, BTreeMap<usize, BranchPoint>)> = None;
    let mut first_result: Option<(BTreeMap<usize, BranchPoint>, BTreeMap<usize, BranchPoint>)> =
        None;
    for config in 0..n_configs {
        let mut u_start: BTreeMap<usize, f64> = BTreeMap::new();
        let mut u_end: BTreeMap<usize, f64> = BTreeMap::new();
        for &h in &hes {
            let emb = &sub.embeddings[&h];
            u_start.insert(h, emb.uv[0].0);
            u_end.insert(h, emb.uv[emb.uv.len() - 1].0);
        }
        for (fi, &ci) in floating.iter().enumerate() {
            if (config >> fi) & 1 == 1 {
                for &h in &components[ci] {
                    if let Some(u) = u_start.get_mut(&h) {
                        *u += period;
                    }
                    if let Some(u) = u_end.get_mut(&h) {
                        *u += period;
                    }
                }
            }
        }
        let clustered = cluster_branches(&hes, ends, &u_start, &u_end, period, &sub.poles);
        if first_result.is_none() {
            first_result = Some(clustered.clone());
        }
        if is_balanced(&clustered.0, &clustered.1) {
            // Persist the winning bumps so directions and diagnostics
            // agree with the clustering.
            for (fi, &ci) in floating.iter().enumerate() {
                if (config >> fi) & 1 == 1 {
                    for &h in &components[ci] {
                        if let Some(e) = sub.embeddings.get_mut(&h) {
                            e.shift_u(period);
                        }
                    }
                }
            }
            chosen = Some(clustered);
            break;
        }
    }
    let (start_bp, end_bp) = match chosen.or(first_result) {
        Some(r) => r,
        None => return,
    };

    // I2: off-seam nodes must be single-branch.
    let mut node_branches: HashMap<usize, BTreeSet<u32>> = HashMap::new();
    for (&h, bp) in start_bp.iter().chain(end_bp.iter()) {
        let _ = h;
        node_branches.entry(bp.node).or_default().insert(bp.branch);
    }
    for (node, branches) in &node_branches {
        if branches.len() > 1 && !seam_nodes.contains(node) {
            sub.report
                .issues
                .push(SubdivisionIssue::BranchSplitOffSeam { node: *node });
        }
    }

    for (&h, &bp) in &start_bp {
        sub.departures.entry(bp).or_default().push(h);
    }
    sub.start_bp = start_bp;
    sub.end_bp = end_bp;
}

/// Cluster endpoint u values per node; branch index ascends with u.
#[allow(clippy::type_complexity)]
pub(super) fn cluster_branches(
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
    for (node, mut entries) in per_node {
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
pub(super) fn is_balanced(
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
pub(super) fn balance_and_duplicates(sub: &mut FaceChartSubdivision, _u_period: Option<f64>) {
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
pub(super) fn boundary_shoelace_sigma(
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
