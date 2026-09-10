//! Floating-component seam-config search and final branch clustering
//! (split from `chart_branches.rs` under the 500-line file budget).

use std::collections::{BTreeMap, BTreeSet, HashMap};

use super::chart_branches::clustering::{self, cluster_branches, is_balanced};
use super::chart_subdivision::{BranchPoint, FaceChartSubdivision, SubdivisionIssue};

/// Resolve floating components' border bumps (grouped bounded search)
/// and commit the final branch clustering into `sub`.
#[allow(clippy::too_many_arguments)]
pub(super) fn resolve_floating_and_cluster(
    sub: &mut FaceChartSubdivision,
    period: f64,
    ends: &dyn Fn(usize) -> (usize, usize),
    hes: &[usize],
    components: &[Vec<usize>],
    floating: &[usize],
    seam_nodes: &BTreeSet<usize>,
    true_seam_nodes: &BTreeSet<usize>,
    wrap_nodes: &BTreeSet<usize>,
    seam_hes: &BTreeSet<usize>,
) {
    sub.report.seam_search_used = !floating.is_empty();

    // Nodes touched by each floating component. A bump changes the
    // clustering ONLY at these nodes, so components interact — for
    // grouping AND for the balance verdict — exactly where their
    // node sets overlap (seam or not: rim junction nodes shared by
    // two components couple their choices the same way).
    let comp_seam_nodes = |ci: usize| -> BTreeSet<usize> {
        components[ci]
            .iter()
            .flat_map(|&h| {
                let (a, b) = ends(h);
                [a, b]
            })
            .collect()
    };
    // Union-find over floating components by shared seam nodes.
    let mut group_of: Vec<usize> = (0..floating.len()).collect();
    fn find(g: &mut Vec<usize>, i: usize) -> usize {
        if g[i] != i {
            let r = find(g, g[i]);
            g[i] = r;
        }
        g[i]
    }
    let seam_sets: Vec<BTreeSet<usize>> = floating.iter().map(|&ci| comp_seam_nodes(ci)).collect();
    for i in 0..floating.len() {
        for j in (i + 1)..floating.len() {
            if !seam_sets[i].is_disjoint(&seam_sets[j]) {
                let (ri, rj) = (find(&mut group_of, i), find(&mut group_of, j));
                if ri != rj {
                    group_of[rj] = ri;
                }
            }
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..floating.len() {
        let r = find(&mut group_of, i);
        groups.entry(r).or_default().push(i);
    }

    // Base u maps (body-side samples — see the comment below) without
    // any bumps.
    let base_maps = |bumped: &BTreeSet<usize>| -> (BTreeMap<usize, f64>, BTreeMap<usize, f64>) {
        let mut u_start: BTreeMap<usize, f64> = BTreeMap::new();
        let mut u_end: BTreeMap<usize, f64> = BTreeMap::new();
        for &h in hes {
            let emb = &sub.embeddings[&h];
            // Cluster by the BODY-SIDE sample, not the exact
            // endpoint: an endpoint landing exactly ON the chart
            // border (u = 0 with the run's body at negative u — the
            // exact-seeded junction class) clusters into the WRONG
            // border's branch by its own u, splitting a record chain
            // away from its seam-copy partners (measured under the
            // exact-torus NURBS arm: the cylinder wall's y<0
            // chain walked as a degenerate out-and-back region and
            // the union declined). The sample one step inside the run
            // carries the body's side; off-border endpoints cluster
            // identically either way.
            let n = emb.uv.len();
            let s_idx = if n > 2 { 1 } else { 0 };
            let e_idx = if n > 2 { n - 2 } else { n - 1 };
            u_start.insert(h, emb.uv[s_idx].0);
            u_end.insert(h, emb.uv[e_idx].0);
        }
        for &fi in bumped {
            for &h in &components[floating[fi]] {
                if let Some(u) = u_start.get_mut(&h) {
                    *u += period;
                }
                if let Some(u) = u_end.get_mut(&h) {
                    *u += period;
                }
            }
        }
        (u_start, u_end)
    };
    let cluster_with = |bumped: &BTreeSet<usize>| {
        let (u_start, u_end) = base_maps(bumped);
        cluster_branches(hes, ends, &u_start, &u_end, period, &sub.poles, wrap_nodes)
    };
    // Imbalance restricted to a node set.
    let balanced_at = |clustered: &(BTreeMap<usize, BranchPoint>, BTreeMap<usize, BranchPoint>),
                       nodes: &BTreeSet<usize>|
     -> bool {
        let mut deg: BTreeMap<BranchPoint, (i64, i64)> = BTreeMap::new();
        for bp in clustered.0.values() {
            deg.entry(*bp).or_default().1 += 1;
        }
        for bp in clustered.1.values() {
            deg.entry(*bp).or_default().0 += 1;
        }
        deg.iter()
            .filter(|(bp, _)| nodes.contains(&bp.node))
            .all(|(_, &(a, d))| a == d)
    };

    let mut exhausted = false;
    let bumped: BTreeSet<usize> = {
        let mut bumped: BTreeSet<usize> = BTreeSet::new();
        for members in groups.values() {
            if members.len() > 12 {
                exhausted = true;
                break;
            }
            let group_nodes: BTreeSet<usize> = members
                .iter()
                .flat_map(|&fi| seam_sets[fi].iter().copied())
                .collect();
            let n_configs = 1usize << members.len();
            let mut group_ok = false;
            for config in 0..n_configs {
                let mut trial = bumped.clone();
                for (bit, &fi) in members.iter().enumerate() {
                    if (config >> bit) & 1 == 1 {
                        trial.insert(fi);
                    }
                }
                let clustered = cluster_with(&trial);
                if balanced_at(&clustered, &group_nodes) {
                    bumped = trial;
                    group_ok = true;
                    break;
                }
            }
            // No balanced config for this group: the historical JOINT
            // search declared the whole chart exhausted in this state
            // (every joint config fails at these nodes), and the
            // decline is load-bearing — silently keeping the
            // unbumped choice walked a balanced-looking garbage
            // clustering on the coarse sphere−cube envelope pin and
            // SERVED an empty body (volume 0.0 vs 0.202).
            if !group_ok {
                exhausted = true;
                break;
            }
        }
        bumped
    };
    // The historical JOINT search accepted a config only when the
    // SEARCH-TIME clustering balanced at EVERY node — group-local
    // acceptance composed can differ (a baseline imbalance at a
    // non-float node fails every joint config, and the joint search
    // declined; the grouped composition sailed past it and the
    // coarse sphere−cube envelope pin SERVED an empty body). Verify
    // the composed config under the same bar.
    if !exhausted {
        let composed = cluster_with(&bumped);
        if !is_balanced(&composed.0, &composed.1) {
            exhausted = true;
        }
    }
    if exhausted {
        sub.report
            .issues
            .push(SubdivisionIssue::SeamSearchExhausted);
        return;
    }
    // Persist the winning bumps so directions and diagnostics agree
    // with the clustering.
    for &fi in &bumped {
        for &h in &components[floating[fi]] {
            if let Some(e) = sub.embeddings.get_mut(&h) {
                e.shift_u(period);
            }
        }
    }
    // Final clustering over the (now-shifted) embeddings.
    let chosen = {
        let mut u_start: BTreeMap<usize, f64> = BTreeMap::new();
        let mut u_end: BTreeMap<usize, f64> = BTreeMap::new();
        for &h in hes {
            let emb = &sub.embeddings[&h];
            let n = emb.uv.len();
            let s_idx = if n > 2 { 1 } else { 0 };
            let e_idx = if n > 2 { n - 2 } else { n - 1 };
            u_start.insert(h, emb.uv[s_idx].0);
            u_end.insert(h, emb.uv[e_idx].0);
        }
        Some(cluster_branches(
            hes, ends, &u_start, &u_end, period, &sub.poles, wrap_nodes,
        ))
    };
    let first_result: Option<(BTreeMap<usize, BranchPoint>, BTreeMap<usize, BranchPoint>)> = None;
    let (mut start_bp, mut end_bp) = match chosen.or(first_result) {
        Some(r) => r,
        None => return,
    };

    // When the component-level config search could not balance,
    // border-endpoint entries at TRUE seam nodes still carry two-copy
    // freedom the component bumps cannot express (a full-period rim
    // chain is anchored, but its endpoint AT the seam node may
    // register at either border). Node-local rebalance (PK-465).
    let pre_balanced = is_balanced(&start_bp, &end_bp);
    if !pre_balanced {
        clustering::rebalance_seam_border_entries(
            &mut start_bp,
            &mut end_bp,
            true_seam_nodes,
            seam_hes,
            &|h, is_start| {
                let emb = &sub.embeddings[&h];
                if is_start {
                    emb.uv[0].0
                } else {
                    emb.uv[emb.uv.len() - 1].0
                }
            },
            period,
        );
    }

    // GLOBAL balance is the acceptance bar the historical JOINT
    // search enforced: a per-group choice can balance every group's
    // own nodes yet leave the whole chart imbalanced, and silently
    // walking that clustering served an EMPTY body on the coarse
    // sphere−cube envelope pin (volume 0.0 vs 0.202 — a serving
    // lie). When even the node-local rebalance cannot close it,
    // decline the chart honestly.
    if !pre_balanced && !is_balanced(&start_bp, &end_bp) {
        sub.report
            .issues
            .push(SubdivisionIssue::SeamSearchExhausted);
        return;
    }

    // I2: off-seam nodes must be single-branch.
    let mut node_branches: HashMap<usize, BTreeSet<u32>> = HashMap::new();
    for (&h, bp) in start_bp.iter().chain(end_bp.iter()) {
        let _ = h;
        node_branches.entry(bp.node).or_default().insert(bp.branch);
    }
    for (node, branches) in &node_branches {
        // Wrap nodes cluster as one branch by construction; any
        // remaining split there is the seam machinery's business,
        // not an off-seam defect.
        if branches.len() > 1 && !seam_nodes.contains(node) && !wrap_nodes.contains(node) {
            sub.report
                .issues
                .push(SubdivisionIssue::BranchSplitOffSeam { node: *node });
        }
    }

    if popkern_core::trace::trace_enabled("branch-pairing") {
        let mut per_bp: BTreeMap<BranchPoint, Vec<String>> = BTreeMap::new();
        for (&h, &bp) in &start_bp {
            per_bp.entry(bp).or_default().push(format!("dep {h}"));
        }
        for (&h, &bp) in &end_bp {
            per_bp.entry(bp).or_default().push(format!("arr {h}"));
        }
        for (bp, members) in &per_bp {
            eprintln!("  [bp] n{} b{}: {:?}", bp.node, bp.branch, members);
        }
    }
    for (&h, &bp) in &start_bp {
        sub.departures.entry(bp).or_default().push(h);
    }
    sub.start_bp = start_bp;
    sub.end_bp = end_bp;
}
