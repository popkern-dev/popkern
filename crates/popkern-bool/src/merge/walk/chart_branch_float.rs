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
    // Cut-domain v keys (PK-459): the same near-end samples as u.
    // Taken once: the border bumps below move runs in u only.
    let v_keys: BTreeMap<usize, (f64, f64)> = sub
        .embeddings
        .iter()
        .map(|(&h, emb)| {
            let (uv, n) = (&emb.uv, emb.uv.len());
            let (s, e) = if n > 2 { (1, n - 2) } else { (0, n - 1) };
            (h, (uv[s].1, uv[e].1))
        })
        .collect();
    let v_key = |h: usize, start: bool| -> f64 {
        let (s, e) = v_keys.get(&h).copied().unwrap_or((0.0, 0.0));
        if start { s } else { e }
    };
    let v_split = sub
        .cut_domain
        .as_ref()
        .map(|d| d.v_period)
        .map(|pv| (&v_key as &dyn Fn(usize, bool) -> f64, pv));
    let cluster_with = |bumped: &BTreeSet<usize>| {
        let (u_start, u_end) = base_maps(bumped);
        cluster_branches(
            hes, ends, &u_start, &u_end, period, &sub.poles, wrap_nodes, v_split,
        )
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

    let float_trace = popkern_core::trace::trace_enabled("branch-pairing");
    let imbalance_of = |clustered: &(
        BTreeMap<usize, BranchPoint>,
        BTreeMap<usize, BranchPoint>,
    )|
     -> Vec<(BranchPoint, i64, i64)> {
        let mut deg: BTreeMap<BranchPoint, (i64, i64)> = BTreeMap::new();
        for bp in clustered.0.values() {
            deg.entry(*bp).or_default().1 += 1;
        }
        for bp in clustered.1.values() {
            deg.entry(*bp).or_default().0 += 1;
        }
        deg.into_iter()
            .filter(|(_, (a, d))| a != d)
            .map(|(bp, (a, d))| (bp, a, d))
            .collect()
    };
    let mut exhausted = false;
    let bumped: BTreeSet<usize> = {
        let mut bumped: BTreeSet<usize> = BTreeSet::new();
        for members in groups.values() {
            if members.len() > 12 {
                if float_trace {
                    eprintln!("  [float] exhausted: group of {} > 12", members.len());
                }
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
                if float_trace {
                    let imb = imbalance_of(&cluster_with(&bumped));
                    eprintln!(
                        "  [float] exhausted: no balanced config for group {members:?}; unbumped imbalance {imb:?}"
                    );
                }
                exhausted = true;
                break;
            }
        }
        bumped
    };
    // The composed config is NOT re-verified here: a residual
    // imbalance after per-group success (or with no floating
    // components at all) is the seam-border ENTRY class — border
    // endpoints at true seam nodes carry two-copy freedom no
    // component bump can express — and the node-local rebalance
    // below owns it, with the same balance bar declining whatever it
    // cannot close (measured, PK-485: a σ = −1 chamfered-puck wall
    // carried zero floating components and every seam-node entry on
    // one border — a pre-rebalance composed check declined the chart
    // to the 3D-Newell walker, whose mega-loop double-claimed the
    // slot chain's twin sides). The group-search failure above keeps
    // its early decline: silently keeping an unbumped choice walked
    // a balanced-looking garbage clustering on the coarse
    // sphere−cube envelope pin and SERVED an empty body.
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
            hes, ends, &u_start, &u_end, period, &sub.poles, wrap_nodes, v_split,
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
        if float_trace {
            eprintln!(
                "  [float] exhausted: imbalance survives the seam-border rebalance {:?}",
                imbalance_of(&(start_bp.clone(), end_bp.clone()))
            );
        }
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
    // A cut domain's v-seam nodes split by v copy (PK-459).
    let v_seam_nodes: BTreeSet<usize> = sub
        .cut_domain
        .iter()
        .flat_map(|d| d.v_seam.iter())
        .flat_map(|&h| {
            let (a, b) = ends(h);
            [a, b]
        })
        .collect();
    for (node, branches) in &node_branches {
        // Wrap nodes cluster as one branch by construction; any
        // remaining split there is the seam machinery's business,
        // not an off-seam defect.
        if branches.len() > 1
            && !seam_nodes.contains(node)
            && !wrap_nodes.contains(node)
            && !v_seam_nodes.contains(node)
        {
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
