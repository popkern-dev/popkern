//! Boundary-chain discovery for KP20 support-cell shared trims.
//!
//! A *chain* is a maximal set of boundary rows glued together by
//! shared trims along one geometric curve. Each trim glues an
//! interval of one row to an interval of another (traversed in
//! opposite directions), so composing the affine interval maps
//! embeds every row of the chain into one common chain
//! parameter. A chain is then split at *clean nodes* — chain
//! parameters that are row endpoints on every row containing
//! them — so each sub-chain can carry a single polynomial curve
//! that every member row restricts exactly. Closed chains (rings
//! around the polar cap) are detected by their embedding period
//! and rebased so the cut lands on a clean node.

use std::collections::HashMap;
use std::collections::VecDeque;
use std::hash::Hash;

/// Relative parameter tolerance for embedding consistency, node
/// deduplication, and clean-node interiority. Chain parameters
/// are in row-length units (each row spans 1.0), so this is a
/// dimensionless parametric guard, not a model tolerance.
const CHAIN_PARAM_TOL: f64 = 1.0e-6;

/// One shared trim, expressed as the glue between two row
/// intervals. `node_start` identifies the shared-edge endpoint
/// at `interval_a.0` (and `interval_b.1`); `node_end` the one at
/// `interval_a.1` (and `interval_b.0`).
#[derive(Debug, Clone, Copy)]
pub(super) struct ChainTrim<K, N> {
    pub(super) row_a: K,
    pub(super) interval_a: (f64, f64),
    pub(super) row_b: K,
    pub(super) interval_b: (f64, f64),
    pub(super) node_start: N,
    pub(super) node_end: N,
}

/// Affine embedding of one row into chain parameters:
/// `chain = scale * row_param + offset`.
#[derive(Debug, Clone, Copy)]
pub(super) struct RowChainEmbedding {
    pub(super) scale: f64,
    pub(super) offset: f64,
}

impl RowChainEmbedding {
    pub(super) fn map(self, row_param: f64) -> f64 {
        self.scale * row_param + self.offset
    }

    pub(super) fn chain_interval(self) -> (f64, f64) {
        let a = self.map(0.0);
        let b = self.map(1.0);
        (a.min(b), a.max(b))
    }
}

/// A fittable run of rows between two clean nodes. `nodes` is
/// sorted by chain parameter and includes both ends.
#[derive(Debug, Clone)]
pub(super) struct SubChain<K, N> {
    pub(super) rows: Vec<(K, RowChainEmbedding)>,
    pub(super) nodes: Vec<(N, f64)>,
}

#[derive(Debug, Clone)]
pub(super) struct ChainDiscovery<K, N> {
    pub(super) sub_chains: Vec<SubChain<K, N>>,
    pub(super) declined_chain_count: usize,
    pub(super) decline_reason: &'static str,
}

/// Discover all chains and split them at clean nodes. Chains
/// that cannot be embedded or split are counted as declined with
/// the first reason retained; callers fall back to the existing
/// per-trim reconciliation for their rows.
pub(super) fn discover_sub_chains<K, N>(trims: &[ChainTrim<K, N>]) -> ChainDiscovery<K, N>
where
    K: Copy + Eq + Hash,
    N: Copy + Eq + Hash,
{
    let mut adjacency: HashMap<K, Vec<usize>> = HashMap::new();
    for (index, trim) in trims.iter().enumerate() {
        adjacency.entry(trim.row_a).or_default().push(index);
        adjacency.entry(trim.row_b).or_default().push(index);
    }
    let mut visited: HashMap<K, RowChainEmbedding> = HashMap::new();
    let mut discovery = ChainDiscovery {
        sub_chains: Vec::new(),
        declined_chain_count: 0,
        decline_reason: "none",
    };
    for trim in trims {
        if visited.contains_key(&trim.row_a) {
            continue;
        }
        match embed_chain(trim.row_a, trims, &adjacency, &mut visited) {
            Ok(chain) => match split_chain_at_clean_nodes(chain) {
                Ok(mut sub_chains) => discovery.sub_chains.append(&mut sub_chains),
                Err(reason) => {
                    discovery.declined_chain_count += 1;
                    if discovery.decline_reason == "none" {
                        discovery.decline_reason = reason;
                    }
                }
            },
            Err(reason) => {
                discovery.declined_chain_count += 1;
                if discovery.decline_reason == "none" {
                    discovery.decline_reason = reason;
                }
            }
        }
    }
    discovery
}

struct EmbeddedChain<K, N> {
    rows: Vec<(K, RowChainEmbedding)>,
    nodes: Vec<(N, f64)>,
    period: Option<f64>,
}

fn embed_chain<K, N>(
    seed: K,
    trims: &[ChainTrim<K, N>],
    adjacency: &HashMap<K, Vec<usize>>,
    visited: &mut HashMap<K, RowChainEmbedding>,
) -> Result<EmbeddedChain<K, N>, &'static str>
where
    K: Copy + Eq + Hash,
    N: Copy + Eq + Hash,
{
    let mut rows = Vec::new();
    let mut node_params: HashMap<N, f64> = HashMap::new();
    let mut period: Option<f64> = None;
    let mut queue = VecDeque::new();
    let seed_embedding = RowChainEmbedding {
        scale: 1.0,
        offset: 0.0,
    };
    visited.insert(seed, seed_embedding);
    rows.push((seed, seed_embedding));
    queue.push_back(seed);
    while let Some(row) = queue.pop_front() {
        let embedding = visited[&row];
        let Some(trim_indices) = adjacency.get(&row) else {
            continue;
        };
        for &trim_index in trim_indices {
            let trim = &trims[trim_index];
            let (this_interval, other_row, other_interval) = if trim.row_a == row {
                (trim.interval_a, trim.row_b, trim.interval_b)
            } else {
                (trim.interval_b, trim.row_a, trim.interval_a)
            };
            let span_this = this_interval.1 - this_interval.0;
            let span_other = other_interval.1 - other_interval.0;
            if span_this.abs() < CHAIN_PARAM_TOL || span_other.abs() < CHAIN_PARAM_TOL {
                return Err("chain_trim_interval_degenerate");
            }
            // Trim correspondence: fraction f on this row pairs
            // with 1 - f on the other, so
            // t_this(t_other) = i0 + s_this * (1 - (t_other - o0)/s_other).
            let scale = -span_this / span_other;
            let offset = this_interval.0 + span_this * (1.0 + other_interval.0 / span_other);
            let other_embedding = RowChainEmbedding {
                scale: embedding.scale * scale,
                offset: embedding.scale * offset + embedding.offset,
            };
            record_trim_nodes(trim, visited, &mut node_params, &mut period)?;
            if let Some(existing) = visited.get(&other_row) {
                consistency_or_period(*existing, other_embedding, &mut period)?;
                continue;
            }
            visited.insert(other_row, other_embedding);
            rows.push((other_row, other_embedding));
            queue.push_back(other_row);
        }
    }
    let mut nodes: Vec<(N, f64)> = node_params.into_iter().collect();
    nodes.sort_by(|a, b| a.1.total_cmp(&b.1));
    Ok(EmbeddedChain {
        rows,
        nodes,
        period,
    })
}

/// Register the trim's two endpoint nodes under the chain
/// parameter of its side-A interval ends (side B agrees by the
/// f ↔ 1 − f correspondence). A node re-registered at a
/// different parameter on a closed chain reveals the period.
fn record_trim_nodes<K, N>(
    trim: &ChainTrim<K, N>,
    visited: &HashMap<K, RowChainEmbedding>,
    node_params: &mut HashMap<N, f64>,
    period: &mut Option<f64>,
) -> Result<(), &'static str>
where
    K: Copy + Eq + Hash,
    N: Copy + Eq + Hash,
{
    let Some(embedding) = visited.get(&trim.row_a).copied().or_else(|| {
        visited.get(&trim.row_b).map(|embedding_b| {
            // Re-derive side A's embedding from side B (compose
            // E_b with the affine map t_b(t_a)) so node
            // parameters are always expressed via side A.
            let span_a = trim.interval_a.1 - trim.interval_a.0;
            let span_b = trim.interval_b.1 - trim.interval_b.0;
            let scale = -span_b / span_a;
            let offset = trim.interval_b.0 + span_b * (1.0 + trim.interval_a.0 / span_a);
            RowChainEmbedding {
                scale: embedding_b.scale * scale,
                offset: embedding_b.scale * offset + embedding_b.offset,
            }
        })
    }) else {
        return Ok(());
    };
    for (node, row_param) in [
        (trim.node_start, trim.interval_a.0),
        (trim.node_end, trim.interval_a.1),
    ] {
        let param = embedding.map(row_param);
        if let Some(&existing) = node_params.get(&node) {
            let delta = (existing - param).abs();
            if delta > CHAIN_PARAM_TOL {
                reconcile_period(delta, period)?;
            }
        } else {
            node_params.insert(node, param);
        }
    }
    Ok(())
}

fn consistency_or_period(
    existing: RowChainEmbedding,
    proposed: RowChainEmbedding,
    period: &mut Option<f64>,
) -> Result<(), &'static str> {
    if (existing.scale - proposed.scale).abs() > CHAIN_PARAM_TOL {
        return Err("chain_embedding_scale_inconsistent");
    }
    let delta = (existing.offset - proposed.offset).abs();
    if delta > CHAIN_PARAM_TOL {
        reconcile_period(delta, period)?;
    }
    Ok(())
}

fn reconcile_period(delta: f64, period: &mut Option<f64>) -> Result<(), &'static str> {
    match period {
        None => {
            *period = Some(delta);
            Ok(())
        }
        Some(existing) if (*existing - delta).abs() <= CHAIN_PARAM_TOL => Ok(()),
        Some(_) => Err("chain_period_inconsistent"),
    }
}

fn split_chain_at_clean_nodes<K, N>(
    chain: EmbeddedChain<K, N>,
) -> Result<Vec<SubChain<K, N>>, &'static str>
where
    K: Copy,
    N: Copy,
{
    let mut rows = chain.rows;
    let mut nodes = chain.nodes;
    if let Some(period) = chain.period {
        rebase_closed_chain(&mut rows, &mut nodes, period)?;
    }
    let mut clean_params = clean_node_params(&rows, &nodes);
    // The chain extremes are clean by construction (no row lies
    // beyond them); include them so rows extending past the
    // outermost trim node remain covered by a window.
    let extremes = rows.iter().fold((f64::MAX, f64::MIN), |acc, (_, e)| {
        let (lo, hi) = e.chain_interval();
        (acc.0.min(lo), acc.1.max(hi))
    });
    if extremes.0 < extremes.1 {
        clean_params.insert(0, extremes.0);
        clean_params.push(extremes.1);
        clean_params.sort_by(f64::total_cmp);
        clean_params.dedup_by(|a, b| (*a - *b).abs() <= CHAIN_PARAM_TOL);
    }
    if clean_params.len() < 2 {
        return Err("chain_without_clean_node_pair");
    }
    let mut sub_chains = Vec::with_capacity(clean_params.len() - 1);
    for window in clean_params.windows(2) {
        let (lo, hi) = (window[0], window[1]);
        let sub_rows = rows
            .iter()
            .copied()
            .filter(|(_, embedding)| {
                let (row_lo, row_hi) = embedding.chain_interval();
                row_lo >= lo - CHAIN_PARAM_TOL && row_hi <= hi + CHAIN_PARAM_TOL
            })
            .collect::<Vec<_>>();
        if sub_rows.is_empty() {
            continue;
        }
        let sub_nodes = nodes
            .iter()
            .copied()
            .filter(|&(_, param)| param >= lo - CHAIN_PARAM_TOL && param <= hi + CHAIN_PARAM_TOL)
            .collect::<Vec<_>>();
        sub_chains.push(SubChain {
            rows: sub_rows,
            nodes: sub_nodes,
        });
    }
    if sub_chains.is_empty() {
        return Err("chain_clean_split_produced_no_rows");
    }
    Ok(sub_chains)
}

/// Shift every row and node of a closed chain into one period
/// window starting at a clean node, so the seam falls where no
/// row spans it.
fn rebase_closed_chain<K, N: Copy>(
    rows: &mut [(K, RowChainEmbedding)],
    nodes: &mut Vec<(N, f64)>,
    period: f64,
) -> Result<(), &'static str> {
    if period <= CHAIN_PARAM_TOL {
        return Err("chain_period_degenerate");
    }
    let wrap = |param: f64| -> f64 { ((param % period) + period) % period };
    // Wrap rows by their interval start so each occupies one
    // contiguous image; find a cut where no wrapped row spans.
    let candidate_cuts: Vec<f64> = nodes.iter().map(|&(_, param)| wrap(param)).collect();
    let cut = candidate_cuts
        .iter()
        .copied()
        .find(|&cut| {
            rows.iter().all(|(_, embedding)| {
                let (lo, hi) = embedding.chain_interval();
                let lo_w = wrap(lo - cut);
                let hi_w = lo_w + (hi - lo);
                hi_w <= period + CHAIN_PARAM_TOL
            })
        })
        .ok_or("closed_chain_without_clean_cut")?;
    for (_, embedding) in rows.iter_mut() {
        let (lo, _) = embedding.chain_interval();
        // Shift so the row's image starts at wrap(lo - cut),
        // i.e. all rows and nodes share the frame [0, period).
        let shift = wrap(lo - cut) - lo;
        embedding.offset += shift;
    }
    for (_, param) in nodes.iter_mut() {
        *param = wrap(*param - cut);
    }
    nodes.sort_by(|a, b| a.1.total_cmp(&b.1));
    // Duplicate the cut node at the far end so the final
    // sub-chain window closes back onto it.
    if let Some(&first) = nodes.first() {
        nodes.push((first.0, first.1 + period));
    }
    Ok(())
}

fn clean_node_params<K, N>(rows: &[(K, RowChainEmbedding)], nodes: &[(N, f64)]) -> Vec<f64> {
    let mut params: Vec<f64> = nodes.iter().map(|&(_, param)| param).collect();
    params.sort_by(f64::total_cmp);
    params.dedup_by(|a, b| (*a - *b).abs() <= CHAIN_PARAM_TOL);
    params
        .into_iter()
        .filter(|&param| {
            rows.iter().all(|(_, embedding)| {
                let (lo, hi) = embedding.chain_interval();
                !(lo + CHAIN_PARAM_TOL < param && param < hi - CHAIN_PARAM_TOL)
            })
        })
        .collect()
}

#[cfg(test)]
#[path = "kp20_polar_support_cell_chain_discovery_tests.rs"]
mod tests;
