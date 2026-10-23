//! Consensus flooding and overwrap repair for periodic chart
//! anchoring (split from `chart_branches.rs` under the 500-line file
//! budget). See [`anchor_components`].

use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};

use super::super::chart_subdivision::FaceChartSubdivision;

fn u_at(
    sub: &FaceChartSubdivision,
    h: usize,
    n: usize,
    ends: &dyn Fn(usize) -> (usize, usize),
) -> f64 {
    let (s, _) = ends(h);
    let emb = &sub.embeddings[&h];
    if s == n {
        emb.uv[0].0
    } else {
        emb.uv[emb.uv.len() - 1].0
    }
}

/// Anchor value for the BFS translation. A CLOSED self-loop (both
/// ends on one node — a full pierce ring closing at its hit node)
/// has TWO u values at that node, one period apart, so endpoint
/// matching is ill-defined: aligning its twin start-to-start
/// stacks the pair over TWO periods and the overwrap check
/// declines the face (measured: the cross-drill's second pierce
/// loop). Self-loops anchor by their span MIDPOINT instead — the
/// twin's reversed unwrap then lands on the same period cell.
fn anchor_u(
    sub: &FaceChartSubdivision,
    h: usize,
    n: usize,
    ends: &dyn Fn(usize) -> (usize, usize),
) -> f64 {
    let (s, e) = ends(h);
    if s == e {
        let emb = &sub.embeddings[&h];
        let (lo, hi) = emb
            .uv
            .iter()
            .fold((f64::MAX, f64::MIN), |(lo, hi), &(u, _)| {
                (lo.min(u), hi.max(u))
            });
        0.5 * (lo + hi)
    } else {
        u_at(sub, h, n, ends)
    }
}

/// Consensus flooding against a given barrier set. Records the
/// period shifts it applies so an overwrapped attempt can be undone
/// and re-run with more barriers (see [`anchor_components`]).
fn run_flood(
    sub: &mut FaceChartSubdivision,
    interior: &[usize],
    barrier: &BTreeSet<usize>,
    period: f64,
    ends: &dyn Fn(usize) -> (usize, usize),
) -> (Vec<Vec<usize>>, HashMap<usize, f64>) {
    let mut node_to_hes: HashMap<usize, Vec<usize>> = HashMap::new();
    for &h in interior {
        let (s, e) = ends(h);
        for n in [s, e] {
            // Pole nodes are excluded like seam nodes: their
            // endpoint u is the incident meridian's azimuth, not a
            // shared chart coordinate — translating a neighbour to
            // match it would smear the frame (the boundary-first
            // unwrap learned the same lesson with degenerate
            // anchors).
            if !barrier.contains(&n) && !sub.poles.contains_key(&n) {
                node_to_hes.entry(n).or_default().push(h);
            }
        }
    }
    let mut applied: HashMap<usize, f64> = HashMap::new();
    let mut comp_of: HashMap<usize, usize> = HashMap::new();
    let mut components: Vec<Vec<usize>> = Vec::new();
    for &seed in interior {
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
                if barrier.contains(&n) {
                    continue;
                }
                for &g in node_to_hes.get(&n).into_iter().flatten() {
                    if comp_of.contains_key(&g) {
                        continue;
                    }
                    let delta = anchor_u(sub, h, n, ends) - anchor_u(sub, g, n, ends);
                    let k = (delta / period).round();
                    if k != 0.0
                        && let Some(emb) = sub.embeddings.get_mut(&g)
                    {
                        emb.shift_u(k * period);
                        *applied.entry(g).or_insert(0.0) += k * period;
                    }
                    comp_of.insert(g, comp_id);
                    members.push(g);
                    queue.push_back(g);
                }
            }
        }
        components.push(members);
    }
    (components, applied)
}

/// Index of the reference component (largest; ties: lowest index).
pub(super) fn ref_comp_of(components: &[Vec<usize>]) -> Option<usize> {
    components
        .iter()
        .enumerate()
        .max_by_key(|(i, m)| (m.len(), usize::MAX - i))
        .map(|(i, _)| i)
}

/// Reference-frame u extent of the largest component.
fn ref_span(sub: &FaceChartSubdivision, components: &[Vec<usize>]) -> (f64, f64) {
    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    if let Some(rc) = ref_comp_of(components) {
        for &h in &components[rc] {
            for &(u, _) in &sub.embeddings[&h].uv {
                lo = lo.min(u);
                hi = hi.max(u);
            }
        }
    }
    (lo, hi)
}

/// Consensus-anchor the interior half-edges: flood with the seam-node
/// barrier, and on OVERWRAP retry once with period-lattice junction
/// barriers added. The chart seam is a COORDINATE line even where no
/// seam edge crosses a junction, and the record splitter terminates
/// seam-crossing records exactly on it — nodes whose every incident
/// endpoint sits ON the period lattice are seam junctions in all but
/// topology. The tuned flood translates straight through them, and a
/// non-contractible groove spiral on a tip-chamfer cone unrolled to
/// 1.5 periods (PK-460: I2 declined the face to the 3D-Newell
/// fallback and the cut stranded 10 faces). Scoped as a RETRY: every
/// face the tuned flood anchors within one period keeps its shifts
/// byte-identically.
///
/// Returns the components, the effective barrier node set, the
/// period-lattice junction nodes (computed unconditionally — the
/// branch clustering merges their border copies whether or not the
/// overwrap retry ran, PK-464), and the reference component's
/// (lo, hi, span).
#[allow(clippy::type_complexity)]
pub(super) fn anchor_components(
    sub: &mut FaceChartSubdivision,
    interior: &[usize],
    seam_nodes: BTreeSet<usize>,
    period: f64,
    ends: &dyn Fn(usize) -> (usize, usize),
    u0: f64,
) -> (
    Vec<Vec<usize>>,
    BTreeSet<usize>,
    BTreeSet<usize>,
    (f64, f64, f64),
) {
    let (mut components, applied) = run_flood(sub, interior, &seam_nodes, period, ends);
    let (mut lo, mut hi) = ref_span(sub, &components);
    let mut span = if hi > lo { hi - lo } else { 0.0 };
    let mut seam_nodes = seam_nodes;
    let lattice_nodes: BTreeSet<usize> = {
        let mut node_us: BTreeMap<usize, Vec<f64>> = BTreeMap::new();
        for &h in interior {
            let (s, e) = ends(h);
            for n in [s, e] {
                node_us.entry(n).or_default().push(u_at(sub, h, n, ends));
            }
        }
        // Lattice epsilon covers record-endpoint fit slop: a crest
        // circle's pcurve endpoint landed 1.2e-3 off the period
        // lattice (M5 countersunk tip chamfer) and the junction
        // failed the exact test — no barrier retry, no wrap-node
        // merge, ChartOverwrap 7.85 declined the cone. Real interior
        // nodes sit far off the lattice; the ALL-endpoints
        // requirement keeps the set surgical.
        let eps = period * 5.0e-4;
        // The lattice is the face's own seam line, not u = 0: a part
        // placed rotated carries its seam (and the thread builder's
        // per-turn helix splits) at u = pi (measured, PK-919: the angled
        // fitting's cap screw root found 1 lattice barrier of its 14 and
        // kept a 41.4 span).
        node_us
            .iter()
            .filter(|(_, us)| {
                us.len() > 1
                    && us.iter().all(|u| {
                        let w = u - u0;
                        (w - period * (w / period).round()).abs() < eps
                    })
            })
            .map(|(n, _)| *n)
            .collect()
    };
    if span > period * (1.0 + 1e-6) && popkern_core::trace::trace_enabled("branch-pairing") {
        // Decline observability: which half-edges stretch the
        // reference frame past one period (the frame the ChartOverwrap
        // decline will cite). Printed before the lattice retry so both
        // the tuned and retried frames are attributable.
        if let Some(rc) = ref_comp_of(&components) {
            for &h in &components[rc] {
                let emb = &sub.embeddings[&h];
                let (blo, bhi) = emb
                    .uv
                    .iter()
                    .fold((f64::MAX, f64::MIN), |(l, hh), &(u, _)| {
                        (l.min(u), hh.max(u))
                    });
                let (s, e) = ends(h);
                eprintln!(
                    "  [branch] overwrap member he {h}: u [{blo:.4},{bhi:.4}] nodes ({s},{e}) uv-samples {:?}",
                    sub.embeddings[&h]
                        .uv
                        .iter()
                        .step_by((sub.embeddings[&h].uv.len() / 8).max(1))
                        .map(|&(u, v)| (
                            (u * 1000.0).round() / 1000.0,
                            (v * 1000.0).round() / 1000.0
                        ))
                        .collect::<Vec<_>>()
                );
            }
        }
    }
    if span > period * (1.0 + 1e-6) && !lattice_nodes.is_empty() {
        {
            for (&h, &shift) in &applied {
                if let Some(emb) = sub.embeddings.get_mut(&h) {
                    emb.shift_u(-shift);
                }
            }
            let widened: BTreeSet<usize> = seam_nodes.union(&lattice_nodes).copied().collect();
            let (retry_components, _) = run_flood(sub, interior, &widened, period, ends);
            let (retry_lo, retry_hi) = ref_span(sub, &retry_components);
            let retry_span = if retry_hi > retry_lo {
                retry_hi - retry_lo
            } else {
                0.0
            };
            if popkern_core::trace::trace_enabled("branch-pairing") {
                eprintln!(
                    "  [branch] overwrap {span:.4}: retry with {} lattice barrier(s) -> {retry_span:.4}",
                    lattice_nodes.len()
                );
            }
            components = retry_components;
            span = retry_span;
            lo = retry_lo;
            hi = retry_hi;
            seam_nodes = widened;
        }
    }
    (components, seam_nodes, lattice_nodes, (lo, hi, span))
}
