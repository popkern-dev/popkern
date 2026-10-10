//! Period-branch resolution for chart subdivisions (PK-104).
//!
//! Split from `chart_subdivision.rs` to honour the 500-line budget:
//! consensus anchoring of embedding components, deterministic seam
//! border assignment (material-side rule), the bounded floating-
//! component search, per-node branch clustering, and the boundary
//! shoelace σ calibration.

use std::collections::{BTreeMap, BTreeSet};

#[path = "chart_branch_anchoring.rs"]
mod anchoring;
#[path = "chart_branch_clustering.rs"]
pub(crate) mod clustering;
use anchoring::ref_comp_of;

use clustering::boundary_shoelace_sigma;

use super::chart_subdivision::{FaceChartSubdivision, SubdivisionIssue};
use crate::merge::graph::IntersectionGraph;

/// Periodic charts: consensus anchoring, seam border assignment with a
/// bounded balance search, then per-node u-clustering into branches.
pub(super) fn anchor_periodic(
    graph: &IntersectionGraph,
    face: &crate::merge::graph::FaceRef,
    sub: &mut FaceChartSubdivision,
    period: f64,
    ends: &dyn Fn(usize) -> (usize, usize),
    ring_sigma: Option<f64>,
) {
    let hes: Vec<usize> = sub.embeddings.keys().copied().collect();
    // Seam candidates are OPERAND edges whose two sides both live on
    // this face (a periodic wall's seam). Trim-RESIDUE record pairs
    // are also both-sides-on-one-face (single-face registration, the
    // same seam-like shape) but are interior slits, not border
    // copies: fed to the seam border assignment they get one side
    // bumped a whole period (measured, PK-465: the angled tap's
    // mouth-adjacent owner-A residue embedded at u 6.6-7.8 on a
    // [0,6.28] wall strip, its branch clusters split 1-in/0-out, and
    // the wall declined to the 3D-Newell fallback). Shared sources
    // are excluded.
    // A Shared pair that took over the face's own seam (the
    // duplicate-cover pass removed the operand seam in its favour) IS
    // the seam: PK-750's flush box removed a rod's seam for a cut line
    // lying on it, and the chart, seeing no seam, anchored the rim
    // record a period off.
    let seam_hes: BTreeSet<usize> = hes
        .iter()
        .copied()
        .filter(|&h| {
            graph.seam_substitutes.contains(&(*face, h))
                || (sub.embeddings[&h].seam_like
                    && !matches!(
                        graph.half_edges[h].source,
                        crate::merge::graph::EdgeSource::Shared { .. }
                    ))
        })
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
    let true_seam_nodes = seam_nodes.clone();
    let (components, seam_nodes, lattice_nodes, (lo, hi, span)) =
        anchoring::anchor_components(sub, &interior, seam_nodes, period, ends);
    // Period-lattice junctions (whether or not the overwrap retry
    // barriered them), minus the true seam nodes: their border copies
    // are one azimuth and cluster as one branch (see
    // `cluster_branches`). Computed unconditionally — a small
    // residual band whose flood never overwraps still carries
    // lattice junctions (PK-464: the crest band above the top
    // chamfer ring declined 1-in/0-out at both its per-turn
    // junctions).
    let wrap_nodes: BTreeSet<usize> = seam_nodes
        .difference(&true_seam_nodes)
        .copied()
        .chain(
            lattice_nodes
                .iter()
                .copied()
                .filter(|n| !true_seam_nodes.contains(n)),
        )
        .collect();
    if span > period * (1.0 + 1e-6) {
        sub.report
            .issues
            .push(SubdivisionIssue::ChartOverwrap { span, period });
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
    // re-created the 3D-Newell 16-half-edge mega-loop orbit.)
    let mut seam_pairs: Vec<(usize, usize)> = Vec::new();
    for &h in &seam_hes {
        let tw = graph.half_edges[h].twin;
        if h < tw && seam_hes.contains(&tw) {
            seam_pairs.push((h, tw));
        }
    }
    // Window the low border must land in: [hi - period, lo]. Without
    // an interior reference, anchor on the first seam he's own u.
    // Window epsilon is CHART-scale, not absolute: a full-period
    // reference ring (a pierce loop wrapping the whole wall) pins the
    // window to near-zero width, and a fitted record's endpoint u
    // carries fit-level noise (~1e-8) that a hard 1e-9 rejects
    // (measured: the cross-drill wall's seam declined as a misfit
    // against its own ring's cut azimuth).
    let eps_w = period * 1e-6;
    let (window_lo_min, window_lo_max) = if hi > lo {
        (hi - period - eps_w, lo + eps_w)
    } else {
        let u0 = seam_pairs
            .first()
            .map(|&(h, _)| sub.embeddings[&h].uv[0].0)
            .unwrap_or(0.0);
        (u0 - eps_w, u0 + eps_w)
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
    let seam_trace = popkern_core::trace::trace_enabled("seam");
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
    // The stored trim ring's winding outranks the registered-run sum
    // (see `stored_ring_sigma`).
    let shoelace = ring_sigma.or_else(|| boundary_shoelace_sigma(graph, sub));
    if seam_trace {
        eprintln!("  [seam] shoelace {shoelace:?} vs provisional sigma {sigma:+.0}");
    }
    if let Some(s) = shoelace
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
    // Domain strip: after seam assignment the face's chart occupies
    // [seam_low, seam_low + period] (the low seam copies sit at the
    // strip's low border). A floating component must land INSIDE this
    // strip — the low-border WINDOW is the wrong target (it is the
    // set of possible border positions, which extends a full period
    // below the strip). Placement inside the strip is geometric and
    // unique for any sub-period component, and balance cannot decide
    // it: on the breached-bore drill wall BOTH placements of the
    // wall×side record arc balanced, but the off-strip one (u < 0)
    // corrupted the rotational order at its two seam nodes and every
    // orbit collided (measured: tilted_hole_plate, 822 open edges).
    let seam_low: Option<f64> = if seam_pairs.is_empty() {
        None
    } else {
        seam_pairs
            .iter()
            .flat_map(|&(h, tw)| [h, tw])
            .map(|he| {
                let emb = &sub.embeddings[&he];
                emb.uv.iter().map(|&(u, _)| u).fold(f64::MAX, f64::min)
            })
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
    };
    let mut floating: Vec<usize> = Vec::new();
    let ref_comp = ref_comp_of(&components);
    for ci in (0..components.len()).filter(|&i| Some(i) != ref_comp) {
        let (mut cmin, mut cmax) = (f64::MAX, f64::MIN);
        for h in &components[ci] {
            for &(u, _) in &sub.embeddings[h].uv {
                cmin = cmin.min(u);
                cmax = cmax.max(u);
            }
        }
        let strip_shift = seam_low.and_then(|u_low| {
            if cmax - cmin >= period * (1.0 - 1e-9) {
                return None; // full-period ring: no freedom, leave as anchored
            }
            let k = ((u_low - cmin) / period - 1e-9).ceil();
            let nmax = cmax + k * period;
            // Only shift when the whole component fits the strip; a
            // border-straddling component keeps the search below.
            // CHART-SCALE fit epsilon: an absolute 1e-9 rejected a
            // component whose endpoint sits EXACTLY on the high
            // border (u = 2π + fp jitter ~2e-9 under exact-arm
            // seeding); the low_shift fallback then dropped the whole
            // chain a full period low, its records clustered with the
            // OPPOSITE border's branch, and the wedge pairing walked
            // it as a degenerate out-and-back region (measured:
            // torus × cylinder wall under the exact-torus NURBS arm).
            (nmax <= u_low + period + period * 1.0e-6).then_some(k * period)
        });
        if popkern_core::trace::trace_enabled("branch-pairing") {
            eprintln!(
                "  [comp] ci {ci} hes {:?} span [{cmin:.3},{cmax:.3}] seam_low {seam_low:?} strip_shift {strip_shift:?}",
                components[ci]
            );
        }
        let applied = match strip_shift {
            Some(shift) => {
                if shift != 0.0 {
                    for &h in &components[ci] {
                        if let Some(emb) = sub.embeddings.get_mut(&h) {
                            emb.shift_u(shift);
                        }
                    }
                    cmin += shift;
                    cmax += shift;
                }
                true
            }
            None => false,
        };
        if !applied
            && let Some(shift) = low_shift(cmin)
            && shift != 0.0
        {
            for &h in &components[ci] {
                if let Some(emb) = sub.embeddings.get_mut(&h) {
                    emb.shift_u(shift);
                }
            }
            cmin += shift;
            cmax += shift;
        }
        let _ = cmin;
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
    //
    // The search FACTORS: two floating components interact only
    // through SHARED SEAM NODES (balance is a per-node property and a
    // bump moves only its own component's entries), so the joint
    // exponential search decomposes into independent per-group
    // searches over seam-node-connected groups (measured: the angled
    // tap's residue-carrying bore wall floats more components than
    // the historical joint 2^4 budget while no connected group has
    // more than a couple — the joint search exhausted and the whole
    // wall fell to the 3D-Newell walker, stranding record sides,
    // PK-465).
    super::chart_branch_float::resolve_floating_and_cluster(
        sub,
        period,
        ends,
        &hes,
        &components,
        &floating,
        &seam_nodes,
        &true_seam_nodes,
        &wrap_nodes,
        &seam_hes,
    );
}
