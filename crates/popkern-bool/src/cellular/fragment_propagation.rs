//! A-side fragment verdict propagation planner (PK-496).
//!
//! Thread-class booleans put O(turns) wrap-band fragments on ONE
//! workpiece face (measured: 43 of 61 fragments at 21 turns), and
//! classifying each with its own point-in-solid rays against the
//! multi-wrap tool walls priced the merge stage quadratically
//! (31 s / 100 s / 391 s at 21/41/81 turns). Crossing operand A's own
//! boundary edge cannot change the side of B, and crossing a
//! TRANSVERSAL intersection record flips it — so one seed ray per
//! connected component, flip parity from BFS, and a ray audit on the
//! deepest members classify the whole component. B-side fragments
//! keep the per-fragment path: every B face carries O(1) fragments,
//! and PK-428 walk-phantom semantics live in B-side ray verdicts
//! (the same reason PK-494 scoped witness overrides A-only).
//!
//! Every gate failure (graze record, residue record, parity cycle
//! conflict, seed or audit surprise) falls back to the per-fragment
//! path for the affected component — the plan only ever REMOVES ray
//! casts whose outcome the audit confirms, never a verdict class.
//! Kill switch: `POPKERN_CELLULAR_PROPAGATE=0` (per-call read, PK-239
//! precedent). Design: docs/decisions/
//! pk496-cellular-fragment-verdict-propagation.md.

use std::collections::BTreeMap;

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::BRep;

use crate::merge::graph::{EdgeSource, IntersectionGraph, Operand};
use crate::merge::walk::Fragment;

/// Minimum component size worth propagating: seed + two audits are
/// ray-classified anyway, so below this the plan saves nothing.
const MIN_COMPONENT: usize = 4;
/// Minimum |sin| between the two faces' normals at a record midpoint
/// for the record to count as a transversal crossing of operand B's
/// boundary. Below it the contact is graze/tangent structure (the
/// near-tangent 1/sin family, PK-463..465): no propagation edge —
/// the fragments stay individually sampled via other paths or not
/// at all.
const MIN_TRANSVERSAL_SIN: f64 = 0.05;
/// Ray-audited members per component: the BFS-deepest members carry
/// the most accumulated flips, so they maximize the chance of
/// catching a wrong flip edge.
const AUDITS_PER_COMPONENT: usize = 2;

/// One propagation component over eligible operand-A fragments.
pub(super) struct Component {
    /// Fragment indices with their flip parity from the seed
    /// (`false` = same verdict as seed, `true` = opposite).
    pub members: Vec<(usize, bool)>,
    /// The ray-classified representative (largest boundary loop).
    pub seed: usize,
    /// Ray-audited members: (fragment index, expected parity).
    pub audits: Vec<(usize, bool)>,
}

/// The plan: which fragments defer to which component.
pub(super) struct PropagationPlan {
    /// Per fragment: `Some(component index)` when the verdict comes
    /// from propagation (a non-seed, non-audit member).
    pub deferred: Vec<Option<usize>>,
    pub components: Vec<Component>,
}

/// Whether propagation is enabled (`POPKERN_CELLULAR_PROPAGATE=0`
/// pins the per-fragment path; per-call read).
pub(super) fn propagation_enabled() -> bool {
    std::env::var_os("POPKERN_CELLULAR_PROPAGATE")
        .map(|v| v != "0")
        .unwrap_or(true)
}

/// Build the A-side propagation plan. `None` means no component is
/// worth propagating (every fragment classifies per-fragment).
pub(super) fn build_plan(
    graph: &IntersectionGraph,
    fragments: &[Fragment],
    degenerate: &[bool],
    overlap_a: &std::collections::HashSet<popkern_core::topology::FaceHandle>,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Option<PropagationPlan> {
    let n = fragments.len();
    // Eligibility mirrors the offset-pair routing in the classifier:
    // fragments that would take the OnBoundary/coincidence machinery
    // never propagate.
    let eligible: Vec<bool> = fragments
        .iter()
        .enumerate()
        .map(|(i, f)| {
            f.face.operand == Operand::A
                && !degenerate[i]
                && !overlap_a.contains(&f.face.face)
                && f.inner_loops.is_empty()
                && !f.outer_loop.iter().any(|&h| {
                    let he = &graph.half_edges[h];
                    he.origin == graph.half_edges[he.twin].origin
                })
        })
        .collect();

    // Half-edge -> eligible fragment map (ORDERED: the planner's
    // iteration order decides seed and audit choice, which decides
    // which fragments are ray-sampled and therefore which witness
    // votes exist — a HashMap here made the whole boolean's
    // serve/decline outcome vary run-to-run, PK-512 measured: the
    // same len-110 cut declined non-manifold three times, then
    // served). A half-edge claimed by two eligible fragments would
    // mean a walk defect; abort the plan.
    let mut owner: BTreeMap<usize, usize> = BTreeMap::new();
    for (i, f) in fragments.iter().enumerate() {
        if !eligible[i] {
            continue;
        }
        for &h in &f.outer_loop {
            if owner.insert(h, i).is_some() {
                return None;
            }
        }
    }

    // Adjacency with flip marks.
    let mut adj: Vec<Vec<(usize, bool)>> = vec![Vec::new(); n];
    for (&h, &f1) in &owner {
        let t = graph.half_edges[h].twin;
        if t < h {
            continue; // visit each pair once
        }
        let Some(&f2) = owner.get(&t) else { continue };
        if f1 == f2 {
            continue;
        }
        let flip = match &graph.half_edges[h].source {
            EdgeSource::OriginalA(_) | EdgeSource::SplitA { .. } => false,
            EdgeSource::Shared { face_a, face_b } => {
                if !record_is_transversal(graph, h, t, *face_a, *face_b, brep_a, brep_b) {
                    continue; // graze / residue / unevaluable: no edge
                }
                true
            }
            _ => continue,
        };
        adj[f1].push((f2, flip));
        adj[f2].push((f1, flip));
    }

    // BFS components with parity; a parity conflict (odd flip cycle —
    // e.g. a graze record misread as transversal) poisons the whole
    // component back to per-fragment sampling.
    let mut comp_of: Vec<Option<usize>> = vec![None; n];
    let mut parity: Vec<bool> = vec![false; n];
    let mut depth: Vec<usize> = vec![0; n];
    let mut components: Vec<Component> = Vec::new();
    let mut deferred: Vec<Option<usize>> = vec![None; n];
    for start in 0..n {
        if !eligible[start] || comp_of[start].is_some() {
            continue;
        }
        let ci = components.len();
        let mut members: Vec<(usize, bool)> = Vec::new();
        let mut queue = std::collections::VecDeque::new();
        comp_of[start] = Some(ci);
        parity[start] = false;
        depth[start] = 0;
        queue.push_back(start);
        let mut consistent = true;
        while let Some(u) = queue.pop_front() {
            members.push((u, parity[u]));
            for &(v, flip) in &adj[u] {
                let want = parity[u] ^ flip;
                match comp_of[v] {
                    None => {
                        comp_of[v] = Some(ci);
                        parity[v] = want;
                        depth[v] = depth[u] + 1;
                        queue.push_back(v);
                    }
                    Some(c) if c == ci => {
                        if parity[v] != want {
                            consistent = false;
                        }
                    }
                    Some(_) => consistent = false,
                }
            }
        }
        if !consistent || members.len() < MIN_COMPONENT {
            if popkern_core::trace::trace_enabled("cellular") && members.len() >= MIN_COMPONENT {
                eprintln!(
                    "  [propagate] component of {} INCONSISTENT parity — per-fragment",
                    members.len()
                );
            }
            continue; // fragments stay per-fragment (deferred = None)
        }
        // Seed: largest boundary loop (most robust sample). Audits:
        // deepest members (most accumulated flips), excluding seed.
        let seed = members
            .iter()
            .map(|&(m, _)| m)
            .max_by_key(|&m| (fragments[m].outer_loop.len(), std::cmp::Reverse(m)))?;
        let mut by_depth: Vec<usize> = members
            .iter()
            .map(|&(m, _)| m)
            .filter(|&m| m != seed)
            .collect();
        by_depth.sort_by_key(|&m| std::cmp::Reverse(depth[m]));
        let seed_parity = parity[seed];
        let audits: Vec<(usize, bool)> = by_depth
            .iter()
            .take(AUDITS_PER_COMPONENT)
            .map(|&m| (m, parity[m] ^ seed_parity))
            .collect();
        for &(m, _) in &members {
            if m != seed && !audits.iter().any(|&(a, _)| a == m) {
                deferred[m] = Some(components.len());
            }
        }
        // Re-normalize member parity relative to the seed.
        let members = members
            .into_iter()
            .map(|(m, p)| (m, p ^ seed_parity))
            .collect();
        components.push(Component {
            members,
            seed,
            audits,
        });
    }
    if components.is_empty() {
        return None;
    }
    Some(PropagationPlan {
        deferred,
        components,
    })
}

/// Whether a `Shared` record half-edge pair is a TRANSVERSAL crossing
/// of the two faces: normals evaluated at each side's pcurve midpoint
/// must meet at more than the graze angle. Residue records (both
/// sides carry the owner face, PK-465) and pairs missing a pcurve or
/// surface decline conservatively.
fn record_is_transversal(
    graph: &IntersectionGraph,
    h: usize,
    t: usize,
    face_a: popkern_core::topology::FaceHandle,
    face_b: popkern_core::topology::FaceHandle,
    brep_a: &BRep,
    brep_b: &BRep,
) -> bool {
    let (he_h, he_t) = (&graph.half_edges[h], &graph.half_edges[t]);
    // Normal records carry face_a on one half-edge and face_b (an
    // operand-B ref) on the other; a residue pair carries the owner
    // on BOTH. Anything else is not the shape this model covers.
    let (a_side, b_side) = match (he_h.face, he_t.face) {
        (Some(fa), Some(fb)) if fa.operand == Operand::A && fb.operand == Operand::B => {
            (he_h, he_t)
        }
        (Some(fa), Some(fb)) if fa.operand == Operand::B && fb.operand == Operand::A => {
            (he_t, he_h)
        }
        _ => return false,
    };
    let normal_at_mid = |he: &crate::merge::graph::GraphHalfEdge, brep: &BRep, face| {
        let pc = he.pcurve.as_ref()?;
        let (t0, t1) = pc.parameter_range();
        let uv = pc.point_at(0.5 * (t0 + t1));
        let surface = brep.faces.get(face)?.surface.as_ref()?;
        surface.geometry.normal_at(uv.x, uv.y).normalized()
    };
    let Some(na) = normal_at_mid(a_side, brep_a, face_a) else {
        return false;
    };
    let Some(nb) = normal_at_mid(b_side, brep_b, face_b) else {
        return false;
    };
    na.cross(&nb).length() >= MIN_TRANSVERSAL_SIN
}
