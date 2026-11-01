//! Splitting a non-manifold selected boundary into manifold pieces
//! (plan §4.9, PK-987).
//!
//! **Edge pairing.** A kept face has exactly one selected side, so
//! around an edge the kept fins, in radial order, alternate between
//! opening and closing a selected wedge. The two kept fins bounding
//! one selected wedge are twins: one result edge per selected wedge,
//! m edges for 2m kept fins. The pairing follows from the radial order
//! alone, not from which fragment the assembler meets first (the
//! first-come pairing that mints a duplicate with a free twin).
//!
//! **Vertex fans.** At a vertex, each kept face contributes a corner
//! (the fin arriving and its successor leaving). Corners whose fins
//! are twins lie in one fan, and each fan is one result vertex (Lee &
//! Lee's partial vertices): two bodies touching at a point keep a
//! vertex each.
//!
//! The plan names, in graph terms, every result edge that differs from
//! the assembler's one-edge-per-graph-edge default: each edge of a
//! non-manifold edge's pairs, and each edge at a split vertex (whose
//! ends pick the fan's vertex). A pair whose two fins lift different
//! graph edges (two coincident graph edges the join made one) becomes
//! one result edge as well.

use std::collections::BTreeMap;

use super::ids::{EdgeId, FinId, VertexId};
use super::model::{CellComplex, SideIndex};
use super::radial::RadialSides;
use super::regions::Forest;

/// One result edge the plan builds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PlannedEdge {
    /// The graph half-edge whose curve the edge takes (the opening
    /// fin's).
    pub curve_half_edge: usize,
    /// (graph node, vertex copy) at the complex edge's start.
    pub start: (usize, usize),
    /// (graph node, vertex copy) at its end.
    pub end: (usize, usize),
}

/// The split of a non-manifold selected boundary.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct SplitPlan {
    /// The result edges the plan builds.
    pub edges: Vec<PlannedEdge>,
    /// (fragment, loop half-edge as the fragment lists it) → (planned
    /// edge, whether that half-edge runs start → end).
    pub uses: BTreeMap<(usize, usize), (usize, bool)>,
    /// Vertex copies per graph node split into fans (two or more).
    pub vertex_copies: BTreeMap<usize, usize>,
}

/// Why a selected boundary cannot be split.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum SplitDecline {
    /// The kept fins around an edge do not alternate between opening
    /// and closing selected wedges.
    UnpairedFins {
        /// The edge.
        edge: EdgeId,
    },
    /// A kept fin or its vertex has no graph source to build from.
    Unmapped {
        /// The fin.
        fin: FinId,
    },
}

/// The twins of every edge: pairs (opening fin, closing fin) of kept
/// fins bounding one selected wedge, in radial order.
fn edge_pairs(
    c: &CellComplex,
    sides: &RadialSides,
    kept_face: &[bool],
    selected_shell: &[bool],
) -> Result<Vec<Vec<(FinId, FinId)>>, SplitDecline> {
    let kept = |f: FinId| kept_face[c.fins[f.index()].face.index()];
    // Whether the wedge after `f` (counter-clockwise) is selected.
    let opens = |f: FinId| {
        let fin = &c.fins[f.index()];
        let side = if sides.ahead_faces_next[f.index()] {
            SideIndex::Ahead
        } else {
            SideIndex::Behind
        };
        c.faces[fin.face.index()].sides[side as usize]
            .shell
            .is_some_and(|s| selected_shell[s.index()])
    };
    let mut out = Vec::with_capacity(c.edges.len());
    for (ei, e) in c.edges.iter().enumerate() {
        let ring: Vec<FinId> = e.radial.iter().copied().filter(|&f| kept(f)).collect();
        let mut pairs = Vec::new();
        if ring.len() % 2 == 1 {
            return Err(SplitDecline::UnpairedFins {
                edge: EdgeId::from_index(ei),
            });
        }
        for (i, &f) in ring.iter().enumerate() {
            if !opens(f) {
                continue;
            }
            let g = ring[(i + 1) % ring.len()];
            if opens(g) {
                return Err(SplitDecline::UnpairedFins {
                    edge: EdgeId::from_index(ei),
                });
            }
            pairs.push((f, g));
        }
        if pairs.len() * 2 != ring.len() {
            return Err(SplitDecline::UnpairedFins {
                edge: EdgeId::from_index(ei),
            });
        }
        out.push(pairs);
    }
    Ok(out)
}

/// Per vertex, the fan of each kept fin's corner there: (fin, vertex)
/// → fan index, and the fan count per vertex.
type Fans = (
    BTreeMap<(FinId, VertexId), usize>,
    BTreeMap<VertexId, usize>,
);

fn vertex_fans(c: &CellComplex, kept_face: &[bool], pairs: &[Vec<(FinId, FinId)>]) -> Fans {
    // Corners: (arriving fin, leaving fin) per vertex, numbered.
    let mut corners: Vec<(VertexId, FinId, FinId)> = Vec::new();
    let mut corner_of: BTreeMap<(FinId, VertexId), usize> = BTreeMap::new();
    for (fi, fin) in c.fins.iter().enumerate() {
        let f = FinId::from_index(fi);
        if !kept_face[fin.face.index()] {
            continue;
        }
        let Some(v) = c.fin_end(f) else {
            continue;
        };
        let k = corners.len();
        corners.push((v, f, fin.next));
        corner_of.insert((f, v), k);
        corner_of.entry((fin.next, v)).or_insert(k);
    }
    let mut forest = Forest::new(corners.len());
    // A fin's corner at its arriving vertex and the corner it leaves
    // from are both keyed (fin, vertex); a closed edge's fin reaches
    // one vertex both ways, and both corners join.
    for ring in pairs {
        for &(f, g) in ring {
            for v in [c.fin_start(f), c.fin_end(f)].into_iter().flatten() {
                let ends = |x: FinId| {
                    [
                        corner_of.get(&(x, v)).copied(),
                        corners
                            .iter()
                            .position(|&(cv, _, leave)| cv == v && leave == x),
                    ]
                };
                let mut at: Vec<usize> = ends(f).into_iter().chain(ends(g)).flatten().collect();
                at.dedup();
                for w in at.windows(2) {
                    forest.union(w[0], w[1]);
                }
            }
        }
    }
    let mut fan_index: BTreeMap<usize, usize> = BTreeMap::new();
    let mut fans_at: BTreeMap<VertexId, usize> = BTreeMap::new();
    for (k, &(v, _, _)) in corners.iter().enumerate() {
        let root = forest.find(k);
        if let std::collections::btree_map::Entry::Vacant(slot) = fan_index.entry(root) {
            let n = fans_at.entry(v).or_insert(0);
            slot.insert(*n);
            *n += 1;
        }
    }
    let mut fan_of = BTreeMap::new();
    for (&(f, v), &k) in &corner_of {
        fan_of.insert((f, v), fan_index[&forest.find(k)]);
    }
    for (k, &(v, _, leave)) in corners.iter().enumerate() {
        fan_of.insert((leave, v), fan_index[&forest.find(k)]);
    }
    (fan_of, fans_at)
}

/// The split plan of a selection (`kept_face`, `selected_shell`) over
/// a radially ordered complex.
pub(crate) fn split_plan(
    c: &CellComplex,
    sides: &RadialSides,
    kept_face: &[bool],
    selected_shell: &[bool],
) -> Result<SplitPlan, SplitDecline> {
    let pairs = edge_pairs(c, sides, kept_face, selected_shell)?;
    let (fan_of, fans_at) = vertex_fans(c, kept_face, &pairs);
    let split = |v: VertexId| fans_at.get(&v).is_some_and(|&n| n > 1);
    let mut plan = SplitPlan::default();
    for (&v, &n) in &fans_at {
        if n > 1
            && let Some(node) = c.vertices[v.index()].source_node
        {
            plan.vertex_copies.insert(node, n);
        }
    }
    for (ei, ring) in pairs.iter().enumerate() {
        let e = &c.edges[ei];
        if ring.len() < 2 && !split(e.start) && !split(e.end) {
            continue;
        }
        for &(f, g) in ring {
            let unmapped = |fin| SplitDecline::Unmapped { fin };
            let node = |v: VertexId, fin| c.vertices[v.index()].source_node.ok_or(unmapped(fin));
            let copy = |v: VertexId| fan_of.get(&(f, v)).copied().unwrap_or(0);
            let curve_half_edge = c.fins[f.index()].source_half_edge.ok_or(unmapped(f))?;
            let index = plan.edges.len();
            plan.edges.push(PlannedEdge {
                curve_half_edge,
                start: (node(e.start, f)?, copy(e.start)),
                end: (node(e.end, f)?, copy(e.end)),
            });
            for x in [f, g] {
                let fin = &c.fins[x.index()];
                let he = fin.source_half_edge.ok_or(unmapped(x))?;
                let fragment = *c.faces[fin.face.index()]
                    .source_fragments
                    .first()
                    .ok_or(unmapped(x))?;
                plan.uses.insert((fragment, he), (index, fin.forward));
            }
        }
    }
    Ok(plan)
}

#[cfg(test)]
#[path = "split_plan_tests.rs"]
mod tests;
