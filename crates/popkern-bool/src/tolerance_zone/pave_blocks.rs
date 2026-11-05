//! Pave blocks, their shrunk ranges and common blocks (plan §4.3,
//! §4.3a micro edge consequence).
//!
//! The resolved vertices cut each record and each operand edge they
//! lie on into pave blocks: the pieces between consecutive paves. A
//! block's shrunk range is the part of its parameter range outside both
//! its end vertices' tolerance balls: where the piece is a curve of its
//! own rather than part of a vertex. A block with an empty shrunk range
//! lies within its vertices (Jackson Fig. 2): an edge shorter than
//! tolerance between two vertices the resolution kept apart. Two
//! vertices whose balls the piece covers overlap, so they shared a
//! component and some rule separated them; serving the block would
//! serve an edge no consumer can tell from a point.
//!
//! Two blocks of different owners joining the same two vertices and
//! lying within those vertices' tolerances of each other along their
//! shrunk ranges are a common block: one edge of the result, which
//! served twice would bound a face piece no wider than tolerance
//! (Jackson Fig. 5's sliver). Two operand edges (E-E), a piece of a
//! record along a piece of an operand edge, or two records tracing one
//! piece. `common_blocks` finds only records lying along an edge for
//! their whole length; a record whose piece runs along an edge is found
//! here. A face piece that narrows to an edge between blocks of their
//! own is decided where faces are subdivided, not here.
//!
//! Measured (bool, examples, mesh, blend; table on): no micro block
//! across 1449 shadowed booleans; 33 common block pairs in 10 the
//! whole-record pass does not see, E-E (two cubes overlapping along
//! edges, flush tangent boxes) and record pieces along an operand edge
//! (a near-tangent cap casting, a notch rim).
//!
//! An operand edge is cut only where some candidate lies in its
//! interior; an edge with none is the operand's own, not a block of
//! the boolean.

use std::collections::BTreeMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, closest_point_in_window};
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};
use slotmap::Key;

use super::candidates::{CandidateVertex, Incidence};
use super::resolve::Resolution;
use crate::intersect::IntersectionRecord;

/// Samples along a block for its shrunk range and common-block test.
const SAMPLES: usize = 32;
/// Bisection steps refining a shrunk range end.
const REFINE: usize = 40;

/// What a pave block is a piece of.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum BlockOwner {
    /// A record.
    Record(usize),
    /// An operand edge, by operand and the edge's handle bits.
    Edge(usize, u64),
}

/// One pave block.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PaveBlock {
    /// What it is a piece of.
    pub owner: BlockOwner,
    /// Its parameter range on the owner's curve.
    pub range: (f64, f64),
    /// The vertices it runs between, `None` at an operand edge's end
    /// no candidate resolved.
    pub ends: [Option<usize>; 2],
    /// The part of `range` outside both end balls; `None` when the
    /// block lies within its vertices.
    pub shrunk: Option<(f64, f64)>,
}

/// The curve a block owner runs along, with its parameter range.
fn owner_curve<'a>(
    owner: BlockOwner,
    operands: [&'a BRep; 2],
    records: &'a [IntersectionRecord],
    edges: &BTreeMap<(usize, u64), EdgeHandle>,
) -> Option<(&'a dyn CurveEval, (f64, f64))> {
    match owner {
        BlockOwner::Record(r) => {
            let c = &records.get(r)?.curve.curve_3d;
            Some((c as &dyn CurveEval, c.parameter_range()))
        }
        BlockOwner::Edge(o, bits) => {
            let c = operands[o]
                .edges
                .get(*edges.get(&(o, bits))?)?
                .curve
                .as_ref()?;
            Some((
                &c.geometry as &dyn CurveEval,
                (c.t_start.min(c.t_end), c.t_start.max(c.t_end)),
            ))
        }
    }
}

/// Each owner's paves: (parameter, group), sorted. An operand edge's
/// end vertices pave it where a candidate resolved them.
fn paves(
    cands: &[CandidateVertex],
    res: &Resolution,
    operands: [&BRep; 2],
    edges: &BTreeMap<(usize, u64), EdgeHandle>,
) -> BTreeMap<BlockOwner, Vec<(f64, usize)>> {
    let mut out: BTreeMap<BlockOwner, Vec<(f64, usize)>> = BTreeMap::new();
    for (i, c) in cands.iter().enumerate() {
        for inc in &c.incidence {
            match *inc {
                Incidence::Record { record, s, .. } => {
                    out.entry(BlockOwner::Record(record))
                        .or_default()
                        .push((s, res.group_of[i]));
                }
                Incidence::Edge { operand, edge, t } => {
                    out.entry(BlockOwner::Edge(operand, edge.data().as_ffi()))
                        .or_default()
                        .push((t, res.group_of[i]));
                }
                _ => {}
            }
        }
    }
    // Operand edges' end vertices, where a candidate resolved them.
    let mut vertex_group: BTreeMap<(usize, u64), usize> = BTreeMap::new();
    for (i, c) in cands.iter().enumerate() {
        for inc in &c.incidence {
            if let Incidence::Vertex { operand, vertex } = *inc {
                vertex_group.insert((operand, vertex.data().as_ffi()), res.group_of[i]);
            }
        }
    }
    // An edge both of whose end vertices resolved runs between two
    // vertices of the boolean: a block whole, which a record running on
    // past one of those vertices lies along piece by piece.
    for (o, brep) in operands.iter().enumerate() {
        for (eh, e) in &brep.edges {
            let Some(h) = brep.half_edges.get(e.half_edge) else {
                continue;
            };
            let twin = h
                .mate()
                .and_then(|m| brep.half_edges.get(m))
                .map(|t| t.origin);
            let resolved = |v: VertexHandle| vertex_group.contains_key(&(o, v.data().as_ffi()));
            if resolved(h.origin) && twin.is_some_and(resolved) {
                out.entry(BlockOwner::Edge(o, eh.data().as_ffi()))
                    .or_default();
            }
        }
    }
    for (owner, list) in out.iter_mut() {
        let BlockOwner::Edge(o, bits) = *owner else {
            continue;
        };
        let Some(edge) = edges
            .get(&(o, bits))
            .and_then(|&eh| operands[o].edges.get(eh))
        else {
            continue;
        };
        let (Some(c), Some(h)) = (
            edge.curve.as_ref(),
            operands[o].half_edges.get(edge.half_edge),
        ) else {
            continue;
        };
        let twin = h
            .mate()
            .and_then(|m| operands[o].half_edges.get(m))
            .map(|t| t.origin);
        let ends = [c.t_start, c.t_end].map(|t| c.geometry.point_at(t));
        for v in std::iter::once(h.origin).chain(twin) {
            let Some(&g) = vertex_group.get(&(o, v.data().as_ffi())) else {
                continue;
            };
            let p = operands[o].vertices.get(v).map(|v| v.position);
            let Some(p) = p else {
                continue;
            };
            // Every end the vertex lies at, within its group's
            // tolerance: a closed edge (a rim circle) starts and ends at
            // its one vertex, and its last block ends there too; the
            // nearer end otherwise.
            let (d0, d1) = (p.distance_to(&ends[0]), p.distance_to(&ends[1]));
            let reach = res.groups[g].tolerance;
            if d0 <= d1 || d0 <= reach {
                list.push((c.t_start, g));
            }
            if d1 < d0 || d1 <= reach {
                list.push((c.t_end, g));
            }
        }
    }
    for list in out.values_mut() {
        list.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        list.dedup();
    }
    out
}

/// The first parameter of `(t0, t1)`, walking from `t0`, where `curve`
/// leaves the ball `(centre, radius)`; `None` when it never does.
fn exit(curve: &dyn CurveEval, (t0, t1): (f64, f64), centre: Point3, radius: f64) -> Option<f64> {
    let outside = |t: f64| curve.point_at(t).distance_to(&centre) > radius;
    let k = (1..=SAMPLES).find(|&k| outside(t0 + (t1 - t0) * k as f64 / SAMPLES as f64))?;
    let (mut lo, mut hi) = (
        t0 + (t1 - t0) * (k - 1) as f64 / SAMPLES as f64,
        t0 + (t1 - t0) * k as f64 / SAMPLES as f64,
    );
    for _ in 0..REFINE {
        let mid = 0.5 * (lo + hi);
        if outside(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some(hi)
}

/// The ball of a block end: its vertex's anchor and tolerance, or the
/// curve's own point with `linear` where no candidate resolved it.
fn ball(
    cands: &[CandidateVertex],
    res: &Resolution,
    end: Option<usize>,
    at: Point3,
    linear: f64,
) -> (Point3, f64) {
    end.map_or((at, linear), |g| {
        (
            cands[res.groups[g].anchor].position,
            res.groups[g].tolerance,
        )
    })
}

/// Every pave block of the records and the operand edges the
/// resolution cuts.
pub(crate) fn pave_blocks(
    cands: &[CandidateVertex],
    res: &Resolution,
    operands: [&BRep; 2],
    records: &[IntersectionRecord],
    linear: f64,
) -> Vec<PaveBlock> {
    let edges: BTreeMap<(usize, u64), EdgeHandle> = (0..2)
        .flat_map(|o| {
            operands[o]
                .edges
                .keys()
                .map(move |eh| ((o, eh.data().as_ffi()), eh))
        })
        .collect();
    let mut out = Vec::new();
    for (owner, list) in paves(cands, res, operands, &edges) {
        let Some((curve, (lo, hi))) = owner_curve(owner, operands, records, &edges) else {
            continue;
        };
        // An operand edge is a block where a candidate cuts it, or
        // whole where both its ends resolved (`paves`).
        let interior = list.iter().any(|&(t, _)| t > lo && t < hi);
        let whole = list.first().is_some_and(|&(t, _)| t <= lo)
            && list.last().is_some_and(|&(t, _)| t >= hi);
        if matches!(owner, BlockOwner::Edge(..)) && !interior && !whole {
            continue;
        }
        let mut stops: Vec<(f64, Option<usize>)> = Vec::new();
        if list.first().is_none_or(|&(t, _)| t > lo) {
            stops.push((lo, None));
        }
        stops.extend(list.iter().map(|&(t, g)| (t, Some(g))));
        if list.last().is_none_or(|&(t, _)| t < hi) {
            stops.push((hi, None));
        }
        for w in stops.windows(2) {
            let ((ta, ga), (tb, gb)) = (w[0], w[1]);
            // A piece from a vertex back to itself is no block: a
            // record whose ends met (`consequences`) or one vertex
            // reached twice.
            if tb <= ta || (ga.is_some() && ga == gb) {
                continue;
            }
            let (pa, ra) = ball(cands, res, ga, curve.point_at(ta), linear);
            let (pb, rb) = ball(cands, res, gb, curve.point_at(tb), linear);
            let shrunk = exit(curve, (ta, tb), pa, ra)
                .zip(exit(curve, (tb, ta), pb, rb))
                .filter(|(a, b)| a < b);
            out.push(PaveBlock {
                owner,
                range: (ta, tb),
                ends: [ga, gb],
                shrunk,
            });
        }
    }
    out
}

/// Pairs of blocks of different owners joining the same two distinct
/// vertices and lying within those vertices' tolerances of each other
/// along their shrunk ranges: common blocks (see the module docs).
pub(crate) fn common_block_pairs(
    blocks: &[PaveBlock],
    res: &Resolution,
    operands: [&BRep; 2],
    records: &[IntersectionRecord],
) -> Vec<(usize, usize)> {
    let edges: BTreeMap<(usize, u64), EdgeHandle> = (0..2)
        .flat_map(|o| {
            operands[o]
                .edges
                .keys()
                .map(move |eh| ((o, eh.data().as_ffi()), eh))
        })
        .collect();
    let mut by_ends: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    for (i, b) in blocks.iter().enumerate() {
        if let ([Some(a), Some(c)], Some(_)) = (b.ends, b.shrunk)
            && a != c
        {
            by_ends.entry((a.min(c), a.max(c))).or_default().push(i);
        }
    }
    let mut out = Vec::new();
    for ((ga, gb), list) in by_ends {
        let reach = res.groups[ga].tolerance.max(res.groups[gb].tolerance);
        for (k, &i) in list.iter().enumerate() {
            for &j in &list[k + 1..] {
                let (p, q) = (&blocks[i], &blocks[j]);
                if p.owner == q.owner {
                    continue;
                }
                let (Some((cp, _)), Some((cq, _)), Some((p0, p1)), Some((q0, q1))) = (
                    owner_curve(p.owner, operands, records, &edges),
                    owner_curve(q.owner, operands, records, &edges),
                    p.shrunk,
                    q.shrunk,
                ) else {
                    continue;
                };
                let within = (0..=SAMPLES).all(|s| {
                    let x = cp.point_at(p0 + (p1 - p0) * s as f64 / SAMPLES as f64);
                    let t = closest_point_in_window(cq, &x, q0.min(q1), q0.max(q1));
                    cq.point_at(t).distance_to(&x) <= reach
                });
                if within {
                    out.push((i, j));
                }
            }
        }
    }
    out
}
