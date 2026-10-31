//! Candidate resolution (plan §4.3a): how many real vertices each group
//! of nearby candidates becomes, with the reason recorded.
//!
//! **Components.** Candidates are linked when their tolerance balls
//! overlap (distance below the sum of the radii); a sweep over x keeps
//! that O(n log n) plus the overlapping pairs.
//!
//! **Compatibility.** Two candidates may not be one vertex when that
//! would identify points the operands keep apart:
//! - two distinct vertices of one operand (a valid operand's vertices
//!   are distinct by construction);
//! - two points on one edge, or on one record, separated along it by
//!   more than the sum of their zones (arc length, not chord: an edge
//!   curling back near itself keeps both points; the zone is the radius
//!   elongated along a shallow crossing, `incidence::add_crossing_zones`);
//! - two points on two records of one face pair, one in its record's
//!   interior: two branches passing close, where the SSI saw no branch
//!   point (`branches`).
//!
//! Incompatibility is a hard constraint; proximity never overrides it.
//!
//! **Partition.** Within a component, in canonical order (lowest
//! dimension first, so a group containing an operand vertex is anchored
//! at it), each candidate joins the first group it is compatible with
//! whose members' balls it overlaps and whose grown tolerance stays
//! within the cap; otherwise it starts a group. The vertex sits at its
//! anchor's position, never a midpoint. Its tolerance grows to cover
//! every member's ball (Jackson §3.3), recorded, up to the cap or the
//! two candidates' zones, whichever is larger (growth a shallow crossing
//! explains); a candidate that only that bound kept out of a group it
//! touches is `ToleranceGrowthExceeded`:
//! a chain of points each within tolerance of the next never silently
//! becomes one vertex.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::BRep;

use super::candidates::{CandidateVertex, Incidence};
use crate::intersect::IntersectionRecord;

/// Why two candidates of one component stay distinct vertices.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Incompatibility {
    /// Both are vertices of one operand.
    DistinctOperandVertices,
    /// Both lie on one operand edge, `arc` apart along it.
    SeparatedAlongEdge {
        /// Arc length between them.
        arc: f64,
    },
    /// Both lie on one record, `arc` apart along it.
    SeparatedAlongRecord {
        /// Arc length between them.
        arc: f64,
    },
    /// They lie on two records of one face pair, one of them in its
    /// record's interior: a branch the SSI marched through, seeing no
    /// branch point, passing within tolerance of another (`branches`).
    DistinctBranches,
}

/// One resolved vertex.
#[derive(Debug, Clone)]
pub(crate) struct VertexGroup {
    /// The candidate it is anchored at (its position is the vertex's).
    pub anchor: usize,
    /// Every member, the anchor first.
    pub members: Vec<usize>,
    /// The vertex tolerance: the anchor's radius grown to cover every
    /// member's ball.
    pub tolerance: f64,
}

/// A resolution that cannot commit its vertices.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum ResolutionDecline {
    /// A candidate touches a group only the growth cap keeps it out of.
    ToleranceGrowthExceeded {
        /// The group.
        group: usize,
        /// The candidate.
        candidate: usize,
        /// The tolerance the group would need.
        needed: f64,
        /// The cap.
        cap: f64,
    },
}

/// Resolution parameters.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ResolutionOptions {
    /// The largest vertex tolerance a group may grow to.
    pub growth_cap: f64,
}

/// The committed vertex table.
#[derive(Debug, Clone, Default)]
pub(crate) struct Resolution {
    /// The vertices.
    pub groups: Vec<VertexGroup>,
    /// Each candidate's group.
    pub group_of: Vec<usize>,
    /// Proximity components with more than one candidate.
    pub components: usize,
    /// Pairs of one component kept apart, with the rule that did it.
    pub separations: Vec<(usize, usize, Incompatibility)>,
}

/// The operands and records the candidates' incidences refer to.
pub(crate) struct ResolutionContext<'a> {
    /// Operands A and B.
    pub operands: [&'a BRep; 2],
    /// The records, indexed as the candidates' `Record` incidences.
    pub records: &'a [IntersectionRecord],
}

/// Arc length of `curve` between `t0` and `t1`, by 32 chords.
fn arc_between(point_at: impl Fn(f64) -> Point3, t0: f64, t1: f64) -> f64 {
    const CHORDS: usize = 32;
    let mut prev = point_at(t0);
    let mut sum = 0.0;
    for k in 1..=CHORDS {
        let p = point_at(t0 + (t1 - t0) * k as f64 / CHORDS as f64);
        sum += prev.distance_to(&p);
        prev = p;
    }
    sum
}

/// Arc length between two parameters of a curve over `(lo, hi)`; on a
/// closed curve the shorter way round.
fn separation(point_at: impl Fn(f64) -> Point3 + Copy, (lo, hi): (f64, f64), t: (f64, f64)) -> f64 {
    let (t0, t1) = (t.0.min(t.1), t.0.max(t.1));
    let inner = arc_between(point_at, t0, t1);
    let closed = point_at(lo).distance_to(&point_at(hi)) <= 1e-12 * (1.0 + inner);
    if closed {
        let outer = arc_between(point_at, lo, t0) + arc_between(point_at, t1, hi);
        inner.min(outer)
    } else {
        inner
    }
}

fn incompatibility(
    ctx: &ResolutionContext<'_>,
    p: &CandidateVertex,
    q: &CandidateVertex,
) -> Option<Incompatibility> {
    for ip in &p.incidence {
        for iq in &q.incidence {
            match (*ip, *iq) {
                (
                    Incidence::Vertex {
                        operand: oa,
                        vertex: va,
                    },
                    Incidence::Vertex {
                        operand: ob,
                        vertex: vb,
                    },
                ) if oa == ob && va != vb => {
                    return Some(Incompatibility::DistinctOperandVertices);
                }
                (
                    Incidence::Edge {
                        operand: oa,
                        edge: ea,
                        t: ta,
                    },
                    Incidence::Edge {
                        operand: ob,
                        edge: eb,
                        t: tb,
                    },
                ) if oa == ob && ea == eb => {
                    let Some(c) = ctx.operands[oa]
                        .edges
                        .get(ea)
                        .and_then(|e| e.curve.as_ref())
                    else {
                        continue;
                    };
                    let range = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
                    let arc = separation(|t| c.geometry.point_at(t), range, (ta, tb));
                    if arc > p.zone + q.zone {
                        return Some(Incompatibility::SeparatedAlongEdge { arc });
                    }
                }
                (
                    Incidence::Record {
                        record: ra, s: sa, ..
                    },
                    Incidence::Record {
                        record: rb, s: sb, ..
                    },
                ) if ra == rb => {
                    let Some(r) = ctx.records.get(ra) else {
                        continue;
                    };
                    let curve = &r.curve.curve_3d;
                    let arc = separation(|s| curve.point_at(s), curve.parameter_range(), (sa, sb));
                    if arc > p.zone + q.zone {
                        return Some(Incompatibility::SeparatedAlongRecord { arc });
                    }
                }
                (
                    Incidence::Record {
                        record: ra,
                        end: ea,
                        ..
                    },
                    Incidence::Record {
                        record: rb,
                        end: eb,
                        ..
                    },
                ) if (ea.is_none() || eb.is_none())
                    && super::branches::same_pair(ctx.records, ra, rb) =>
                {
                    return Some(Incompatibility::DistinctBranches);
                }
                _ => {}
            }
        }
    }
    None
}

/// Disjoint-set forest.
struct Forest(Vec<usize>);

impl Forest {
    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]];
            x = self.0[x];
        }
        x
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
            self.0[hi] = lo;
        }
    }
}

/// Whether two candidates' tolerance zones overlap: their balls, or,
/// when both lie on one operand edge, their reaches along it.
fn linked(p: &CandidateVertex, q: &CandidateVertex) -> bool {
    let d = p.position.distance_to(&q.position);
    if d < p.radius + q.radius {
        return true;
    }
    d < p.zone + q.zone
        && p.incidence.iter().any(|ip| {
            q.incidence.iter().any(|iq| {
                matches!(
                    (ip, iq),
                    (
                        Incidence::Edge { operand: oa, edge: ea, .. },
                        Incidence::Edge { operand: ob, edge: eb, .. },
                    ) if oa == ob && ea == eb
                )
            })
        })
}

/// Proximity components: candidates whose zones overlap, transitively.
/// Each component lists its members in candidate order.
fn components(cands: &[CandidateVertex]) -> Vec<Vec<usize>> {
    let n = cands.len();
    let mut by_x: Vec<usize> = (0..n).collect();
    by_x.sort_by(|&i, &j| cands[i].position.x.total_cmp(&cands[j].position.x));
    let r_max = cands.iter().map(|c| c.zone).fold(0.0_f64, f64::max);
    let mut forest = Forest((0..n).collect());
    for (k, &i) in by_x.iter().enumerate() {
        let ci = &cands[i];
        for &j in &by_x[k + 1..] {
            let cj = &cands[j];
            if cj.position.x - ci.position.x > ci.zone + r_max {
                break;
            }
            if linked(ci, cj) {
                forest.union(i, j);
            }
        }
    }
    let mut out: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
    for i in 0..n {
        let root = forest.find(i);
        out.entry(root).or_default().push(i);
    }
    out.into_values().collect()
}

/// Resolve `cands` (sorted by key, as `collect_candidates` returns
/// them) into vertices.
pub(crate) fn resolve(
    cands: &[CandidateVertex],
    ctx: &ResolutionContext<'_>,
    options: ResolutionOptions,
) -> Result<Resolution, ResolutionDecline> {
    let mut res = Resolution {
        group_of: vec![usize::MAX; cands.len()],
        ..Resolution::default()
    };
    for comp in components(cands) {
        if comp.len() > 1 {
            res.components += 1;
        }
        let first_group = res.groups.len();
        for &c in &comp {
            let cand = &cands[c];
            let mut joined = false;
            let mut capped: Option<(usize, f64)> = None;
            for g in first_group..res.groups.len() {
                let group = &res.groups[g];
                let touches = group.members.iter().any(|&m| linked(&cands[m], cand));
                if !touches {
                    continue;
                }
                let clash = group
                    .members
                    .iter()
                    .find_map(|&m| incompatibility(ctx, &cands[m], cand).map(|why| (m, why)));
                if let Some((m, why)) = clash {
                    res.separations.push((m, c, why));
                    continue;
                }
                let anchor = &cands[group.anchor];
                let needed = group
                    .tolerance
                    .max(anchor.position.distance_to(&cand.position) + cand.radius);
                // Growth the crossing angles explain is allowed past the
                // cap; the cap bounds what nothing explains.
                if needed > options.growth_cap.max(anchor.zone + cand.zone) {
                    capped.get_or_insert((g, needed));
                    continue;
                }
                let group = &mut res.groups[g];
                group.members.push(c);
                group.tolerance = needed;
                res.group_of[c] = g;
                joined = true;
                break;
            }
            if joined {
                continue;
            }
            if let Some((group, needed)) = capped {
                return Err(ResolutionDecline::ToleranceGrowthExceeded {
                    group,
                    candidate: c,
                    needed,
                    cap: options.growth_cap,
                });
            }
            res.group_of[c] = res.groups.len();
            res.groups.push(VertexGroup {
                anchor: c,
                members: vec![c],
                tolerance: cand.radius,
            });
        }
    }
    Ok(res)
}
