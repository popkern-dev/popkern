//! The fewest vertices a component can be (plan §4.3a step 2).
//!
//! The canonical-order pass (`resolve`) puts each candidate in the
//! first group that takes it. That is minimal whenever it reaches the
//! lower bound, the size of the largest set of mutually incompatible
//! candidates: no partition puts two of them in one vertex. When it
//! does not, an early join can cost a vertex later: a point compatible
//! with two vertices joins the first, and a point only it could have
//! carried into the second starts a third.
//!
//! A component of at most `exact_limit` candidates is then searched
//! exhaustively for a partition with fewer vertices (or, when the pass
//! left an order conflict, one that keeps every record's order),
//! visiting candidates in canonical order, existing groups before a new
//! one, so the first partition found at the fewest vertices is the
//! canonical one. Every group of it obeys the pass's own joining rules:
//! touching, pairwise compatible, within the growth bound, in record
//! order. A larger component the pass cannot prove minimal declines as
//! `CandidateResolutionAmbiguous`: the table never guesses.
//!
//! Measured (bool, examples, mesh, blend; table on): 2791/2791 pass;
//! across 1449 shadowed booleans the pass was provably minimal in every
//! component (no re-partition, no ambiguity decline). The search is the
//! guarantee behind the pass, not a change to what it serves there.

use super::candidates::CandidateVertex;
use super::resolve::{
    Incompatibility, Resolution, ResolutionContext, ResolutionDecline, ResolutionOptions,
    VertexGroup,
};

/// Pairwise facts of one component, in its local order.
pub(crate) struct Pairwise {
    n: usize,
    linked: Vec<bool>,
    clash: Vec<Option<Incompatibility>>,
}

impl Pairwise {
    /// Whether local members `i` and `j` touch.
    pub(crate) fn linked(&self, i: usize, j: usize) -> bool {
        self.linked[i * self.n + j]
    }

    /// Why local members `i` and `j` may not be one vertex.
    pub(crate) fn clash(&self, i: usize, j: usize) -> Option<Incompatibility> {
        self.clash[i * self.n + j]
    }
}

/// The pairwise facts of component `comp`.
pub(crate) fn pairwise(
    cands: &[CandidateVertex],
    ctx: &ResolutionContext<'_>,
    comp: &[usize],
) -> Pairwise {
    let n = comp.len();
    let mut pw = Pairwise {
        n,
        linked: vec![false; n * n],
        clash: vec![None; n * n],
    };
    for i in 0..n {
        for j in i + 1..n {
            let (p, q) = (&cands[comp[i]], &cands[comp[j]]);
            let linked = super::resolve::linked(p, q);
            let clash = super::resolve::incompatibility(ctx, p, q);
            for (a, b) in [(i, j), (j, i)] {
                pw.linked[a * n + b] = linked;
                pw.clash[a * n + b] = clash;
            }
        }
    }
    pw
}

/// The size of the largest set of mutually incompatible members: no
/// partition has fewer vertices. Exact (Bron–Kerbosch with pivoting
/// over bit sets) up to 64 members; a greedy clique, still a valid
/// bound, above that.
pub(crate) fn lower_bound(pw: &Pairwise) -> usize {
    if pw.n == 0 {
        return 0;
    }
    if pw.n > 64 {
        let mut clique: Vec<usize> = Vec::new();
        for i in 0..pw.n {
            if clique.iter().all(|&j| pw.clash(i, j).is_some()) {
                clique.push(i);
            }
        }
        return clique.len();
    }
    let adj: Vec<u64> = (0..pw.n)
        .map(|i| {
            (0..pw.n)
                .filter(|&j| j != i && pw.clash(i, j).is_some())
                .fold(0u64, |m, j| m | 1 << j)
        })
        .collect();
    let all = if pw.n == 64 {
        u64::MAX
    } else {
        (1u64 << pw.n) - 1
    };
    let mut best = 0;
    bron_kerbosch(&adj, 0, all, 0, &mut best);
    best
}

fn bron_kerbosch(adj: &[u64], size: usize, mut p: u64, mut x: u64, best: &mut usize) {
    if p == 0 {
        *best = (*best).max(size);
        return;
    }
    if size + p.count_ones() as usize <= *best {
        return;
    }
    let pivot = (p | x).trailing_zeros() as usize;
    let mut todo = p & !adj[pivot];
    while todo != 0 {
        let v = todo.trailing_zeros() as usize;
        todo &= todo - 1;
        bron_kerbosch(adj, size + 1, p & adj[v], x & adj[v], best);
        p &= !(1 << v);
        x |= 1 << v;
    }
}

struct Search<'a> {
    cands: &'a [CandidateVertex],
    ctx: &'a ResolutionContext<'a>,
    comp: &'a [usize],
    pw: &'a Pairwise,
    options: ResolutionOptions,
    floor: usize,
    best: Option<Vec<VertexGroup>>,
    below: usize,
}

impl Search<'_> {
    fn visit(&mut self, k: usize, groups: &mut Vec<VertexGroup>, local: &mut Vec<Vec<usize>>) {
        if groups.len() >= self.below || self.below <= self.floor {
            return;
        }
        if k == self.comp.len() {
            self.below = groups.len();
            self.best = Some(groups.clone());
            return;
        }
        let c = self.comp[k];
        let cand = &self.cands[c];
        for g in 0..groups.len() {
            let fits = local[g].iter().any(|&m| self.pw.linked(m, k))
                && local[g].iter().all(|&m| self.pw.clash(m, k).is_none());
            if !fits {
                continue;
            }
            let anchor = &self.cands[groups[g].anchor];
            let needed = groups[g]
                .tolerance
                .max(anchor.position.distance_to(&cand.position) + cand.radius);
            if needed > self.options.growth_cap.max(anchor.zone + cand.zone)
                || super::record_order::reversal(self.cands, self.ctx.records, groups, Some(g), c)
                    .is_some()
            {
                continue;
            }
            let before = groups[g].tolerance;
            groups[g].members.push(c);
            groups[g].tolerance = needed;
            local[g].push(k);
            self.visit(k + 1, groups, local);
            local[g].pop();
            groups[g].members.pop();
            groups[g].tolerance = before;
        }
        if super::record_order::reversal(self.cands, self.ctx.records, groups, None, c).is_none() {
            groups.push(VertexGroup {
                anchor: c,
                members: vec![c],
                tolerance: cand.radius,
            });
            local.push(vec![k]);
            self.visit(k + 1, groups, local);
            local.pop();
            groups.pop();
        }
    }
}

/// The canonical partition of component `comp` into fewer than `below`
/// vertices, every record's order kept; `None` when there is none.
/// `floor` is a lower bound on any partition: the search stops once it
/// reaches it.
pub(crate) fn exact(
    cands: &[CandidateVertex],
    ctx: &ResolutionContext<'_>,
    comp: &[usize],
    pw: &Pairwise,
    options: ResolutionOptions,
    (floor, below): (usize, usize),
) -> Option<Vec<VertexGroup>> {
    let mut search = Search {
        cands,
        ctx,
        comp,
        pw,
        options,
        floor,
        best: None,
        below,
    };
    search.visit(0, &mut Vec::new(), &mut Vec::new());
    search.best
}

/// Settle component `comp`, the `index`th, whose canonical-order groups
/// start at `starts.0` (its separations at `starts.1`, its order
/// conflicts at `starts.2`): kept when provably minimal and in order,
/// replaced by the exact partition when that has fewer vertices (or
/// keeps an order the pass could not), declined when the component is
/// too large to search and the pass cannot be proved minimal.
pub(crate) fn settle(
    cands: &[CandidateVertex],
    ctx: &ResolutionContext<'_>,
    options: ResolutionOptions,
    (index, comp): (usize, &[usize]),
    starts: (usize, usize, usize),
    res: &mut Resolution,
) -> Result<(), ResolutionDecline> {
    let (first_group, first_separation, first_conflict) = starts;
    let vertices = res.groups.len() - first_group;
    let conflicted = res.order_conflicts.len() > first_conflict;
    if vertices <= 1 && !conflicted {
        return Ok(());
    }
    let pw = pairwise(cands, ctx, comp);
    let floor = lower_bound(&pw);
    if vertices <= floor && !conflicted {
        return Ok(());
    }
    if comp.len() > options.exact_limit {
        return if vertices > floor {
            Err(ResolutionDecline::CandidateResolutionAmbiguous {
                component: index,
                members: comp.to_vec(),
                vertices,
                lower_bound: floor,
            })
        } else {
            Ok(())
        };
    }
    let below = if conflicted { vertices + 1 } else { vertices };
    let Some(groups) = exact(cands, ctx, comp, &pw, options, (floor, below)) else {
        return Ok(());
    };
    res.groups.truncate(first_group);
    for g in groups {
        for &m in &g.members {
            res.group_of[m] = res.groups.len();
        }
        res.groups.push(g);
    }
    res.separations.truncate(first_separation);
    for i in 0..comp.len() {
        for j in i + 1..comp.len() {
            if res.group_of[comp[i]] != res.group_of[comp[j]]
                && pw.linked(i, j)
                && let Some(why) = pw.clash(i, j)
            {
                res.separations.push((comp[i], comp[j], why));
            }
        }
    }
    res.order_conflicts.truncate(first_conflict);
    res.exact_partitions += 1;
    Ok(())
}
