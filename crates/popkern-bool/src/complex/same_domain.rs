//! Same-domain groups: the fragments that cover one coincident patch.
//!
//! A coincident patch is walked once on each operand's face. When both
//! walks cut it the same way, the two fragments bound the same edge set
//! and group directly. One operand can also cut it where the other does
//! not: its own seam (a cylinder wall's seam line crossing a coaxial
//! flush wall that has none) splits its copy into pieces joined by
//! edges no fragment of the partner uses. Those pieces group with the
//! partner's single fragment when their joined boundary equals it, and
//! that fragment leads: its loops bound the whole patch, so the joining
//! edges never enter the complex.

use std::collections::{BTreeMap, BTreeSet};

use super::lift::operand_index;
use super::spans::SpanRep;
use crate::merge::graph::{FacePair, FaceRef, IntersectionGraph};
use crate::merge::walk::Fragment;

/// The edge keys a fragment's loops use, degenerate keys removed:
/// the same-domain signature.
fn signature(
    graph: &IntersectionGraph,
    frag: &Fragment,
    reps: &BTreeMap<usize, SpanRep>,
) -> Vec<usize> {
    let mut keys: Vec<usize> = frag
        .outer_loop
        .iter()
        .chain(frag.inner_loops.iter().flatten())
        .filter_map(|&h| reps.get(&h.min(graph.half_edges[h].twin)))
        .map(|rep| rep.key)
        .collect();
    keys.sort_unstable();
    keys
}

/// Whether faces `x` (of `x_operand`) and `y` are a coincident pair.
fn coincident(overlaps: &BTreeSet<FacePair>, x: FaceRef, y: FaceRef) -> bool {
    overlaps.contains(&(x, y)) || overlaps.contains(&(y, x))
}

/// Group `fragments` into same-domain faces. Each group's first
/// fragment is its primary (its loops become the complex face's).
/// Groups follow their lowest fragment index; within a directly
/// matched group, operand A's fragment leads (the cellular rules keep
/// A's copy of a shared skin).
pub(crate) fn same_domain_groups(
    graph: &IntersectionGraph,
    fragments: &[Fragment],
    reps: &BTreeMap<usize, SpanRep>,
    overlaps: &[FacePair],
) -> Vec<Vec<usize>> {
    let overlap_set: BTreeSet<FacePair> = overlaps.iter().copied().collect();
    let sigs: Vec<Vec<usize>> = fragments
        .iter()
        .map(|f| signature(graph, f, reps))
        .collect();
    let mut by_signature: BTreeMap<&[usize], Vec<usize>> = BTreeMap::new();
    for (i, sig) in sigs.iter().enumerate() {
        by_signature.entry(sig.as_slice()).or_default().push(i);
    }
    let mut group_of = vec![usize::MAX; fragments.len()];
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for i in 0..fragments.len() {
        if group_of[i] != usize::MAX {
            continue;
        }
        let gi = groups.len();
        let mut group = vec![i];
        group_of[i] = gi;
        if !sigs[i].is_empty() {
            for &j in &by_signature[sigs[i].as_slice()] {
                if j == i || group_of[j] != usize::MAX {
                    continue;
                }
                let (fa, fb) = (fragments[i].face, fragments[j].face);
                if fa.operand != fb.operand && coincident(&overlap_set, fa, fb) {
                    group.push(j);
                    group_of[j] = gi;
                }
            }
        }
        group.sort_by_key(|&k| (operand_index(fragments[k].face.operand), k));
        groups.push(group);
    }
    join_split_copies(fragments, &sigs, &overlap_set, groups)
}

/// Merge each unmatched fragment with the partner face's unmatched
/// pieces whose joined boundary equals its own.
fn join_split_copies(
    fragments: &[Fragment],
    sigs: &[Vec<usize>],
    overlaps: &BTreeSet<FacePair>,
    mut groups: Vec<Vec<usize>>,
) -> Vec<Vec<usize>> {
    if overlaps.is_empty() {
        return groups;
    }
    let single = |g: &Vec<usize>| g.len() == 1;
    let unmatched: Vec<usize> = groups.iter().filter(|g| single(g)).map(|g| g[0]).collect();
    let mut absorbed: BTreeSet<usize> = BTreeSet::new();
    let mut joined: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for &j in &unmatched {
        if absorbed.contains(&j) || sigs[j].is_empty() {
            continue;
        }
        let fj = fragments[j].face;
        // The partner faces' unmatched pieces, by face.
        let mut pieces: BTreeMap<FaceRef, Vec<usize>> = BTreeMap::new();
        for &k in &unmatched {
            let fk = fragments[k].face;
            if k != j
                && !absorbed.contains(&k)
                && !joined.contains_key(&k)
                && fk.operand != fj.operand
                && coincident(overlaps, fj, fk)
            {
                pieces.entry(fk).or_default().push(k);
            }
        }
        let target: BTreeSet<usize> = sigs[j].iter().copied().collect();
        if target.len() != sigs[j].len() {
            continue;
        }
        for piece_set in pieces.values() {
            if let Some(component) = joined_component(piece_set, sigs, &target) {
                absorbed.extend(component.iter().copied());
                joined.insert(j, component);
                break;
            }
        }
    }
    if joined.is_empty() {
        return groups;
    }
    groups.retain(|g| !(single(g) && absorbed.contains(&g[0])));
    for g in &mut groups {
        if single(g)
            && let Some(component) = joined.get(&g[0])
        {
            g.extend(component.iter().copied());
        }
    }
    groups
}

/// A set of two or more of `pieces`, connected through edges both
/// bound, whose remaining edges are exactly `target`. Edges bounding
/// two pieces are interior to their union; a piece is a candidate only
/// when every edge it shares with no other piece lies on `target`.
fn joined_component(
    pieces: &[usize],
    sigs: &[Vec<usize>],
    target: &BTreeSet<usize>,
) -> Option<Vec<usize>> {
    let mut uses: BTreeMap<usize, usize> = BTreeMap::new();
    for &p in pieces {
        for &k in &sigs[p] {
            *uses.entry(k).or_default() += 1;
        }
    }
    // Pieces whose every edge is a target edge or an interior one.
    let fits = |p: usize| {
        sigs[p]
            .iter()
            .all(|k| target.contains(k) || (uses[k] == 2 && !target.contains(k)))
    };
    let candidates: Vec<usize> = pieces.iter().copied().filter(|&p| fits(p)).collect();
    // Connected components over interior edges.
    let mut seen: BTreeSet<usize> = BTreeSet::new();
    for &start in &candidates {
        if !seen.insert(start) {
            continue;
        }
        let mut component = vec![start];
        let mut cursor = 0;
        while cursor < component.len() {
            let p = component[cursor];
            cursor += 1;
            for &q in &candidates {
                if !seen.contains(&q)
                    && sigs[p]
                        .iter()
                        .any(|k| !target.contains(k) && sigs[q].contains(k))
                {
                    seen.insert(q);
                    component.push(q);
                }
            }
        }
        if component.len() < 2 {
            continue;
        }
        let mut boundary: BTreeMap<usize, usize> = BTreeMap::new();
        for &p in &component {
            for &k in &sigs[p] {
                *boundary.entry(k).or_default() += 1;
            }
        }
        let outer: BTreeSet<usize> = boundary
            .iter()
            .filter(|&(_, &n)| n == 1)
            .map(|(&k, _)| k)
            .collect();
        let interior_ok = boundary
            .iter()
            .all(|(k, &n)| n == 1 || (n == 2 && !target.contains(k)));
        if interior_ok && outer == *target {
            component.sort_unstable();
            return Some(component);
        }
    }
    None
}

#[cfg(test)]
#[path = "same_domain_tests.rs"]
mod tests;
