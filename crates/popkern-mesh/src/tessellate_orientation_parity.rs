//! Face-pair parity votes resolved into one orientation per face: a
//! maximum-weight spanning forest over the pair constraints (split from
//! `tessellate_orientation.rs` under the file-size rule).
//!
//! Each face pair sharing an edge votes whether their meshes run the
//! edge the same way (one of them inside out) or opposite ways. The
//! constraints are taken strongest first, by the margin of the pair's
//! majority, and one contradicting a cycle the stronger ones already
//! closed is dropped: the weaker reading is the one a thin strip beside
//! the edge contaminated. Taken in index order instead, the first
//! assignment won whatever its weight (measured: two cones touching
//! beside a rim, the union's lens edges 1e-6 apart along 1.2e-4; the
//! cones' shared 3e-4 edge read `same` at 1.1e-4 against the cycle's
//! `opposite` at 1.3e-4 and 4e-3, three faces meshed inside out and the
//! volume read 7.18e-9 for 6.13e-9).

use super::PairWindingVotes;

/// Root of `f` and `f`'s parity relative to it (true: the two face
/// meshes' windings differ).
fn find(parent: &mut [usize], parity: &mut [bool], f: usize) -> (usize, bool) {
    let mut path = Vec::new();
    let mut x = f;
    while parent[x] != x {
        path.push(x);
        x = parent[x];
    }
    // Compress: walk the path from the root's side, each node's parity
    // to the root the sum of its own and its parent's.
    let mut acc = false;
    for &y in path.iter().rev() {
        acc ^= parity[y];
        parity[y] = acc;
        parent[y] = x;
    }
    (x, if path.is_empty() { false } else { parity[f] })
}

/// Every face's flip (its mesh to be reversed relative to its
/// component's first face) and the components, each listed in face
/// order and the components in order of their first face. Faces of no
/// constraint are components of their own.
pub(super) fn resolve_parity(
    n: usize,
    votes: &PairWindingVotes,
) -> (Vec<Option<bool>>, Vec<Vec<usize>>) {
    let trace = popkern_core::trace::trace_enabled("mesh-orientation");
    let mut constraints: Vec<(usize, usize, bool, f64)> = Vec::new();
    for (&(f1, f2), &(same, opposite)) in votes {
        if same == opposite {
            // Tied evidence: no trustworthy constraint.
            if trace {
                eprintln!("  [orient] pair ({f1},{f2}) tied votes {same}:{opposite} — unlinked");
            }
            continue;
        }
        let same_dir = same > opposite;
        if same.min(opposite) > 0.0 && trace {
            eprintln!(
                "  [orient] pair ({f1},{f2}) conflicted votes same={same:.3} opposite={opposite:.3} — majority {}",
                if same_dir { "same" } else { "opposite" }
            );
        }
        constraints.push((f1, f2, same_dir, (same - opposite).abs()));
    }
    // Strongest first; the pair key breaks ties so the order never rests
    // on the vote map's iteration.
    constraints.sort_by(|a, b| b.3.total_cmp(&a.3).then((a.0, a.1).cmp(&(b.0, b.1))));
    let mut parent: Vec<usize> = (0..n).collect();
    let mut parity = vec![false; n];
    for (f1, f2, same_dir, margin) in constraints {
        let (r1, p1) = find(&mut parent, &mut parity, f1);
        let (r2, p2) = find(&mut parent, &mut parity, f2);
        if r1 == r2 {
            if (p1 ^ p2) != same_dir && trace {
                eprintln!(
                    "  [orient] pair ({f1},{f2}) margin {margin:.5} contradicts a stronger cycle — dropped"
                );
            }
            continue;
        }
        // Root the merged tree at the lower face so each component's
        // first face keeps its winding.
        let (keep, moved, rel) = if r1 < r2 {
            (r1, r2, p1 ^ p2 ^ same_dir)
        } else {
            (r2, r1, p1 ^ p2 ^ same_dir)
        };
        parent[moved] = keep;
        parity[moved] = rel;
    }
    let mut flip = Vec::with_capacity(n);
    let mut components: Vec<Vec<usize>> = Vec::new();
    let mut slot: Vec<Option<usize>> = vec![None; n];
    for f in 0..n {
        let (root, p) = find(&mut parent, &mut parity, f);
        flip.push(Some(p));
        let k = *slot[root].get_or_insert_with(|| {
            components.push(Vec::new());
            components.len() - 1
        });
        components[k].push(f);
    }
    (flip, components)
}

#[cfg(test)]
#[path = "tessellate_orientation_parity_tests.rs"]
mod tests;
