//! Mate-graph decomposition (design `six-axis-arm.md` §5.3): the order
//! in which root instances are solved, and which solve together.
//!
//! One Gauss–Newton over every free root is `6·N` unknowns with a dense
//! SVD per iteration — a 32-root arm is 192 columns, too slow for a
//! slider. Mechanisms are mostly trees: each link is fixed by its mates
//! to the link before it. So, from the ground out, a root whose mates
//! to roots already placed fix it alone (rank 6 in its own six
//! columns) — or which has no mate to any root still pending — is
//! solved alone with those mates, and joins the placed set. When no
//! root can be peeled (a loop: a four-bar, two fingers on one rail),
//! the pending roots connected to the first one through mates among
//! pending roots are solved together as one cluster.
//!
//! Solving a fully fixed root alone gives the same placement the joint
//! solve would; a root with no pending coupling has a block-diagonal
//! Jacobian, so its minimum-norm step is the joint solve's. Every mate
//! is re-checked at the final placements by the caller.

use crate::error::AssemblyError;
use crate::mate_linalg::rank_of;
use crate::mate_solve::Problem;
use crate::mates::{Mate, Resolved};
use crate::placement::Placement;

const RANK_REL: f64 = 1e-9;

/// The clusters (root indices), in solve order.
///
/// # Errors
///
/// A mate whose entities cannot take it (from the rank probe).
pub(crate) fn plan(
    resolved: &[(Mate, Resolved, Resolved, Option<Resolved>)],
    touches: &[Vec<usize>],
    free: &[usize],
    current: &[Placement],
    size: f64,
    roots: usize,
) -> Result<Vec<Vec<usize>>, AssemblyError> {
    let mut known: Vec<bool> = vec![true; roots];
    for &f in free {
        known[f] = false;
    }
    let mut pending: Vec<usize> = free.to_vec();
    let mut out = Vec::new();
    while !pending.is_empty() {
        let mut chosen = None;
        for (pi, &r) in pending.iter().enumerate() {
            let own: Vec<usize> = (0..resolved.len())
                .filter(|&i| {
                    touches[i].contains(&r) && touches[i].iter().all(|t| *t == r || known[*t])
                })
                .collect();
            if own.is_empty() {
                continue;
            }
            let coupled = (0..resolved.len()).any(|i| {
                touches[i].contains(&r) && touches[i].iter().any(|t| *t != r && !known[*t])
            });
            if !coupled || fixed_alone(resolved, &own, r, current, size)? {
                chosen = Some(pi);
                break;
            }
        }
        let cluster = match chosen {
            Some(pi) => vec![pending.remove(pi)],
            None => {
                // The pending roots linked to the first through mates
                // among pending roots, in pending order.
                let mut members = vec![pending[0]];
                let mut grew = true;
                while grew {
                    grew = false;
                    for t in touches {
                        if t.iter().any(|r| members.contains(r)) {
                            for &r in t {
                                if !known[r] && !members.contains(&r) {
                                    members.push(r);
                                    grew = true;
                                }
                            }
                        }
                    }
                }
                let ordered: Vec<usize> = pending
                    .iter()
                    .copied()
                    .filter(|r| members.contains(r))
                    .collect();
                pending.retain(|r| !members.contains(r));
                ordered
            }
        };
        for &r in &cluster {
            known[r] = true;
        }
        out.push(cluster);
    }
    Ok(out)
}

/// Whether mates `own` fix root `r` completely (rank 6 in its columns)
/// at the current placements.
fn fixed_alone(
    resolved: &[(Mate, Resolved, Resolved, Option<Resolved>)],
    own: &[usize],
    r: usize,
    current: &[Placement],
    size: f64,
) -> Result<bool, AssemblyError> {
    let sub: Vec<_> = own.iter().map(|&i| resolved[i].clone()).collect();
    let free = [r];
    let problem = Problem {
        resolved: &sub,
        free: &free,
        base: current,
        size,
    };
    let jacobian = problem.jacobian(&[0.0; 6])?;
    Ok(rank_of(&jacobian, 6, RANK_REL) == 6)
}
