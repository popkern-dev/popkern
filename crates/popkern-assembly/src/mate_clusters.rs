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
//! root can be peeled, the pending roots that only *follow* others are
//! set aside first, from the outside in: a root whose mates fix it
//! alone and ask nothing of the roots they reach (the rank of their
//! rows in its own six columns is 6 and is the rank in all of theirs)
//! — a reducer's rotor, discs and cage on a released joint (PK-808),
//! a finger on its rail. Each is solved alone after the roots it
//! follows. What remains (a loop: a four-bar, a released chain held
//! by a target) is solved as one cluster of the pending roots
//! connected to the first one through mates among pending roots.
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
    // Roots set aside to follow, and whether a mate reaches one (its
    // mates are the follower's own, solved with it).
    let mut followers: Vec<usize> = Vec::new();
    let mut set_aside = vec![false; roots];
    while !pending.is_empty() {
        let live = |i: usize| !touches[i].iter().any(|t| set_aside[*t]);
        let mut chosen = None;
        for (pi, &r) in pending.iter().enumerate() {
            let own: Vec<usize> = (0..resolved.len())
                .filter(|&i| {
                    live(i)
                        && touches[i].contains(&r)
                        && touches[i].iter().all(|t| *t == r || known[*t])
                })
                .collect();
            if own.is_empty() {
                continue;
            }
            let coupled = (0..resolved.len()).any(|i| {
                live(i)
                    && touches[i].contains(&r)
                    && touches[i].iter().any(|t| *t != r && !known[*t])
            });
            if !coupled || fixed_alone(resolved, &own, r, current, size)? {
                chosen = Some(pi);
                break;
            }
        }
        if chosen.is_none() {
            let mut any = false;
            for pi in (0..pending.len()).rev() {
                let r = pending[pi];
                let own: Vec<usize> = (0..resolved.len())
                    .filter(|&i| live(i) && touches[i].contains(&r))
                    .collect();
                if follows(resolved, touches, &own, r, &known, current, size)? {
                    set_aside[r] = true;
                    followers.push(pending.remove(pi));
                    any = true;
                    break;
                }
            }
            if any {
                continue;
            }
        }
        let live = |i: usize| !touches[i].iter().any(|t| set_aside[*t]);
        let cluster = match chosen {
            Some(pi) => vec![pending.remove(pi)],
            None => {
                // The pending roots linked to the first through mates
                // among pending roots, in pending order.
                let mut members = vec![pending[0]];
                let mut grew = true;
                while grew {
                    grew = false;
                    for (i, t) in touches.iter().enumerate() {
                        if live(i) && t.iter().any(|r| members.contains(r)) {
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
    // Followers after what they follow: the last set aside first.
    out.extend(followers.into_iter().rev().map(|r| vec![r]));
    Ok(out)
}

/// Whether mates `own` fix root `r` alone and ask nothing of the
/// pending roots they reach: rank 6 in `r`'s columns, and no more in
/// all of theirs. Such a root moves to satisfy them wherever the others
/// land.
fn follows(
    resolved: &[(Mate, Resolved, Resolved, Option<Resolved>)],
    touches: &[Vec<usize>],
    own: &[usize],
    r: usize,
    known: &[bool],
    current: &[Placement],
    size: f64,
) -> Result<bool, AssemblyError> {
    if own.is_empty() {
        return Ok(false);
    }
    let mut free = vec![r];
    for &i in own {
        for &t in &touches[i] {
            if !known[t] && !free.contains(&t) {
                free.push(t);
            }
        }
    }
    let sub: Vec<_> = own.iter().map(|&i| resolved[i].clone()).collect();
    let problem = Problem {
        resolved: &sub,
        free: &free,
        base: current,
        size,
        limits: &[],
    };
    let jacobian = problem.jacobian(&vec![0.0; 6 * free.len()])?;
    let mine: Vec<Vec<f64>> = jacobian.iter().map(|row| row[..6].to_vec()).collect();
    let rank = rank_of(&mine, 6, RANK_REL);
    Ok(rank == 6 && rank_of(&jacobian, 6 * free.len(), RANK_REL) == 6)
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
        limits: &[],
    };
    let jacobian = problem.jacobian(&[0.0; 6])?;
    Ok(rank_of(&jacobian, 6, RANK_REL) == 6)
}
