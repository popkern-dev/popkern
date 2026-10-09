//! Exact motion prediction before the mate solve (six-axis-arm §5.3).
//!
//! A slider moves one parameter; the solve then has to find every root
//! that moves with it. Gauss–Newton from the previous placements sees a
//! revolute's whole turn at once — a reducer's rotor geared 29:1 turns
//! 145° for a 5° slider step — and its angle residual is an `atan2`
//! that wraps, so a large turn needed continuation: the driven angles
//! stepped by at most 30° with a solve per step (five solves for one
//! slider step, hundreds for a full revolution).
//!
//! But a revolute's motion is known exactly. In the solve's own order
//! (from the ground out, see `mate_clusters`), each root solved alone
//! is first carried by the motion of the roots that fix it, then — when
//! it is a joint's follower whose axis is still coincident — turned
//! about the joint axis by the angle its mate asks for. Roots solved
//! together (a loop) are carried rigidly when everything fixing them
//! moved alike, which keeps their assembly mode. The solve then starts
//! at the answer and only polishes. A loop whose own driven mate jumps
//! is the caller's: it keeps continuation for those.

use crate::joints::Joint;
use crate::mates::{Mate, MateKind, Primitive, Resolved};
use crate::placement::Placement;

/// Carried placements of equal motions agree to this (rotation entries;
/// translations relative to the assembly's size).
const SAME_MOTION: f64 = 1e-12;

/// A follower whose axis is off its base axis by more than this (sine,
/// or distance relative to the size) is left to the solve.
const ON_AXIS: f64 = 1e-9;

/// Move `current` (the previous solution) toward the solution of the
/// mates in `resolved`, cluster by cluster in `clusters` order; roots in
/// no cluster are held. Returns, per cluster, whether its motion was
/// predicted (see `ClusterReport::predicted`).
pub(crate) fn predict(
    joints: &[Joint],
    resolved: &[(Mate, Resolved, Resolved, Option<Resolved>)],
    touches: &[Vec<usize>],
    clusters: &[Vec<usize>],
    current: &mut [Placement],
    size: f64,
) -> Vec<bool> {
    let identity = Placement::identity();
    // Each root's motion this pass; held roots do not move.
    let mut motion: Vec<Option<Placement>> = vec![Some(identity); current.len()];
    for cluster in clusters {
        for &r in cluster {
            motion[r] = None;
        }
    }
    let same = |a: &Placement, b: &Placement| {
        let (ma, mb) = (&a.transform().matrix, &b.transform().matrix);
        (0..3).all(|i| {
            (0..4).all(|j| {
                let s = if j == 3 { size } else { 1.0 };
                (ma[i][j] - mb[i][j]).abs() <= SAME_MOTION * s
            })
        })
    };
    let mut predicted = Vec::with_capacity(clusters.len());
    for cluster in clusters {
        // The mates fixing the cluster from outside, and the motions of
        // the roots they reach that are already placed.
        let outside: Vec<Placement> = (0..resolved.len())
            .filter(|&i| touches[i].iter().any(|t| cluster.contains(t)))
            .flat_map(|i| touches[i].iter().copied())
            .filter(|t| !cluster.contains(t))
            .filter_map(|t| motion[t])
            .collect();
        let common = match outside.split_first() {
            Some((first, rest)) if rest.iter().all(|m| same(m, first)) => Some(*first),
            _ => None,
        };
        // A single root that is a joint's follower rides its joint's
        // base, then turns to the joint's angle.
        let joint = match cluster.as_slice() {
            [r] => joints.iter().find(|j| {
                j.mates
                    .iter()
                    .all(|&m| touches.get(m).is_some_and(|t| t.contains(r)))
                    && resolved.get(j.mates[0]).is_some_and(|m| m.2.root == *r)
            }),
            _ => None,
        };
        let carry = match joint {
            Some(j) => resolved.get(j.mates[0]).and_then(|m| motion[m.1.root]),
            // Fixed by nothing already placed: it starts where it was.
            None if outside.is_empty() => Some(identity),
            None => common,
        };
        predicted.push(carry.is_some());
        for &r in cluster {
            let before = current[r];
            if let Some(d) = carry {
                current[r] = before.then(&d);
            }
            if let Some(j) = joint
                && let Some(turned) = turn_to_angle(j, resolved, current, size)
            {
                current[r] = turned;
            }
            motion[r] = Some(before.inverse().then(&current[r]));
        }
    }
    predicted
}

/// `current[follower]` turned about joint `j`'s axis so its angle mate
/// holds; `None` when the follower is off the axis (the solve places it).
fn turn_to_angle(
    j: &Joint,
    resolved: &[(Mate, Resolved, Resolved, Option<Resolved>)],
    current: &[Placement],
    size: f64,
) -> Option<Placement> {
    let world = |r: &Resolved| r.local.placed(&current[r.root]);
    let (_, axis_a, axis_b, _) = resolved.get(j.mates[0])?;
    let (Primitive::Axis(pa, ka), Primitive::Axis(pb, kb)) = (world(axis_a), world(axis_b)) else {
        return None;
    };
    // Coincident axes only: turning about the base axis is then the
    // joint's own motion.
    if ka.cross(&kb).length() > ON_AXIS || (pb - pa).cross(&ka).length() > ON_AXIS * size {
        return None;
    }
    let (mate, zero_a, zero_b, about) = resolved.get(j.mates[2])?;
    let MateKind::Angle { degrees, .. } = mate.kind else {
        return None;
    };
    let k = about.as_ref().map(world)?.direction()?;
    let (ua, ub) = (world(zero_a).direction()?, world(zero_b).direction()?);
    use popkern_core::geometry::PortableMath;
    let signed = ua.cross(&ub).dot(&k).patan2(ua.dot(&ub)).to_degrees();
    let delta = (degrees - signed + 180.0).rem_euclid(360.0) - 180.0;
    let turn = Placement::rotation(pa, k, delta).ok()?;
    Some(current[axis_b.root].then(&turn))
}
