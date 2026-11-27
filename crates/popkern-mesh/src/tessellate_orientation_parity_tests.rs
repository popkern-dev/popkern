use super::resolve_parity;
use crate::tessellate::PairWindingVotes;

/// A face pair and its (same, opposite) vote weights.
type PairVote = ((usize, usize), (f64, f64));

fn votes(pairs: &[PairVote]) -> PairWindingVotes {
    pairs.iter().copied().collect()
}

/// The union of two cones touching beside a rim, as measured: each
/// cone and its disc vote opposite, and the cones' shared edge reads
/// `same` weaker than the cycle through B's disc reads `opposite`. Taken
/// in index order the `same` vote came first and turned B's cone inside
/// out (flip[2] = true); strongest first it is the one dropped.
#[test]
fn a_weak_vote_contradicting_a_stronger_cycle_is_dropped() {
    let v = votes(&[
        ((0, 1), (0.0, 0.00565)),
        ((2, 3), (0.0, 0.00804)),
        ((0, 2), (0.00011, 0.0)),
        ((0, 3), (0.0, 0.00013)),
    ]);
    let (flip, components) = resolve_parity(4, &v);
    assert_eq!(flip, vec![Some(false); 4]);
    assert_eq!(components, vec![vec![0, 1, 2, 3]]);
}

#[test]
fn a_same_direction_vote_flips_the_second_face() {
    let v = votes(&[((0, 1), (0.3, 0.1)), ((1, 2), (0.0, 0.2))]);
    let (flip, components) = resolve_parity(3, &v);
    assert_eq!(flip, vec![Some(false), Some(true), Some(true)]);
    assert_eq!(components, vec![vec![0, 1, 2]]);
}

#[test]
fn tied_and_absent_votes_leave_faces_their_own_components() {
    let v = votes(&[((0, 2), (0.5, 0.5)), ((1, 3), (0.0, 1.0))]);
    let (flip, components) = resolve_parity(4, &v);
    assert_eq!(flip, vec![Some(false); 4]);
    assert_eq!(components, vec![vec![0], vec![1, 3], vec![2]]);
}
