//! The shadow comparison names how each differing fragment differs:
//! only a membership difference changes the served point set.

use super::difference;
use crate::complex_shadow::ShadowDifference;
use crate::merge::fragment_classify::FragmentSelection;

fn sel(keep: bool, flip: bool) -> FragmentSelection {
    FragmentSelection { keep, flip }
}

#[test]
fn a_flip_on_a_fragment_both_keep_is_orientation() {
    let ours = [sel(true, true)];
    let theirs = [sel(true, false)];
    let got = difference(0, Some(&[0]), &ours, &theirs);
    assert_eq!(got, ShadowDifference::Orientation);
}

#[test]
fn a_fragment_no_complex_face_lifts_is_a_hole_filler() {
    let ours = [sel(false, false)];
    let theirs = [sel(true, false)];
    assert_eq!(
        difference(0, None, &ours, &theirs),
        ShadowDifference::HoleFiller
    );
}

/// Fragments 0 and 1 are one operand's seam-split copy of a patch,
/// fragment 2 the partner's whole copy: each classifier keeps one copy.
#[test]
fn another_copy_of_a_kept_patch_is_a_coincident_copy() {
    let ours = [sel(false, false), sel(false, false), sel(true, false)];
    let theirs = [sel(true, false), sel(true, false), sel(false, false)];
    for i in 0..3 {
        assert_eq!(
            difference(i, Some(&[2, 0, 1]), &ours, &theirs),
            ShadowDifference::CoincidentCopy
        );
    }
}

/// The patch kept by one classifier and by neither copy in the other,
/// and a lone fragment kept by one only: membership.
#[test]
fn a_patch_or_fragment_only_one_keeps_is_membership() {
    let ours = [sel(false, false), sel(true, false), sel(true, false)];
    let theirs = [sel(false, false), sel(false, false), sel(false, false)];
    assert_eq!(
        difference(1, Some(&[1, 0]), &ours, &theirs),
        ShadowDifference::Membership
    );
    assert_eq!(
        difference(2, Some(&[2]), &ours, &theirs),
        ShadowDifference::Membership
    );
}
