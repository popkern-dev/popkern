//! De Casteljau subdivision helpers for KP20 bi-6 patches.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::Kp20Bi6Patch;

pub(super) fn split_bi6_patch_u(patch: &Kp20Bi6Patch, t: f64) -> (Kp20Bi6Patch, Kp20Bi6Patch) {
    let mut left = Kp20Bi6Patch::zero();
    let mut right = Kp20Bi6Patch::zero();
    for column in 0..7 {
        let mut curve = [Point3::new(0.0, 0.0, 0.0); 7];
        for (row, point) in curve.iter_mut().enumerate() {
            *point = patch.cps[row][column];
        }
        let (left_curve, right_curve) = split_bi6_curve(curve, t);
        for row in 0..7 {
            left.cps[row][column] = left_curve[row];
            right.cps[row][column] = right_curve[row];
        }
    }
    (left, right)
}

pub(super) fn split_bi6_patch_v(patch: &Kp20Bi6Patch, t: f64) -> (Kp20Bi6Patch, Kp20Bi6Patch) {
    let mut lower = Kp20Bi6Patch::zero();
    let mut upper = Kp20Bi6Patch::zero();
    for row in 0..7 {
        let (lower_curve, upper_curve) = split_bi6_curve(patch.cps[row], t);
        lower.cps[row] = lower_curve;
        upper.cps[row] = upper_curve;
    }
    (lower, upper)
}

fn split_bi6_curve(mut levels: [Point3; 7], t: f64) -> ([Point3; 7], [Point3; 7]) {
    let mut left = [Point3::new(0.0, 0.0, 0.0); 7];
    let mut right = [Point3::new(0.0, 0.0, 0.0); 7];
    left[0] = levels[0];
    right[6] = levels[6];
    for degree in 1..=6 {
        for idx in 0..=(6 - degree) {
            levels[idx] = lerp(levels[idx], levels[idx + 1], t);
        }
        left[degree] = levels[0];
        right[6 - degree] = levels[6 - degree];
    }
    (left, right)
}

fn lerp(a: Point3, b: Point3, t: f64) -> Point3 {
    a + (b - a) * t
}
