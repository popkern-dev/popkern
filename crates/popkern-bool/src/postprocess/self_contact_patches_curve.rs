//! An edge curve proved apart from a NURBS face's patches, so the
//! contact search need not walk it (PK-1172).

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::nurbs_curve::NurbsCurve;

use super::{DEPTH, Patch};
use crate::parallel::ff_parallel::ff_parallel_reach::box_gap;

/// Whether `curve`, of positive weights, provably lies farther than
/// `pad` from the patches `theirs`: each piece of it, halved while
/// near, has a control box (holding the piece) farther than `pad` from
/// each patch, halved while near.
pub(in super::super) fn curve_apart(
    curve: &NurbsCurve,
    theirs: &[Patch],
    pad: f64,
    depth: usize,
) -> bool {
    if curve.weights.iter().any(|&w| w <= 0.0) {
        return false;
    }
    let Some(bx) = Aabb::from_points(&curve.control_points) else {
        return false;
    };
    if theirs.iter().all(|b| box_apart(&bx, b, pad, DEPTH)) {
        return true;
    }
    if depth == 0 {
        return false;
    }
    let p = curve.degree;
    let (Some(&t0), Some(&t1)) = (
        curve.knots.get(p),
        curve
            .knots
            .len()
            .checked_sub(p + 1)
            .and_then(|i| curve.knots.get(i)),
    ) else {
        return false;
    };
    curve.split_at(0.5 * (t0 + t1)).is_ok_and(|(l, r)| {
        curve_apart(&l, theirs, pad, depth - 1) && curve_apart(&r, theirs, pad, depth - 1)
    })
}

/// Whether the box `bx` lies farther than `pad` from patch `b`, halving
/// the patch while its box is near.
fn box_apart(bx: &Aabb, b: &Patch, pad: f64, depth: usize) -> bool {
    if box_gap(bx, &b.bx) > pad {
        return true;
    }
    if depth == 0 {
        return false;
    }
    let (eu, ev) = b.extents();
    b.halves(eu >= ev).is_some_and(|(x, y)| {
        box_apart(bx, &x, pad, depth - 1) && box_apart(bx, &y, pad, depth - 1)
    })
}

#[cfg(test)]
#[path = "self_contact_patches_curve_tests.rs"]
mod tests;
