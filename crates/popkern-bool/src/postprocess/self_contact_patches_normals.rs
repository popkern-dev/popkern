//! Whether a NURBS face's normals can be perpendicular to a direction
//! (PK-1175): a cylinder's normals all are to its axis, so where none of
//! a face's are, the two have no parallel normals, and their nearest
//! points lie on a boundary of one (faces of a valid shell do not
//! cross), found by the edge walks.
//!
//! On a polynomial Bézier patch of `nu x nv` points, `S_u` is a positive
//! combination of the net's u-differences `d = P(i+1, j) - P(i, j)`, and
//! `S_v` of its v-differences `e`, so `det(S_u, S_v, a) = (S_u x S_v) · a`
//! is a positive combination of the products `det(d, e, a)`. When those
//! all have one strict sign, the normal is never perpendicular to `a`
//! on the patch. A patch whose products differ in sign is halved along
//! its longer direction (its differences close on the derivatives), to
//! a depth. A rational patch (weights not all equal) is not decided.

use popkern_core::geometry::Vec3;

use super::Patch;

/// Halvings of a patch whose difference products differ in sign.
const NORMAL_DEPTH: usize = 6;

/// Whether no normal of the patches is perpendicular to `axis` (unit).
pub(in crate::postprocess) fn normals_never_across(patches: &[Patch], axis: Vec3) -> bool {
    patches
        .iter()
        .all(|p| patch_never_across(p, axis, NORMAL_DEPTH))
}

fn patch_never_across(p: &Patch, axis: Vec3, depth: usize) -> bool {
    match product_sign(p, axis) {
        Some(true) => true,
        Some(false) => false,
        None if depth == 0 => false,
        None => {
            let (eu, ev) = p.extents();
            p.halves(eu >= ev).is_some_and(|(x, y)| {
                patch_never_across(&x, axis, depth - 1) && patch_never_across(&y, axis, depth - 1)
            })
        }
    }
}

/// `Some(true)` when the difference products of `p` share one strict
/// sign, `Some(false)` when `p` is rational (not decided by them), and
/// `None` when they differ in sign.
fn product_sign(p: &Patch, axis: Vec3) -> Option<bool> {
    let w0 = p.net[0][3];
    if p.net.iter().any(|h| (h[3] - w0).abs() > 1e-12 * w0.abs()) {
        return Some(false);
    }
    let at = |a: usize, b: usize| p.points[a * p.nv + b];
    let du: Vec<Vec3> = (0..p.nu - 1)
        .flat_map(|a| (0..p.nv).map(move |b| (a, b)))
        .map(|(a, b)| at(a + 1, b) - at(a, b))
        .collect();
    let dv: Vec<Vec3> = (0..p.nu)
        .flat_map(|a| (0..p.nv - 1).map(move |b| (a, b)))
        .map(|(a, b)| at(a, b + 1) - at(a, b))
        .collect();
    let (mut pos, mut neg) = (false, false);
    for d in &du {
        let da = d.cross(&axis);
        for e in &dv {
            // det(d, e, a) = (d x e) · a = (a x d) · e = -(d x a) · e.
            let x = -da.dot(e);
            if x > 0.0 {
                pos = true;
            } else if x < 0.0 {
                neg = true;
            } else {
                return None;
            }
            if pos && neg {
                return None;
            }
        }
    }
    Some(true)
}

#[cfg(test)]
#[path = "self_contact_patches_normals_tests.rs"]
mod tests;
