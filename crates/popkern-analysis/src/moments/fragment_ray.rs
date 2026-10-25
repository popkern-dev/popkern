//! Rays against a cell's surface pieces: where material beside a
//! surface ends (its own opposite surface), and whether the straight way
//! between two points inside the material leaves it.

use crate::clip::FacetPiece;

/// The smallest distance `t` in (`tol`, `max`) at which
/// `origin + t dir` (`dir` a unit vector) meets one of `pieces` (planar
/// convex polygons, met on their boundary too, to `tol`), skipping index
/// `skip`, and that piece; `max` when none. Pieces whose plane holds the
/// ray are skipped.
pub(super) fn first_hit(
    pieces: &[FacetPiece],
    skip: Option<usize>,
    origin: [f64; 3],
    dir: [f64; 3],
    (tol, max): (f64, f64),
) -> (f64, Option<usize>) {
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let (mut best, mut hit) = (max, None);
    for (i, piece) in pieces.iter().enumerate() {
        if Some(i) == skip {
            continue;
        }
        let v = piece.polygon.vertices();
        let Some(&a) = v.first() else {
            continue;
        };
        let n = piece.normal;
        let along = dot(n, dir);
        if along.abs() <= 1.0e-12 {
            continue;
        }
        let t = dot(n, sub(a, origin)) / along;
        if t <= tol || t >= best {
            continue;
        }
        let x = [0, 1, 2].map(|k| origin[k] + t * dir[k]);
        // Inside (or on the edge of) the convex polygon: x on the same
        // side of every edge, to a tolerance on the edge's length.
        let side = |e: usize| {
            let (p, q) = (v[e], v[(e + 1) % v.len()]);
            let edge = sub(q, p);
            dot(cross(edge, sub(x, p)), n) / dot(edge, edge).sqrt().max(f64::MIN_POSITIVE)
        };
        if (0..v.len()).all(|e| side(e) >= -tol) || (0..v.len()).all(|e| side(e) <= tol) {
            (best, hit) = (t, Some(i));
        }
    }
    (best, hit)
}
