//! Knot refinement: insert many knots in one pass.
//!
//! Repeated single insertion (`knot_insertion::insert_knots_1d`)
//! rebuilds the whole control polygon for every knot, so refining
//! every span of an n-span curve costs O(n²). The refinement below
//! (Piegl & Tiller, *The NURBS Book*, A5.4) builds the refined polygon
//! once, right to left, in O((n + r)·p) for r inserted knots. Measured
//! (PK-786): the radial no-crossing certificate refines a thread wall
//! up to three times; by single insertion that made the rod wall x
//! thread flank intersections quadratic in turns (4.1 s at 40 turns,
//! 16.6 s at 80, single-threaded), and the M6 x 80 rod's intersection
//! stage drops from 15.4 s to 9.9 s with this refinement.

use super::{Homogeneous, find_span};
use crate::geometry::Point3;

/// Insert the sorted knots `inserts` into a degree-`degree` rational
/// curve given by `knots` and homogeneous control points `points`.
/// Returns the refined knot vector and control points; the curve is
/// unchanged. `None` when the inputs are not a valid curve (too few
/// points for the degree, a knot vector of the wrong length) or an
/// insert lies outside the knot domain.
pub(crate) fn refine_knot_vector(
    degree: usize,
    knots: &[f64],
    points: &[Homogeneous],
    inserts: &[f64],
) -> Option<(Vec<f64>, Vec<Homogeneous>)> {
    let p = degree;
    let count = points.len();
    if count <= p || knots.len() != count + p + 1 {
        return None;
    }
    if inserts.is_empty() {
        return Some((knots.to_vec(), points.to_vec()));
    }
    let (lo, hi) = (knots[p], knots[count]);
    if inserts.iter().any(|x| !(*x >= lo && *x <= hi)) || inserts.windows(2).any(|w| w[1] < w[0]) {
        return None;
    }
    let n = count - 1;
    let m = n + p + 1;
    let r = inserts.len() - 1;
    let a = find_span(p, knots, count, inserts[0]);
    let b = find_span(p, knots, count, inserts[r]) + 1;

    let zero = Homogeneous {
        wx: 0.0,
        wy: 0.0,
        wz: 0.0,
        w: 0.0,
    };
    let mut q = vec![zero; count + r + 1];
    let mut u_bar = vec![0.0; m + r + 2];
    q[..=(a - p)].copy_from_slice(&points[..=(a - p)]);
    for j in (b - 1)..=n {
        q[j + r + 1] = points[j];
    }
    u_bar[..=a].copy_from_slice(&knots[..=a]);
    for j in (b + p)..=m {
        u_bar[j + r + 1] = knots[j];
    }

    let mut i = b + p - 1;
    let mut k = b + p + r;
    for &x in inserts.iter().rev() {
        while x <= knots[i] && i > a {
            q[k - p - 1] = points[i - p - 1];
            u_bar[k] = knots[i];
            k -= 1;
            i -= 1;
        }
        q[k - p - 1] = q[k - p];
        for l in 1..=p {
            let ind = k - p + l;
            let num = u_bar[k + l] - x;
            if num == 0.0 {
                q[ind - 1] = q[ind];
            } else {
                let alpha = num / (u_bar[k + l] - knots[i - p + l]);
                q[ind - 1] = Homogeneous::lerp(&q[ind], &q[ind - 1], alpha);
            }
        }
        u_bar[k] = x;
        k -= 1;
    }
    Some((u_bar, q))
}

/// [`refine_knot_vector`] on Euclidean points with weights.
pub(crate) fn refine_knot_vector_points(
    degree: usize,
    knots: &[f64],
    points: &[Point3],
    weights: &[f64],
    inserts: &[f64],
) -> Option<(Vec<f64>, Vec<Point3>, Vec<f64>)> {
    if points.len() != weights.len() {
        return None;
    }
    let hom: Vec<Homogeneous> = points
        .iter()
        .zip(weights)
        .map(|(pt, &w)| Homogeneous::from_point_weight(pt, w))
        .collect();
    let (knots, refined) = refine_knot_vector(degree, knots, &hom, inserts)?;
    let (points, weights) = refined
        .into_iter()
        .map(|h| h.to_point3_and_weight())
        .unzip();
    Some((knots, points, weights))
}

#[cfg(test)]
#[path = "knot_refinement_tests.rs"]
mod tests;
