//! Knot removal with tolerance control (Piegl & Tiller A5.8).

use super::Homogeneous;
use crate::geometry::Point3;

/// Attempt to remove the knot at `knot_index` from a 1D NURBS.
/// Returns Some((new_knots, new_points, new_weights)) if removal deviation < tolerance.
/// Returns None if removal would exceed tolerance or index is invalid.
/// Use tolerance = 0.0 for exact removal only.
#[allow(dead_code)]
pub(crate) fn remove_knot_1d(
    degree: usize,
    knots: &[f64],
    control_points: &[Point3],
    weights: &[f64],
    knot_index: usize,
    tolerance: f64,
) -> Option<(Vec<f64>, Vec<Point3>, Vec<f64>)> {
    let n = control_points.len();
    let p = degree;
    let m = knots.len();

    // Validate: knot_index must refer to an interior knot
    if knot_index < p + 1 || knot_index >= m - p - 1 {
        return None;
    }
    if n <= p + 1 {
        return None;
    }

    let u = knots[knot_index];

    // Find the last occurrence of this knot value (Piegl & Tiller convention)
    let mut r = knot_index;
    while r + 1 < m && (knots[r + 1] - u).abs() < 1e-15 {
        r += 1;
    }
    if r >= m - p - 1 {
        return None;
    }

    // Build the post-removal knot vector (removing one instance of u at index r).
    // This equals the pre-insertion knot vector, so we use it to compute the
    // insertion alphas needed to reverse the blending.
    let mut new_knots = Vec::with_capacity(m - 1);
    new_knots.extend_from_slice(&knots[..r]);
    new_knots.extend_from_slice(&knots[r + 1..]);

    // Convert to homogeneous coordinates
    let pw: Vec<Homogeneous> = control_points
        .iter()
        .zip(weights.iter())
        .map(|(pt, &w)| Homogeneous::from_point_weight(pt, w))
        .collect();

    // Find the span k in the post-removal knot vector such that
    // new_knots[k] <= u < new_knots[k+1].
    let k = {
        let mut span = p;
        for i in p..n - 1 {
            if new_knots[i + 1] > u + 1e-15 {
                span = i;
                break;
            }
            span = i;
        }
        span
    };

    // The insertion affected current control points at indices first..=last.
    let first = if k >= p { k - p + 1 } else { 1 };
    let last = k;

    if first == 0 || last + 1 >= n {
        return None;
    }

    let num = last - first + 1; // number of affected current control points

    // Left sweep: recover pre-insertion points by solving insertion equations
    // from left to right. current[i] = (1-a_i)*old[i-1] + a_i*old[i]
    // => old[i] = (current[i] - (1-a_i)*old[i-1]) / a_i
    let mut left_cand = vec![hom_zero(); num];
    for j in 0..num {
        let idx = first + j;
        let a = insertion_alpha(u, idx, p, &new_knots);
        let prev = if j == 0 {
            pw[first - 1]
        } else {
            left_cand[j - 1]
        };
        if a.abs() < 1e-30 {
            left_cand[j] = pw[idx];
        } else {
            left_cand[j] = hom_solve_right(&pw[idx], &prev, a);
        }
    }

    // Right sweep: solve from right to left.
    // current[i] = (1-a_i)*old[i-1] + a_i*old[i]
    // => old[i-1] = (current[i] - a_i*old[i]) / (1-a_i)
    // The right anchor is pw[last+1] which equals old[k] (shifted).
    let mut right_cand = vec![hom_zero(); num];
    for j in 0..num {
        let idx = last - j;
        let a = insertion_alpha(u, idx, p, &new_knots);
        let one_a = 1.0 - a;
        let next = if j == 0 {
            pw[last + 1]
        } else {
            right_cand[j - 1]
        };
        if one_a.abs() < 1e-30 {
            right_cand[j] = pw[idx];
        } else {
            right_cand[j] = hom_solve_left(&pw[idx], &next, a);
        }
    }

    // Removability check.
    // left_cand[j] = candidate for old[first+j] (forward order)
    // right_cand[j] = candidate for old[last-j-1] (reverse order)
    // n_repl = num - 1 replacement points needed.
    let n_repl = num - 1;

    let removable = if n_repl == 0 {
        // Single affected point: check it's reconstructable from its neighbors.
        let a = insertion_alpha(u, first, p, &new_knots);
        let recon = Homogeneous::lerp(&pw[first - 1], &pw[last + 1], a);
        pw[first].to_point3().distance_to(&recon.to_point3()) <= tolerance
    } else {
        // Compare left and right candidates at the center replacement point.
        let center = n_repl / 2;
        let lp = left_cand[center].to_point3();
        let rp = right_cand[n_repl - 1 - center].to_point3();
        lp.distance_to(&rp) <= tolerance
    };

    if !removable {
        return None;
    }

    // Build new control points (n-1 total).
    let mut new_hw = Vec::with_capacity(n - 1);

    // Unaffected left
    new_hw.extend_from_slice(&pw[..first]);

    // Replacement points from the two sweeps.
    if n_repl > 0 {
        if n_repl % 2 == 1 {
            // Odd: left half, averaged center, right half
            let center = n_repl / 2;
            new_hw.extend_from_slice(&left_cand[..center]);
            let lp = &left_cand[center];
            let rp = &right_cand[n_repl - 1 - center];
            new_hw.push(Homogeneous {
                wx: 0.5 * (lp.wx + rp.wx),
                wy: 0.5 * (lp.wy + rp.wy),
                wz: 0.5 * (lp.wz + rp.wz),
                w: 0.5 * (lp.w + rp.w),
            });
            for j in (center + 1)..n_repl {
                new_hw.push(right_cand[n_repl - 1 - j]);
            }
        } else {
            // Even: left half, right half
            let half = n_repl / 2;
            new_hw.extend_from_slice(&left_cand[..half]);
            for j in half..n_repl {
                new_hw.push(right_cand[n_repl - 1 - j]);
            }
        }
    }

    // Unaffected right
    new_hw.extend_from_slice(&pw[last + 1..]);

    if new_hw.len() != n - 1 {
        return None;
    }

    // Convert back from homogeneous
    let mut new_points = Vec::with_capacity(new_hw.len());
    let mut new_weights = Vec::with_capacity(new_hw.len());
    for h in &new_hw {
        let (pt, w) = h.to_point3_and_weight();
        new_points.push(pt);
        new_weights.push(w);
    }

    Some((new_knots, new_points, new_weights))
}

/// Compute the insertion alpha: (u - knots[i]) / (knots[i+p] - knots[i])
/// using the post-removal (= pre-insertion) knot vector.
fn insertion_alpha(u: f64, i: usize, p: usize, new_knots: &[f64]) -> f64 {
    if i + p >= new_knots.len() {
        return 0.0;
    }
    let denom = new_knots[i + p] - new_knots[i];
    if denom.abs() < 1e-30 {
        0.0
    } else {
        (u - new_knots[i]) / denom
    }
}

fn hom_zero() -> Homogeneous {
    Homogeneous {
        wx: 0.0,
        wy: 0.0,
        wz: 0.0,
        w: 0.0,
    }
}

/// Solve for Q[i] from: P[i] = (1-a)*Q[i-1] + a*Q[i]
fn hom_solve_right(p_old: &Homogeneous, prev: &Homogeneous, a: f64) -> Homogeneous {
    let one_a = 1.0 - a;
    Homogeneous {
        wx: (p_old.wx - one_a * prev.wx) / a,
        wy: (p_old.wy - one_a * prev.wy) / a,
        wz: (p_old.wz - one_a * prev.wz) / a,
        w: (p_old.w - one_a * prev.w) / a,
    }
}

/// Solve for Q[i-1] from: P[i] = (1-a)*Q[i-1] + a*Q[i]
fn hom_solve_left(p_old: &Homogeneous, next: &Homogeneous, a: f64) -> Homogeneous {
    let one_a = 1.0 - a;
    Homogeneous {
        wx: (p_old.wx - a * next.wx) / one_a,
        wy: (p_old.wy - a * next.wy) / one_a,
        wz: (p_old.wz - a * next.wz) / one_a,
        w: (p_old.w - a * next.w) / one_a,
    }
}

/// Remove knot up to `times` times, returning result and actual count.
#[allow(dead_code)]
pub(crate) fn remove_knot_times_1d(
    degree: usize,
    knots: &[f64],
    control_points: &[Point3],
    weights: &[f64],
    knot_index: usize,
    times: usize,
    tolerance: f64,
) -> (Vec<f64>, Vec<Point3>, Vec<f64>, usize) {
    let mut cur_knots = knots.to_vec();
    let mut cur_pts = control_points.to_vec();
    let mut cur_wts = weights.to_vec();
    let mut count = 0;
    let cur_index = knot_index;

    for _ in 0..times {
        match remove_knot_1d(degree, &cur_knots, &cur_pts, &cur_wts, cur_index, tolerance) {
            Some((nk, np, nw)) => {
                cur_knots = nk;
                cur_pts = np;
                cur_wts = nw;
                count += 1;
                if cur_index >= cur_knots.len() {
                    break;
                }
            }
            None => break,
        }
    }

    (cur_knots, cur_pts, cur_wts, count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;
    use crate::geometry::nurbs_ops::knot_insertion::insert_knot_1d;

    #[test]
    fn remove_inserted_knot_roundtrip() {
        let pts = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        let wts = vec![1.0, 1.0];
        let knots = vec![0.0, 0.0, 1.0, 1.0];

        let (k2, p2, w2) = insert_knot_1d(1, &knots, &pts, &wts, 0.5);
        assert_eq!(p2.len(), 3);

        let result = remove_knot_1d(1, &k2, &p2, &w2, 2, 0.0);
        assert!(result.is_some());
        let (k3, p3, _w3) = result.unwrap();
        assert_eq!(p3.len(), 2);
        assert_eq!(k3.len(), 4);
    }

    #[test]
    fn cannot_remove_boundary_knot() {
        let pts = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        let wts = vec![1.0, 1.0];
        let knots = vec![0.0, 0.0, 1.0, 1.0];
        assert!(remove_knot_1d(1, &knots, &pts, &wts, 0, 0.0).is_none());
        assert!(remove_knot_1d(1, &knots, &pts, &wts, 3, 0.0).is_none());
    }

    #[test]
    fn remove_knot_with_tolerance() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.3, 0.5, 0.0),
            Point3::new(0.7, 0.5, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0, 1.0];
        let knots = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];

        let result = remove_knot_1d(2, &knots, &pts, &wts, 3, 10.0);
        assert!(result.is_some());
    }

    #[test]
    fn remove_knot_times() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0];
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];

        let (k2, p2, w2) = insert_knot_1d(2, &knots, &pts, &wts, 0.5);
        let (k3, p3, w3) = insert_knot_1d(2, &k2, &p2, &w2, 0.5);
        let (_kr, pr, _wr, count) = remove_knot_times_1d(2, &k3, &p3, &w3, 3, 3, 0.0);
        assert_eq!(count, 2);
        assert_eq!(pr.len(), 3);
    }
}
