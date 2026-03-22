//! Shared NURBS helper types and algorithms reused across curve, surface, and
//! modification operations (knot insertion, removal, degree elevation, splitting).

pub mod degree_elevation;
pub mod knot_insertion;
pub mod knot_removal;
pub mod splitting;

use crate::geometry::Point3;

/// A point in 4D homogeneous coordinates (wx, wy, wz, w).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Homogeneous {
    pub wx: f64,
    pub wy: f64,
    pub wz: f64,
    pub w: f64,
}

impl Homogeneous {
    pub fn from_point_weight(p: &Point3, w: f64) -> Self {
        Self {
            wx: w * p.x,
            wy: w * p.y,
            wz: w * p.z,
            w,
        }
    }

    pub fn lerp(a: &Self, b: &Self, alpha: f64) -> Self {
        let one_minus = 1.0 - alpha;
        Self {
            wx: one_minus * a.wx + alpha * b.wx,
            wy: one_minus * a.wy + alpha * b.wy,
            wz: one_minus * a.wz + alpha * b.wz,
            w: one_minus * a.w + alpha * b.w,
        }
    }

    pub fn to_point3(self) -> Point3 {
        if self.w.abs() < 1e-30 {
            Point3::new(0.0, 0.0, 0.0)
        } else {
            Point3::new(self.wx / self.w, self.wy / self.w, self.wz / self.w)
        }
    }

    #[allow(dead_code)]
    pub fn to_point3_and_weight(self) -> (Point3, f64) {
        (self.to_point3(), self.w)
    }
}

/// Find the knot span index `k` such that `knots[k] <= t < knots[k+1]`.
///
/// `degree`: spline degree
/// `knots`: knot vector
/// `n_controls`: number of control points
/// `t`: parameter value
///
/// Valid spans are `[degree, n_controls - 1]`.
pub(crate) fn find_span(degree: usize, knots: &[f64], n_controls: usize, t: f64) -> usize {
    let n = n_controls - 1;

    // Clamp to valid range
    if t >= knots[n + 1] {
        return n;
    }
    if t <= knots[degree] {
        return degree;
    }

    // Binary search
    let mut low = degree;
    let mut high = n + 1;
    let mut mid = (low + high) / 2;
    while t < knots[mid] || t >= knots[mid + 1] {
        if t < knots[mid] {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    mid
}

/// 1D de Boor evaluation in homogeneous coordinates.
///
/// `degree`: spline degree
/// `knots`: knot vector
/// `homogeneous_points`: homogeneous control points
/// `t`: parameter value
pub(crate) fn de_boor_1d(
    degree: usize,
    knots: &[f64],
    homogeneous_points: &[Homogeneous],
    t: f64,
) -> Homogeneous {
    let n = homogeneous_points.len();
    let k = find_span(degree, knots, n, t);

    // Initialize with the degree+1 relevant control points
    let mut d: Vec<Homogeneous> = (0..=degree)
        .map(|j| homogeneous_points[k - degree + j])
        .collect();

    // De Boor recursion
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let idx = k - degree + j;
            let denom = knots[idx + degree - r + 1] - knots[idx];
            let alpha = if denom.abs() < 1e-30 {
                0.0
            } else {
                (t - knots[idx]) / denom
            };
            d[j] = Homogeneous::lerp(&d[j - 1], &d[j], alpha);
        }
    }

    d[degree]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;

    #[test]
    fn homogeneous_zero_weight_returns_origin() {
        let h = Homogeneous {
            wx: 1.0,
            wy: 2.0,
            wz: 3.0,
            w: 0.0, // zero weight
        };
        let p = h.to_point3();
        assert_eq!(p.x, 0.0);
        assert_eq!(p.y, 0.0);
        assert_eq!(p.z, 0.0);
    }

    #[test]
    fn de_boor_with_repeated_knots() {
        // Degree 2, knots = [0,0,0,0,0,0] (all same) => denom = 0 => alpha = 0
        // This is degenerate but shouldn't panic
        let pts = vec![
            Homogeneous::from_point_weight(&Point3::new(0.0, 0.0, 0.0), 1.0),
            Homogeneous::from_point_weight(&Point3::new(1.0, 0.0, 0.0), 1.0),
            Homogeneous::from_point_weight(&Point3::new(2.0, 0.0, 0.0), 1.0),
        ];
        let knots = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let result = de_boor_1d(2, &knots, &pts, 0.0);
        // Should not panic; result is some point
        assert!(result.w.is_finite());
    }
}
