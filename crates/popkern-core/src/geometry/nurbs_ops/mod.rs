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

/// Maximum NURBS degree supported by the stack-allocated de Boor
/// scratch buffer. Practical NURBS curves and surfaces in popkern
/// stay well below this — sweep spines are typically degree 3, G2
/// fillet cross-sections are degree 5, and STEP exchange rarely
/// exceeds degree 7. Raising this bumps every de Boor call's stack
/// footprint by 32 bytes per slot.
pub(crate) const MAX_DE_BOOR_DEGREE: usize = 32;

/// 1D de Boor evaluation in homogeneous coordinates.
///
/// `degree`: spline degree (must be < [`MAX_DE_BOOR_DEGREE`])
/// `knots`: knot vector
/// `homogeneous_points`: homogeneous control points
/// `t`: parameter value
///
/// Uses a fixed-size on-stack scratch buffer so point evaluation is
/// allocation-free. Prefer [`de_boor_1d_points`] when the caller has
/// `&[Point3]` and `&[f64]` — it avoids the `Vec<Homogeneous>` temp
/// the caller would otherwise have to build.
pub(crate) fn de_boor_1d(
    degree: usize,
    knots: &[f64],
    homogeneous_points: &[Homogeneous],
    t: f64,
) -> Homogeneous {
    debug_assert!(
        degree < MAX_DE_BOOR_DEGREE,
        "de_boor_1d: degree {degree} exceeds MAX_DE_BOOR_DEGREE {MAX_DE_BOOR_DEGREE}"
    );
    let n = homogeneous_points.len();
    let k = find_span(degree, knots, n, t);

    // Stack scratch for the active degree+1 control points.
    let mut d = [Homogeneous {
        wx: 0.0,
        wy: 0.0,
        wz: 0.0,
        w: 0.0,
    }; MAX_DE_BOOR_DEGREE];
    for j in 0..=degree {
        d[j] = homogeneous_points[k - degree + j];
    }

    // De Boor recursion (in-place on the first degree+1 slots).
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

/// 1D de Boor evaluation directly from Cartesian control points
/// plus weights.
///
/// Does the homogeneous conversion on-the-fly into the stack scratch
/// buffer, so the caller never materialises a `Vec<Homogeneous>`.
/// Called from the NurbsCurve hot path (`point_at`) where every
/// avoided allocation shows up in profiling.
pub(crate) fn de_boor_1d_points(
    degree: usize,
    knots: &[f64],
    control_points: &[Point3],
    weights: &[f64],
    t: f64,
) -> Homogeneous {
    debug_assert!(
        degree < MAX_DE_BOOR_DEGREE,
        "de_boor_1d_points: degree {degree} exceeds MAX_DE_BOOR_DEGREE {MAX_DE_BOOR_DEGREE}"
    );
    debug_assert_eq!(
        control_points.len(),
        weights.len(),
        "de_boor_1d_points: control_points and weights length mismatch"
    );
    let n = control_points.len();
    let k = find_span(degree, knots, n, t);

    let mut d = [Homogeneous {
        wx: 0.0,
        wy: 0.0,
        wz: 0.0,
        w: 0.0,
    }; MAX_DE_BOOR_DEGREE];
    for (j, slot) in d.iter_mut().take(degree + 1).enumerate() {
        let idx = k - degree + j;
        *slot = Homogeneous::from_point_weight(&control_points[idx], weights[idx]);
    }

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
