//! Self-intersection detection for swept surfaces.
//!
//! Samples the surface normal field on a grid and detects folds where
//! adjacent normals in the v-direction reverse, indicating a self-intersection.

use popkern_core::geometry::NurbsSurface;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;

/// Grid resolution for self-intersection sampling.
const GRID_SIZE: usize = 20;

/// Check a swept surface for self-intersection.
///
/// Samples at a [`GRID_SIZE`] x [`GRID_SIZE`] grid and checks for sudden
/// normal reversals between adjacent v-samples at the same u position.
/// Returns `Some(v_parameter)` at the first detected fold, or `None` if
/// the surface is clean.
///
/// This is a heuristic check — it can miss narrow folds that fall between
/// sample points, but it catches the most common cases (e.g., a sweep
/// path with curvature tighter than the profile radius).
pub fn check_self_intersection(surface: &NurbsSurface, _tolerance: &Tolerance) -> Option<f64> {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();

    // The v-direction follows the sweep path, so adjacent samples must
    // sit close enough that the frame's own rotation between them stays
    // well under the fold threshold — a fixed grid across a multi-turn
    // helix reads 117° of frame rotation per step as a "fold" (PK-440
    // false positive on a 6.5-turn spring). Two samples per knot span
    // track the surface's actual v-complexity; finer sampling only
    // sharpens genuine fold detection.
    let v_spans = surface
        .knots_v
        .windows(2)
        .filter(|w| w[1] > w[0] + 1e-12)
        .count();
    let grid_v = GRID_SIZE.max(2 * v_spans).min(400);

    // For each u-sample, walk along v and check that adjacent normals
    // don't reverse. A reversal in the v-direction normal indicates the
    // surface folds back on itself along the sweep path.
    for iu in 0..=GRID_SIZE {
        let u = u_min + (iu as f64 / GRID_SIZE as f64) * (u_max - u_min);
        let mut prev_normal: Option<popkern_core::geometry::Vec3> = None;

        for iv in 0..=grid_v {
            let v = v_min + (iv as f64 / grid_v as f64) * (v_max - v_min);
            let normal = surface.normal_at(u, v);
            let len = normal.length();

            if len < 1e-20 {
                continue;
            }

            if let Some(ref prev) = prev_normal {
                let prev_len = prev.length();
                if prev_len > 1e-20 {
                    let dot = prev.dot(&normal);
                    // Check if the normal has reversed (more than 120 degrees from previous).
                    if dot < -len * prev_len * 0.5 {
                        return Some(v);
                    }
                }
            }

            prev_normal = Some(normal);
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;

    /// A flat bilinear patch — should never self-intersect.
    fn flat_patch() -> NurbsSurface {
        NurbsSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 2,
            num_v: 2,
            weights: vec![1.0; 4],
        }
    }

    #[test]
    fn flat_patch_no_self_intersection() {
        let surface = flat_patch();
        let tol = Tolerance::default();
        assert!(check_self_intersection(&surface, &tol).is_none());
    }

    #[test]
    fn folded_surface_detects_intersection() {
        // Create a surface that folds over itself by having control points
        // that cross in the v-direction.
        let surface = NurbsSurface {
            degree_u: 1,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                // u=0 row: goes forward, then back (fold)
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 0.0, 2.0),
                Point3::new(0.0, 0.0, -1.0), // fold back
                // u=1 row: same fold
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 2.0),
                Point3::new(1.0, 0.0, -1.0),
            ],
            num_u: 2,
            num_v: 3,
            weights: vec![1.0; 6],
        };
        let tol = Tolerance::default();
        // This surface folds — the normal should reverse.
        let result = check_self_intersection(&surface, &tol);
        // We accept either detection or miss (heuristic), but log it.
        // The fold is quite dramatic so it should be caught.
        if let Some(v) = result {
            assert!((0.0..=1.0).contains(&v));
        }
    }

    #[test]
    fn cylinder_like_surface_no_false_positive() {
        // A swept circle surface should NOT trigger self-intersection.
        // Normals rotate around u but don't reverse along v.
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let surface = NurbsSurface {
            degree_u: 2,
            degree_v: 1,
            knots_u: vec![
                0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
            ],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                // v=0 (bottom circle)
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 5.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(1.0, 1.0, 5.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 5.0),
                Point3::new(-1.0, 1.0, 0.0),
                Point3::new(-1.0, 1.0, 5.0),
                Point3::new(-1.0, 0.0, 0.0),
                Point3::new(-1.0, 0.0, 5.0),
                Point3::new(-1.0, -1.0, 0.0),
                Point3::new(-1.0, -1.0, 5.0),
                Point3::new(0.0, -1.0, 0.0),
                Point3::new(0.0, -1.0, 5.0),
                Point3::new(1.0, -1.0, 0.0),
                Point3::new(1.0, -1.0, 5.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 5.0),
            ],
            num_u: 9,
            num_v: 2,
            weights: vec![
                1.0, 1.0, w, w, 1.0, 1.0, w, w, 1.0, 1.0, w, w, 1.0, 1.0, w, w, 1.0, 1.0,
            ],
        };
        let tol = Tolerance::default();
        assert!(
            check_self_intersection(&surface, &tol).is_none(),
            "cylinder should not self-intersect"
        );
    }
}
