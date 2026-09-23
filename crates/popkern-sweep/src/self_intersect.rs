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
    // The span-tracking grid is cost-capped, so on a long coil a
    // single coarse step can legitimately carry more than the fold
    // threshold of SMOOTH frame rotation (a 150-turn thread tool at
    // the historical 400-sample cap stepped 135° — a false
    // SelfIntersection decline, PK-495). A flagged step is therefore
    // refined before it convicts: smooth rotation halves with the
    // interval, a genuine fold's normal flip does not. The cap keeps
    // the un-aliased band: below it a coarse step stays under 240°
    // (the dot test flags [120°, 240°]), so every fold-carrying step
    // is flagged and handed to refinement — 2 per span holds that to
    // ~1365 helix turns at the skinning density, far past any served
    // class.
    let grid_v = GRID_SIZE.max(2 * v_spans).min(2048);

    // For each u-sample, walk along v and check that adjacent normals
    // don't reverse. A reversal in the v-direction normal indicates the
    // surface folds back on itself along the sweep path.
    for iu in 0..=GRID_SIZE {
        let u = u_min + (iu as f64 / GRID_SIZE as f64) * (u_max - u_min);
        let mut prev: Option<(f64, popkern_core::geometry::Vec3)> = None;

        for iv in 0..=grid_v {
            let v = v_min + (iv as f64 / grid_v as f64) * (v_max - v_min);
            let normal = surface.normal_at(u, v);
            let len = normal.length();

            if len < 1e-20 {
                continue;
            }

            if let Some((prev_v, prev_n)) = prev {
                let prev_len = prev_n.length();
                // Reversal (more than 120° between adjacent normals)
                // marks the step as suspicious; refinement decides
                // whether it is a fold or accumulated smooth rotation.
                if prev_len > 1e-20
                    && prev_n.dot(&normal) < -len * prev_len * 0.5
                    && fold_survives_refinement(surface, u, prev_v, v)
                {
                    return Some(v);
                }
            }

            prev = Some((v, normal));
        }
    }

    None
}

/// Whether a normal reversal across `[v0, v1]` at `u` persists under
/// interval bisection. Smooth frame rotation spreads evenly, so each
/// half of a bisected interval carries about half the angle and the
/// reversal dissolves within a few levels; a genuine fold is a
/// normal FLIP at a point, and the half containing it keeps the
/// reversal at every level. Depth 48 halvings put the interval at
/// the f64 resolution floor of any practical knot domain — a
/// reversal still alive there is a fold.
fn fold_survives_refinement(surface: &NurbsSurface, u: f64, v0: f64, v1: f64) -> bool {
    let (mut lo, mut hi) = (v0, v1);
    for _ in 0..48 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            return true;
        }
        let n_lo = surface.normal_at(u, lo);
        let n_mid = surface.normal_at(u, mid);
        let n_hi = surface.normal_at(u, hi);
        let reversed =
            |a: &popkern_core::geometry::Vec3, b: &popkern_core::geometry::Vec3| -> bool {
                let (la, lb) = (a.length(), b.length());
                la > 1e-20 && lb > 1e-20 && a.dot(b) < -la * lb * 0.5
            };
        if n_mid.length() < 1e-20 {
            // Degenerate normal mid-interval: cannot separate the
            // halves — keep the conservative verdict.
            return true;
        }
        if reversed(&n_lo, &n_mid) {
            hi = mid;
        } else if reversed(&n_mid, &n_hi) {
            lo = mid;
        } else {
            // Neither half reverses: the coarse step's angle was
            // accumulated smooth rotation, not a flip.
            return false;
        }
    }
    true
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
        // The fold is a genuine normal FLIP, so it must survive the
        // refinement pass (PK-495) — a miss here would mean the
        // refinement dissolves real folds, not just smooth rotation.
        let v = check_self_intersection(&surface, &tol).expect("fold must be detected");
        assert!((0.0..=1.0).contains(&v));
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
