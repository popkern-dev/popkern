//! Constructor for [`RectangularGregoryTwoStripe`] —
//! promotion from [`super::RectangularGregoryBiseptic`].
//!
//! Seeds the 4 newly-free corner-adjacent positions
//! `(1,1), (1,6), (6,1), (6,6)` by evaluating the biseptic's
//! rational-blend at a representative `(u, v)` corresponding
//! to each position's Bezier index — this is the value the
//! biseptic surface has AT the corner-adjacent position's
//! parameter when the rational blend is active.
//!
//! For the `(i, j)` Bezier index, the representative
//! parameter is `(i/7, j/7)` — the "natural" location of
//! that CP in parameter space. The rational blend evaluated
//! there returns a `Point3` that the two-stripe seeds into
//! `grid[i][j]`. The rest of the grid carries over verbatim.
//!
//! After this seeding the two-stripe patch agrees with the
//! biseptic patch ONLY at the corner and boundary positions
//! (not in the interior — the rational blend varies with
//! (u, v) while the plain Bezier does not). The spoke-G1
//! solver's later two-stripe refinement passes refine the interior.

use super::rectangular_biseptic::RectangularGregoryBiseptic;
use super::rectangular_gregory_two_stripe::RectangularGregoryTwoStripe;

impl RectangularGregoryTwoStripe {
    /// Promote a [`RectangularGregoryBiseptic`] patch to
    /// two-stripe form: take its `grid` verbatim and replace
    /// the 4 corner-adjacent rational-blend placeholders
    /// `(1,1), (1,6), (6,1), (6,6)` with the blend's value
    /// evaluated at the corresponding Bezier-index parameter
    /// `(i/7, j/7)`.
    pub fn from_biseptic(src: &RectangularGregoryBiseptic) -> Self {
        let mut grid = src.grid;
        for (i, j) in [(1usize, 1usize), (1, 6), (6, 1), (6, 6)] {
            let u = (i as f64) / 7.0;
            let v = (j as f64) / 7.0;
            grid[i][j] = evaluate_corner_adjacent_blend(src, i, j, u, v);
        }
        RectangularGregoryTwoStripe { grid }
    }
}

/// Reconstructs what the biseptic's rational blend evaluates
/// to at `(u, v)` for corner-adjacent index `(i, j)`. Mirrors
/// the private `corner_adjacent_blend` helper in
/// [`super::rectangular_biseptic`]; duplicating the tiny
/// formula here avoids widening that module's public surface.
fn evaluate_corner_adjacent_blend(
    src: &RectangularGregoryBiseptic,
    i: usize,
    j: usize,
    u: f64,
    v: f64,
) -> crate::geometry::Point3 {
    let (k, wu, wv) = match (i, j) {
        (1, 1) => (0, u, v),
        (1, 6) => (1, u, 1.0 - v),
        (6, 6) => (2, 1.0 - u, 1.0 - v),
        (6, 1) => (3, 1.0 - u, v),
        _ => unreachable!(),
    };
    let (p_u, p_v) = src.corner_inner_cps[k];
    let denom = wu + wv;
    if denom < 1e-15 {
        return crate::geometry::Point3::new(
            0.5 * (p_u.x + p_v.x),
            0.5 * (p_u.y + p_v.y),
            0.5 * (p_u.z + p_v.z),
        );
    }
    let c_u = wv / denom;
    let c_v = wu / denom;
    crate::geometry::Point3::new(
        c_u * p_u.x + c_v * p_v.x,
        c_u * p_u.y + c_v * p_v.y,
        c_u * p_u.z + c_v * p_v.z,
    )
}
