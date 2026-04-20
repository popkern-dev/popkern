//! Phase H — pairwise outer-corner twist (K-P 2004 §4.2 analog).
//!
//! At each outer hex vertex `V_{2k + 2}` two guides meet —
//! sub-patch `k` at `c2 = (u = 1, v = 1)` and sub-patch `(k + 1)
//! % N` at `c0 = (u = 0, v = 0)`. For G1 to hold along the full
//! shared spoke between these two guides (not just at the spoke
//! endpoint `V`), the **inner CPs adjacent to the shared edge at
//! `V`** must satisfy a pairwise Shirman-Séquin constraint:
//!
//! ```text
//!   (Y_k − Y_{k+1}) · p_spoke = 0
//! ```
//!
//! where `Y_k` is sub-patch `k`'s `grid[2][2]` inner CP, `Y_{k+1}`
//! is sub-patch `(k + 1)`'s `grid[1][1]` inner CP, and `p_spoke`
//! is the 2-D direction in the input face tangent plane at `V`
//! perpendicular to the shared spoke.
//!
//! # Closed-form solution
//!
//! With only **one** constraint between 2 unknowns the
//! Lagrangian QP degenerates to an explicit algebraic formula.
//! Let the parallelogram-rule seeds (in 2-D tangent-plane coords
//! at `V`) be `R_k`, `R_{k+1}` and the unit spoke-perpendicular
//! be `p_spoke`. Minimising
//! `‖Y_k − R_k‖² + ‖Y_{k+1} − R_{k+1}‖²` subject to the
//! constraint gives
//!
//! ```text
//!   λ = (R_k − R_{k+1}) · p_spoke
//!   Y_k       = R_k      − (λ / 2) · p_spoke
//!   Y_{k + 1} = R_{k + 1} + (λ / 2) · p_spoke
//! ```
//!
//! the projections onto `p_spoke` come out exactly equal by
//! construction, so the G1 constraint is satisfied and the
//! perturbation from the seeds is minimal.
//!
//! # Derived twists
//!
//! The bicubic Hermite inner-CP formulas at the two corners are:
//!
//! ```text
//!   sub-patch k at c2 (u = 1, v = 1):
//!     grid[2][2] = V − tu_k / 3 − tv_k / 3 + twist_k / 9
//!     ⟹ twist_k = 9 · (grid[2][2] − V) + 3 · tu_k + 3 · tv_k
//!
//!   sub-patch k + 1 at c0 (u = 0, v = 0):
//!     grid[1][1] = V + tu_{k + 1} / 3 + tv_{k + 1} / 3 + twist_{k + 1} / 9
//!     ⟹ twist_{k + 1} = 9 · (grid[1][1] − V) − 3 · tu_{k + 1} − 3 · tv_{k + 1}
//! ```

use super::direct_input::DirectBoundaryInput;
use super::direct_linalg::orthonormal_tangent_basis;
use super::kp_guide_build::GuideCornerData;
use super::kp_guide_hex_hermite::{SubPatchHexCorners, extract_hex_corners};
use crate::geometry::{Point3, Vec3};

/// Return [`SubPatchHexCorners`] for every sub-patch with
/// **properly pairwise-solved twists** at the `c0` and `c2`
/// outer corners. The `c3` (hex midpoint) twist is still zero
/// — H.4.c owns the hex-edge boundary and no G1-along-edge
/// constraint exists to pin its value.
///
/// Position and tangents are identical to the baseline
/// [`extract_hex_corners`]; only the `twist` field changes.
pub fn extract_hex_corners_with_twist(input: &DirectBoundaryInput) -> Vec<SubPatchHexCorners> {
    let base = extract_hex_corners(input);
    let n = base.len();
    let mut out = base.clone();

    for k in 0..n {
        let k_next = (k + 1) % n;
        // Shared outer vertex V_{2k + 2} = input.outer_corners[k].
        let v = input.outer_corners[k].position;
        let n_v = input.outer_corners[k].tangent_plane.normal();
        let (t1, t2) = orthonormal_tangent_basis(n_v);

        let c2_k = base[k].c2; // sub-patch k at c2 (u = 1, v = 1)
        let c0_kn = base[k_next].c0; // sub-patch (k + 1) at c0 (u = 0, v = 0)

        // Parallelogram seeds in 3-D:
        //   seed_k = V − tu_k/3 − tv_k/3   (for grid[2][2])
        //   seed_kn = V + tu_kn/3 + tv_kn/3  (for grid[1][1])
        let seed_k = Point3::new(
            v.x - c2_k.tangent_u.x / 3.0 - c2_k.tangent_v.x / 3.0,
            v.y - c2_k.tangent_u.y / 3.0 - c2_k.tangent_v.y / 3.0,
            v.z - c2_k.tangent_u.z / 3.0 - c2_k.tangent_v.z / 3.0,
        );
        let seed_kn = Point3::new(
            v.x + c0_kn.tangent_u.x / 3.0 + c0_kn.tangent_v.x / 3.0,
            v.y + c0_kn.tangent_u.y / 3.0 + c0_kn.tangent_v.y / 3.0,
            v.z + c0_kn.tangent_u.z / 3.0 + c0_kn.tangent_v.z / 3.0,
        );

        // Project seeds into the 2-D face tangent plane at V.
        let seed_k_2d = project_to_2d(v, seed_k, t1, t2);
        let seed_kn_2d = project_to_2d(v, seed_kn, t1, t2);

        // Spoke direction in 2-D: use c2.tangent_v (the along-
        // spoke tangent at V for sub-patch k). Equivalently,
        // minus c0_kn.tangent_u — both point toward / away from
        // `C`. Project into the (t1, t2) basis.
        let sigma = c2_k.tangent_v;
        let sx = sigma.dot(&t1);
        let sy = sigma.dot(&t2);
        let mag = (sx * sx + sy * sy).sqrt();
        if mag < 1e-12 {
            // Degenerate — spoke tangent vanishes in the face
            // tangent plane (e.g. spoke exactly perpendicular to
            // the face). Leave twists at their baseline (zero)
            // values and continue.
            continue;
        }
        let p_spoke = [-sy / mag, sx / mag]; // 90° CCW perp

        // Closed-form pairwise Lagrangian:
        //   λ = (R_k − R_{k+1}) · p_spoke
        //   Y_k = R_k − (λ / 2) · p_spoke
        //   Y_{k+1} = R_{k+1} + (λ / 2) · p_spoke
        let diff = [seed_k_2d[0] - seed_kn_2d[0], seed_k_2d[1] - seed_kn_2d[1]];
        let lambda = diff[0] * p_spoke[0] + diff[1] * p_spoke[1];
        let y_k_2d = [
            seed_k_2d[0] - (lambda / 2.0) * p_spoke[0],
            seed_k_2d[1] - (lambda / 2.0) * p_spoke[1],
        ];
        let y_kn_2d = [
            seed_kn_2d[0] + (lambda / 2.0) * p_spoke[0],
            seed_kn_2d[1] + (lambda / 2.0) * p_spoke[1],
        ];

        let y_k = lift_to_3d(v, y_k_2d, t1, t2);
        let y_kn = lift_to_3d(v, y_kn_2d, t1, t2);

        // Derive twists.
        let twist_k = Vec3::new(
            9.0 * (y_k.x - v.x) + 3.0 * c2_k.tangent_u.x + 3.0 * c2_k.tangent_v.x,
            9.0 * (y_k.y - v.y) + 3.0 * c2_k.tangent_u.y + 3.0 * c2_k.tangent_v.y,
            9.0 * (y_k.z - v.z) + 3.0 * c2_k.tangent_u.z + 3.0 * c2_k.tangent_v.z,
        );
        let twist_kn = Vec3::new(
            9.0 * (y_kn.x - v.x) - 3.0 * c0_kn.tangent_u.x - 3.0 * c0_kn.tangent_v.x,
            9.0 * (y_kn.y - v.y) - 3.0 * c0_kn.tangent_u.y - 3.0 * c0_kn.tangent_v.y,
            9.0 * (y_kn.z - v.z) - 3.0 * c0_kn.tangent_u.z - 3.0 * c0_kn.tangent_v.z,
        );

        out[k].c2 = GuideCornerData {
            twist: twist_k,
            ..c2_k
        };
        out[k_next].c0 = GuideCornerData {
            twist: twist_kn,
            ..c0_kn
        };
    }
    out
}

#[inline]
fn project_to_2d(origin: Point3, p: Point3, t1: Vec3, t2: Vec3) -> [f64; 2] {
    let off = Vec3::new(p.x - origin.x, p.y - origin.y, p.z - origin.z);
    [off.dot(&t1), off.dot(&t2)]
}

#[inline]
fn lift_to_3d(origin: Point3, xy: [f64; 2], t1: Vec3, t2: Vec3) -> Point3 {
    Point3::new(
        origin.x + xy[0] * t1.x + xy[1] * t2.x,
        origin.y + xy[0] * t1.y + xy[1] * t2.y,
        origin.z + xy[0] * t1.z + xy[1] * t2.z,
    )
}

#[cfg(test)]
#[path = "kp_guide_outer_twist_tests.rs"]
mod tests;
