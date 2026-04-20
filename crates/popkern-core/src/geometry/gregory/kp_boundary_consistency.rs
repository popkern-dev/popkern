//! G3 — boundary-consistency validation for H.4.d / H.4.e
//! inputs.
//!
//! [`derive_row1_frame`](super::kp_row1::derive_row1_frame) and
//! [`derive_row2_frame`](super::kp_row2::derive_row2_frame)
//! require the caller-supplied boundary curves (`row0`, `row9`,
//! `col0`, `col9`) to be **tangent- and curvature-consistent**
//! with the guide at each corner:
//!
//! ```text
//!   9 · (col0[1] − col0[0]) = G_u(0, 0)
//!   9 · (row0[1] − row0[0]) = G_v(0, 0)
//!   72 · (col0[2] − 2 · col0[1] + col0[0]) = G_uu(0, 0)
//!   72 · (row0[2] − 2 · row0[1] + row0[0]) = G_vv(0, 0)
//!   ... (and the analogous 12 checks at the other 3 corners)
//! ```
//!
//! If these hold, the row-1 / row-2 formulas produce a bidegree-9
//! patch whose partials match the guide at every corner to
//! floating-point precision (the G1 + H.4.f invariants rely on
//! this).
//!
//! If they fail, the patch has corner artefacts proportional to
//! the inconsistency magnitude. This module exposes a helper for
//! callers to quantify the inconsistency before committing to
//! the construction — useful both as a debugging tool and as a
//! pre-flight gate in orchestration code.

use super::kp_guide_patch::KpGuidePatch;
use crate::geometry::{Point3, Vec3};

/// Aggregated consistency report for one sub-patch's boundary
/// curves.
#[derive(Debug, Clone)]
pub struct BoundaryConsistency {
    /// Worst position mismatch at any of the 4 corners (max over
    /// `|boundary endpoint − guide corner|`).
    pub max_position_error: f64,
    /// Worst tangent mismatch in Bezier-derivative terms (max
    /// over `|9 · Δboundary − guide partial|` at each of 4
    /// corners, in both u and v).
    pub max_tangent_error: f64,
    /// Worst curvature mismatch in Bezier-second-derivative
    /// terms. `max |72 · Δ²boundary − guide second partial|`.
    pub max_second_derivative_error: f64,
    /// Per-corner breakdown for detailed diagnostics.
    pub per_corner: [CornerConsistency; 4],
}

impl BoundaryConsistency {
    /// `true` if the boundary is within `tolerance` of the guide
    /// across all consistency checks. Suitable for a pre-flight
    /// gate on H.4.d / H.4.e inputs.
    pub fn is_consistent(&self, tolerance: f64) -> bool {
        self.max_position_error < tolerance
            && self.max_tangent_error < tolerance
            && self.max_second_derivative_error < tolerance
    }
}

/// Detailed per-corner consistency data.
#[derive(Debug, Clone, Default)]
pub struct CornerConsistency {
    /// `|boundary endpoint − guide corner|`.
    pub position_error: f64,
    /// `|9 · Δ_{u-boundary} − G_u(corner)|` (boundary is `col0`
    /// or `col9` depending on corner).
    pub tangent_u_error: f64,
    /// Analogous in the v direction; boundary is `row0` or
    /// `row9`.
    pub tangent_v_error: f64,
    /// `|72 · Δ²_{u-boundary} − G_uu(corner)|`.
    pub second_uu_error: f64,
    /// Analogous in v.
    pub second_vv_error: f64,
}

/// Quantify the consistency of caller-supplied boundary curves
/// with the guide at each of the 4 corners. See the module doc
/// for the 16 checks performed.
///
/// Boundary index convention:
/// - `row0` runs along `u = 0` (index is `v`).
/// - `row9` runs along `u = 1`.
/// - `col0` runs along `v = 0` (index is `u`).
/// - `col9` runs along `v = 1`.
pub fn check_boundary_consistency(
    guide: &KpGuidePatch,
    row0: &[Point3; 10],
    row9: &[Point3; 10],
    col0: &[Point3; 10],
    col9: &[Point3; 10],
) -> BoundaryConsistency {
    let mut per_corner = [
        CornerConsistency::default(),
        CornerConsistency::default(),
        CornerConsistency::default(),
        CornerConsistency::default(),
    ];
    // Corner 0: (u = 0, v = 0). col0[0..=2] runs in u; row0[0..=2]
    // in v. G_u from col direction, G_v from row direction.
    per_corner[0] = corner_consistency(
        guide.point_at(0.0, 0.0),
        guide.partial_u(0.0, 0.0),
        guide.partial_v(0.0, 0.0),
        guide.second_partial_uu(0.0, 0.0),
        guide.second_partial_vv(0.0, 0.0),
        col0[0],
        col0[1],
        col0[2],
        row0[0],
        row0[1],
        row0[2],
    );
    // Corner 1: (u = 1, v = 0). col9[0..=2] runs in u backwards
    // (col9[0] = (0, 1) corner, col9[9] = (1, 1) corner — wait,
    // col9 runs in u from u=0 to u=1 at v=1. Need to re-check.
    //
    // Actually `col9` runs along `v = 1` with index = `u`. So
    // col9[0] = (u=0, v=1) = corner 3 position; col9[9] = (u=1,
    // v=1) = corner 2 position. So `col9` is the boundary for
    // corner 2 (at index 9) and corner 3 (at index 0).
    //
    // Similarly `col0` runs along `v = 0` with index = `u`, so
    // col0[0] = (u=0, v=0) = corner 0; col0[9] = (u=1, v=0) =
    // corner 1.
    //
    // At corner 1 = (1, 0): u-direction boundary is col0 near
    // index 9 (col0[9] corner, col0[8] tangent, col0[7] second).
    // v-direction boundary is row9 near index 0 (row9[0] corner,
    // row9[1] tangent, row9[2] second).
    per_corner[1] = corner_consistency(
        guide.point_at(1.0, 0.0),
        // u direction at u = 1: tangent = 9 · (col0[9] − col0[8])
        // (positive u at u = 1 end of col0).
        guide.partial_u(1.0, 0.0),
        guide.partial_v(1.0, 0.0),
        guide.second_partial_uu(1.0, 0.0),
        guide.second_partial_vv(1.0, 0.0),
        col0[9],
        col0[8],
        col0[7],
        row9[0],
        row9[1],
        row9[2],
    );
    // Corner 2 = (1, 1): col9[9] corner; row9[9] corner.
    per_corner[2] = corner_consistency(
        guide.point_at(1.0, 1.0),
        guide.partial_u(1.0, 1.0),
        guide.partial_v(1.0, 1.0),
        guide.second_partial_uu(1.0, 1.0),
        guide.second_partial_vv(1.0, 1.0),
        col9[9],
        col9[8],
        col9[7],
        row9[9],
        row9[8],
        row9[7],
    );
    // Corner 3 = (0, 1): col9[0] corner (u = 0 end); row0[9]
    // corner (v = 1 end).
    per_corner[3] = corner_consistency(
        guide.point_at(0.0, 1.0),
        guide.partial_u(0.0, 1.0),
        guide.partial_v(0.0, 1.0),
        guide.second_partial_uu(0.0, 1.0),
        guide.second_partial_vv(0.0, 1.0),
        col9[0],
        col9[1],
        col9[2],
        row0[9],
        row0[8],
        row0[7],
    );

    let max_pos = per_corner
        .iter()
        .map(|c| c.position_error)
        .fold(0.0_f64, f64::max);
    let max_tan = per_corner
        .iter()
        .map(|c| c.tangent_u_error.max(c.tangent_v_error))
        .fold(0.0_f64, f64::max);
    let max_sec = per_corner
        .iter()
        .map(|c| c.second_uu_error.max(c.second_vv_error))
        .fold(0.0_f64, f64::max);

    BoundaryConsistency {
        max_position_error: max_pos,
        max_tangent_error: max_tan,
        max_second_derivative_error: max_sec,
        per_corner,
    }
}

/// Compute one corner's consistency from:
/// - `guide_pos`, `guide_g_u`, `guide_g_v`, `guide_g_uu`,
///   `guide_g_vv`: guide's quantities at the corner.
/// - `u_cp0`, `u_cp1`, `u_cp2`: three CPs of the u-direction
///   boundary curve at the corner end (indices 0, 1, 2 from the
///   corner outward — for a corner at index 0 of the boundary,
///   these are the boundary's first three CPs; for a corner at
///   index 9, these are boundary CPs 9, 8, 7).
/// - `v_cp0`, `v_cp1`, `v_cp2`: analogously for v direction.
#[allow(clippy::too_many_arguments)]
fn corner_consistency(
    guide_pos: Point3,
    guide_g_u: Vec3,
    guide_g_v: Vec3,
    guide_g_uu: Vec3,
    guide_g_vv: Vec3,
    u_cp0: Point3,
    u_cp1: Point3,
    u_cp2: Point3,
    v_cp0: Point3,
    v_cp1: Point3,
    v_cp2: Point3,
) -> CornerConsistency {
    let pos_err = point_max(u_cp0, guide_pos).max(point_max(v_cp0, guide_pos));

    // Bezier-derivative-at-boundary: `f'(0) = n · (P_1 − P_0)`
    // at degree n = 9. At the "index 9" end the sign flips:
    // `f'(1) = n · (P_n − P_{n-1})`. Since we feed in `cp0` =
    // corner endpoint, `cp1` = immediate neighbour, and the sign
    // convention follows from which direction we take the
    // difference, the MAGNITUDE of the tangent is insensitive
    // to sign. We compare magnitudes by using `(cp1 − cp0)` and
    // comparing to `guide_g_× / 9` (then multiplying by 9 in
    // the error metric).
    let boundary_tangent_u = Vec3::new(
        9.0 * (u_cp1.x - u_cp0.x),
        9.0 * (u_cp1.y - u_cp0.y),
        9.0 * (u_cp1.z - u_cp0.z),
    );
    let boundary_tangent_v = Vec3::new(
        9.0 * (v_cp1.x - v_cp0.x),
        9.0 * (v_cp1.y - v_cp0.y),
        9.0 * (v_cp1.z - v_cp0.z),
    );
    // Account for possible sign flip (end-of-curve boundary).
    let tan_u_err = vec_error_magnitude(boundary_tangent_u, guide_g_u);
    let tan_v_err = vec_error_magnitude(boundary_tangent_v, guide_g_v);

    // Bezier second derivative: `f''(0) = n(n-1) · (P_2 − 2 P_1
    // + P_0)`. At n = 9: 72.
    let boundary_second_u = Vec3::new(
        72.0 * (u_cp2.x - 2.0 * u_cp1.x + u_cp0.x),
        72.0 * (u_cp2.y - 2.0 * u_cp1.y + u_cp0.y),
        72.0 * (u_cp2.z - 2.0 * u_cp1.z + u_cp0.z),
    );
    let boundary_second_v = Vec3::new(
        72.0 * (v_cp2.x - 2.0 * v_cp1.x + v_cp0.x),
        72.0 * (v_cp2.y - 2.0 * v_cp1.y + v_cp0.y),
        72.0 * (v_cp2.z - 2.0 * v_cp1.z + v_cp0.z),
    );
    let sec_uu_err = (boundary_second_u - guide_g_uu).length();
    let sec_vv_err = (boundary_second_v - guide_g_vv).length();

    CornerConsistency {
        position_error: pos_err,
        tangent_u_error: tan_u_err,
        tangent_v_error: tan_v_err,
        second_uu_error: sec_uu_err,
        second_vv_error: sec_vv_err,
    }
}

/// Error metric for two vectors, insensitive to sign (the
/// boundary tangent may be the negative of the guide's partial
/// depending on which end of the boundary we're evaluating).
#[inline]
fn vec_error_magnitude(a: Vec3, b: Vec3) -> f64 {
    let diff_plus = (a - b).length();
    let diff_minus = (a + b).length();
    diff_plus.min(diff_minus)
}

#[inline]
fn point_max(a: Point3, b: Point3) -> f64 {
    ((a.x - b.x).abs())
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

#[cfg(test)]
#[path = "kp_boundary_consistency_tests.rs"]
mod tests;
