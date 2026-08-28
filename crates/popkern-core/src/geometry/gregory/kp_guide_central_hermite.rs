//! Extract the C-corner Hermite data for every sub-
//! patch's guide surface.
//!
//! Given a [`DirectBoundaryInput`] describing an N-sided vertex
//! blend, produce per-sub-patch [`GuideCornerData`] at the
//! central-vertex corner `c1 = (u = 1, v = 0)` of each guide.
//!
//! # Guide corner layout recap
//!
//! Per ADR §2.2, guide patch `k` has
//!
//! ```text
//!   G_k(0, 0) = V_{2k}             (outer hex vertex, left)
//!   G_k(1, 0) = C                  (central vertex)
//!   G_k(1, 1) = V_{2k + 2}         (outer hex vertex, right)
//!   G_k(0, 1) = M_k                (hex midpoint)
//! ```
//!
//! so `C` is at the `c1` corner of every guide. This module
//! computes just that corner's Hermite data; the other three
//! corners are computed in H.4.c.3 (outer corner) and H.4.c.4
//! (hex midpoint), and the full 4-tuple is assembled in
//! H.4.c.5.
//!
//! # Tangent direction derivation
//!
//! At the `c1` corner of guide `k`:
//!
//! - `tangent_u` = `∂G_k/∂u` at `(u = 1, v = 0)`, which is the
//!   +u-direction tangent to guide's `v = 0` boundary at its
//!   `u = 1` end. The `v = 0` boundary runs from `V_{2k}` at
//!   `u = 0` to `C` at `u = 1`; in the input, this corresponds
//!   to spoke `k − 1` (reversed). Since spokes are parameterised
//!   with `t = 0` at `C`, the reversal flips the sign, and the
//!   result is `tangent_u = +7 · (C − spoke_{k−1}.cps[1])` at
//!   the biseptic input degree.
//!
//! - `tangent_v` = `∂G_k/∂v` at `(u = 1, v = 0)`, the +v-direction
//!   tangent at `C` along guide's `u = 1` boundary, which IS
//!   spoke `k` (forward). Straightforward: `tangent_v = +7 ·
//!   (spoke_k.cps[1] − C)`.
//!
//! # Shirman-Séquin N-way projection
//!
//! The raw spoke tangents do NOT in general lie in a common
//! tangent plane — each pair of adjacent spokes defines its own
//! plane at `C`. For every sub-patch's guide to have the **same**
//! tangent plane at `C` (required for G1 across shared spokes),
//! we project each sub-patch's two tangents onto the common
//! plane with normal `n_C`, where `n_C` is produced by
//! [`solve_central_vertex_corner`]. The projection is
//!
//! ```text
//!   v' = v − (v · n_C) · n_C
//! ```
//!
//! i.e. orthogonal projection that zeros the out-of-plane
//! component. After projection, every sub-patch's `tangent_u ∧
//! tangent_v` spans the plane orthogonal to `n_C` — so all N
//! tangent planes coincide by construction.
//!
//! # Twist
//!
//! The mixed-partial `∂²G/∂u∂v` at `C` ("twist") is set to
//! [`Vec3::zero`] for the bicubic guide in this phase. The
//! Shirman-Séquin classical construction uses a non-zero twist
//! at `C` (derived from the N-way inner-CP solve) to achieve
//! proper G1 at `C` for the WHOLE bicubic patch, not just the
//! tangent plane. Integrating that step requires cross-checking
//! the outer-corner solve results, so it is deferred to H.4.c.5
//! (ring-wide assembly).

use super::kp_guide_build::GuideCornerData;
use crate::geometry::gregory::direct_corner_c::solve_central_vertex_corner;
use crate::geometry::gregory::direct_errors::DirectConstructionError;
use crate::geometry::gregory::direct_input::DirectBoundaryInput;
use crate::geometry::gregory::direct_options::DirectConstructionOptions;
use crate::geometry::{Point3, Vec3};

/// Produce [`GuideCornerData`] for the `c1 = (u = 1, v = 0)`
/// corner of every sub-patch's guide, with tangents projected
/// into the common tangent plane at `C` (normal `n_C`).
///
/// Length of the returned `Vec` equals `input.n()`. Entry `k`
/// is the Hermite data for guide `k`'s `c1` corner.
///
/// Errors surface from the underlying [`solve_central_vertex_corner`] —
/// specifically if the N input face tangent planes are too
/// mutually inconsistent to agree on a common normal.
pub fn extract_central_vertex_hermite(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Result<Vec<GuideCornerData>, DirectConstructionError> {
    let central = solve_central_vertex_corner(input, options)?;
    let n_c = central.common_normal;
    let c = input.central_vertex;
    let n = input.n();

    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let k_prev = (k + n - 1) % n;

        let t_u_raw = tangent_u_raw(c, input.spokes[k_prev].boundary_cps[1]);
        let t_v_raw = tangent_v_raw(c, input.spokes[k].boundary_cps[1]);

        let t_u = project_onto_plane(t_u_raw, n_c);
        let t_v = project_onto_plane(t_v_raw, n_c);

        out.push(GuideCornerData {
            position: c,
            tangent_u: t_u,
            tangent_v: t_v,
            twist: Vec3::zero(),
        });
    }
    Ok(out)
}

/// Raw `tangent_u` at `C` for guide `k` — the degree-7 input
/// boundary tangent of spoke `k − 1` (reversed) at its `C` end.
///
/// `7 · (C − spoke[k-1].cps[1])` — the sign reflects that the
/// spoke is reversed when used as the guide's `v = 0` edge.
#[inline]
fn tangent_u_raw(c: Point3, spoke_prev_cp1: Point3) -> Vec3 {
    Vec3::new(
        7.0 * (c.x - spoke_prev_cp1.x),
        7.0 * (c.y - spoke_prev_cp1.y),
        7.0 * (c.z - spoke_prev_cp1.z),
    )
}

/// Raw `tangent_v` at `C` for guide `k` — the degree-7 input
/// boundary tangent of spoke `k` (forward) at its `C` end.
///
/// `7 · (spoke[k].cps[1] − C)`.
#[inline]
fn tangent_v_raw(c: Point3, spoke_cp1: Point3) -> Vec3 {
    Vec3::new(
        7.0 * (spoke_cp1.x - c.x),
        7.0 * (spoke_cp1.y - c.y),
        7.0 * (spoke_cp1.z - c.z),
    )
}

/// Orthogonal projection of `v` onto the plane with unit `normal`.
/// `v' = v − (v · normal) · normal`. Expects `normal` to be a
/// unit vector (callers upstream of [`solve_central_vertex_corner`]
/// guarantee this).
#[inline]
fn project_onto_plane(v: Vec3, normal: Vec3) -> Vec3 {
    let dot = v.x * normal.x + v.y * normal.y + v.z * normal.z;
    Vec3::new(
        v.x - dot * normal.x,
        v.y - dot * normal.y,
        v.z - dot * normal.z,
    )
}

#[cfg(test)]
#[path = "kp_guide_central_hermite_tests.rs"]
mod tests;
