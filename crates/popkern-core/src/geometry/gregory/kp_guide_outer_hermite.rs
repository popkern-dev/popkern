//! Extract the outer-hex-vertex Hermite data for
//! guides' `c0` and `c2` corners.
//!
//! Each outer hex vertex `V_{2k+2}` (= `input.outer_corners[k]`)
//! appears as
//!
//! * sub-patch `k`'s `c2` corner at `(u = 1, v = 1)`
//! * sub-patch `(k + 1) mod N`'s `c0` corner at `(u = 0, v = 0)`
//!
//! The two adjacent guides meeting at this vertex must share the
//! SAME tangent plane there for G1 to hold across the shared
//! spoke at its outer end. Per ADR §2.2 and K-P 2004 §4.2, the
//! input provides the target tangent plane via
//! [`OuterCornerData::tangent_plane`] (derived from the two
//! adjacent input faces' normals). Both sub-patches' Hermite
//! tangents at `V_{2k+2}` are projected into this plane.
//!
//! # Along-spoke tangents
//!
//! At `V_{2k+2}` the spoke `k` (from `C` outward) ends. Its
//! tangent at the `t = 1` (V) end is `7 · (cps[7] − cps[6])` for
//! the degree-7 biseptic input. In the guide geometry:
//!
//! * For sub-patch `k` at `c2 = (u = 1, v = 1)`: spoke `k` is
//!   the `u = 1` edge, parameterised forward. `∂G_k/∂v` at
//!   `(1, 1)` = `+ 7 · (V_{2k+2} − cps[6])` (tangent at the `v =
//!   1` end pointing **out** of the patch toward extensions
//!   beyond `V`).
//! * For sub-patch `(k + 1) % N` at `c0 = (u = 0, v = 0)`: spoke
//!   `k` is the `v = 0` edge, reversed. `∂G_{k+1}/∂u` at
//!   `(0, 0)` = `+ 7 · (cps[6] − V_{2k+2})` (tangent at `u = 0`
//!   pointing **into** the patch toward `C`).
//!
//! The two tangents are exact negatives of each other in raw
//! form, and remain negatives after projection onto a common
//! tangent plane — the shared-spoke parameter reversal is a
//! linear sign flip that commutes with orthogonal projection.
//!
//! # Hex-edge tangents deferred
//!
//! The *other* Hermite tangent at each of these corners (`∂G_k/∂u`
//! at `c2`, `∂G_{k+1}/∂v` at `c0`) describes the hex-external
//! boundary curve from `V_{2k+2}` toward the hex midpoint `M_k`
//! (patch `k`) or `M_{k+1}` (patch `k+1`). Those require the hex
//! midpoint positions / tangents, which come in **H.4.c.4**.
//! Here we fill them with [`Vec3::zero`] as placeholders.
//!
//! # Twist
//!
//! Set to [`Vec3::zero`] pending ring-wide assembly (H.4.c.5).

use super::kp_guide_build::GuideCornerData;
use crate::geometry::gregory::direct_input::DirectBoundaryInput;
use crate::geometry::{Point3, Vec3};

/// Hermite data for the pair of guide corners meeting at a single
/// outer hex vertex.
#[derive(Debug, Clone)]
pub struct OuterVertexHermite {
    /// Position of the outer vertex `V_{2k+2}`.
    pub outer_vertex: Point3,
    /// Index of the spoke whose `t = 1` end lands at this vertex.
    /// The two adjacent sub-patches are `spoke_index` and
    /// `(spoke_index + 1) % N`.
    pub spoke_index: usize,
    /// Hermite data for sub-patch `k = spoke_index` at its `c2`
    /// corner `(u = 1, v = 1)`. `tangent_v` is the along-spoke
    /// tangent; `tangent_u` is a placeholder (hex-edge, H.4.c.4).
    pub patch_k_c2: GuideCornerData,
    /// Hermite data for sub-patch `(spoke_index + 1) % N` at its
    /// `c0` corner `(u = 0, v = 0)`. `tangent_u` is the along-
    /// spoke tangent; `tangent_v` is a placeholder.
    pub patch_kp1_c0: GuideCornerData,
}

/// Extract outer-vertex Hermite pairs for every shared outer hex
/// vertex. Output length equals `input.n()`.
///
/// Each entry's along-spoke tangent is **projected** into the
/// input face tangent plane at that outer vertex (from
/// [`OuterCornerData::tangent_plane`](crate::geometry::gregory::direct_input::OuterCornerData::tangent_plane)),
/// guaranteeing the two adjacent guides have coplanar along-spoke
/// tangents and therefore share a tangent plane at `V_{2k+2}`
/// once the hex-edge tangents land in the same plane (H.4.c.4).
pub fn extract_outer_corner_hermite(input: &DirectBoundaryInput) -> Vec<OuterVertexHermite> {
    let n = input.n();
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let v = input.outer_corners[k].position;
        let n_v = input.outer_corners[k].tangent_plane.normal();
        let cps = &input.spokes[k].boundary_cps;

        // Along-spoke tangent at V for sub-patch k at c2:
        //   ∂G_k/∂v(1, 1) = 7 · (cps[7] − cps[6]) = 7 · (V − cps[6])
        let t_v_k_raw = Vec3::new(
            7.0 * (v.x - cps[6].x),
            7.0 * (v.y - cps[6].y),
            7.0 * (v.z - cps[6].z),
        );
        // Along-spoke tangent at V for sub-patch k+1 at c0:
        //   ∂G_{k+1}/∂u(0, 0) = 7 · (cps[6] − cps[7]) = − t_v_k_raw
        let t_u_kp1_raw = Vec3::new(-t_v_k_raw.x, -t_v_k_raw.y, -t_v_k_raw.z);

        let t_v_k = project_onto_plane(t_v_k_raw, n_v);
        let t_u_kp1 = project_onto_plane(t_u_kp1_raw, n_v);

        out.push(OuterVertexHermite {
            outer_vertex: v,
            spoke_index: k,
            patch_k_c2: GuideCornerData {
                position: v,
                tangent_u: Vec3::zero(),
                tangent_v: t_v_k,
                twist: Vec3::zero(),
            },
            patch_kp1_c0: GuideCornerData {
                position: v,
                tangent_u: t_u_kp1,
                tangent_v: Vec3::zero(),
                twist: Vec3::zero(),
            },
        });
    }
    out
}

/// Orthogonal projection of `v` onto the plane with unit
/// `normal`. `v' = v − (v · normal) · normal`.
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
#[path = "kp_guide_outer_hermite_tests.rs"]
mod tests;
