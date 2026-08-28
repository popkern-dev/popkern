//! Extract Hermite data for the hex-edge corners of
//! every sub-patch's guide surface.
//!
//! Each sub-patch `k` has three non-central corners on the hex
//! boundary:
//!
//! * `c0 = (u = 0, v = 0)` at `V_{2k}` (= `outer_corners[(k + N
//!   − 1) % N]`).
//! * `c2 = (u = 1, v = 1)` at `V_{2k + 2}` (= `outer_corners[k]`).
//! * `c3 = (u = 0, v = 1)` at `M_k` (the hex-edge midpoint).
//!
//! This module computes all three, filling in the hex-edge
//! tangents that H.4.c.3 left as placeholders.
//!
//! # M_k — tangent-plane-consistent hex midpoint (G5 + G6)
//!
//! The hex edge between `V_{2k}` and `V_{2k + 2}` is nominally a
//! curve on the input face surface. We approximate it with a
//! **cubic Bezier** whose endpoint tangents lie in the respective
//! input face tangent planes at `V_{2k}` and `V_{2k + 2}`, with
//! magnitudes set so `M_k` is the midpoint of the Bezier at
//! `t = 1/2`. The tangent directions are perpendicular to the
//! corresponding spoke tangents (projected onto the face plane),
//! giving the "natural" hex-edge departure direction at each V.
//!
//! Concretely, at `V_{2k}`: the spoke tangent is the direction
//! of spoke `k − 1` at its `t = 1` end; the hex-edge departure
//! direction is a unit vector in the face plane at `V_{2k}`,
//! perpendicular to that spoke tangent, oriented so the hex
//! edge heads roughly toward `V_{2k + 2}`. Same at `V_{2k + 2}`
//! with spoke `k`.
//!
//! For the **face tangent plane at M_k** (needed to project the
//! hex-edge tangents at `c3`), we blend the two adjacent outer-
//! vertex face tangent planes' normals and normalise. For planar
//! inputs this gives the same plane back; for non-planar it's a
//! reasonable interpolation without requiring an API extension
//! to [`DirectBoundaryInput`].
//!
//! # Straight-line hex edges
//!
//! With `M_k` at the midpoint, the sub-patch's left edge (`u =
//! 0`, `V_{2k} → M_k`) and top edge (`v = 1`, `M_k → V_{2k+2}`)
//! are both straight lines. For a cubic Bezier curve from `A` to
//! `B` with uniformly-spaced CPs (the straight-line case), the
//! tangent at both endpoints equals `B − A`:
//!
//! ```text
//!   tangent along left edge  = M_k − V_{2k}        (constant)
//!   tangent along top edge   = V_{2k + 2} − M_k    (constant)
//! ```
//!
//! Hex-edge tangents at the outer vertices (`c0.tangent_v` at
//! `V_{2k}` and `c2.tangent_u` at `V_{2k+2}`) are projected onto
//! the respective input-face tangent plane to enforce the G1
//! outer-corner constraint. At `M_k` (the `c3` corner) we don't
//! have an input plane, so tangents are used raw.
//!
//! # Along-spoke tangents
//!
//! `c0.tangent_u` (along-spoke at `V_{2k}`) and `c2.tangent_v`
//! (along-spoke at `V_{2k+2}`) are reused from
//! [`extract_outer_corner_hermite`]. This keeps the two
//! extractors composable without recomputing spoke tangents.
//!
//! # Twist
//!
//! Zero everywhere in this phase; proper twists are resolved in
//! H.4.c.5 ring-wide assembly.

use super::kp_guide_build::GuideCornerData;
use super::kp_guide_outer_hermite::extract_outer_corner_hermite;
use crate::geometry::gregory::direct_input::DirectBoundaryInput;
use crate::geometry::{Point3, Vec3};

/// Hermite data for the three non-central corners (`c0`, `c2`,
/// `c3`) of one sub-patch's guide.
///
/// The central-vertex corner `c1` is **not** included here; that
/// comes from
/// [`extract_central_vertex_hermite`](super::extract_central_vertex_hermite).
/// Ring-wide assembly in H.4.c.5 combines both.
#[derive(Debug, Clone)]
pub struct SubPatchHexCorners {
    /// `c0 = (u = 0, v = 0)` at `V_{2k}`.
    pub c0: GuideCornerData,
    /// `c2 = (u = 1, v = 1)` at `V_{2k + 2}`.
    pub c2: GuideCornerData,
    /// `c3 = (u = 0, v = 1)` at `M_k` (hex-edge midpoint).
    pub c3: GuideCornerData,
    /// Hex midpoint position `M_k`, exposed for higher-level
    /// callers (e.g. visualisation, validation).
    pub hex_midpoint: Point3,
}

/// Produce [`SubPatchHexCorners`] for every sub-patch in the
/// ring. Length of the returned `Vec` equals `input.n()`.
pub fn extract_hex_corners(input: &DirectBoundaryInput) -> Vec<SubPatchHexCorners> {
    let n = input.n();
    let outer = extract_outer_corner_hermite(input);
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let k_prev = (k + n - 1) % n;
        let v_prev = input.outer_corners[k_prev].position; // V_{2k}
        let v_next = input.outer_corners[k].position; // V_{2k + 2}
        let n_prev = input.outer_corners[k_prev].tangent_plane.normal();
        let n_next = input.outer_corners[k].tangent_plane.normal();

        // G5: tangent-plane-consistent hex midpoint. Build a
        // cubic Bezier from V_prev to V_next with endpoint
        // tangents perpendicular to the respective spoke tangents
        // (in each face plane), magnitudes chosen so the Bezier
        // evaluates to the chord midpoint at t = 1/2 for a
        // reasonable default shape. `M_k` is the Bezier at
        // t = 1/2.
        let chord = Vec3::new(
            v_next.x - v_prev.x,
            v_next.y - v_prev.y,
            v_next.z - v_prev.z,
        );
        let hex_tangent_at_v_prev =
            hex_edge_departure(outer[k_prev].patch_kp1_c0.tangent_u, chord, n_prev);
        let hex_tangent_at_v_next = hex_edge_departure(
            outer[k].patch_k_c2.tangent_v * -1.0,
            Vec3::new(-chord.x, -chord.y, -chord.z),
            n_next,
        );
        let m_k = hex_edge_midpoint(v_prev, hex_tangent_at_v_prev, hex_tangent_at_v_next, v_next);

        // G6: derive a face tangent plane at M_k by blending the
        // two adjacent V tangent-plane normals. For planar inputs
        // this reduces to the same plane; for non-planar inputs
        // it's a geometrically reasonable interpolation.
        let n_m_k = blend_normals(n_prev, n_next);

        // Left-edge (u = 0) and top-edge (v = 1) tangents — at
        // M_k we now project onto the blended plane n_m_k rather
        // than using raw straight-line vectors.
        let tangent_left_at_v_prev =
            Vec3::new(m_k.x - v_prev.x, m_k.y - v_prev.y, m_k.z - v_prev.z);
        let tangent_top_at_v_next = Vec3::new(v_next.x - m_k.x, v_next.y - m_k.y, v_next.z - m_k.z);

        let left_proj_at_v_prev = project_onto_plane(tangent_left_at_v_prev, n_prev);
        let top_proj_at_v_next = project_onto_plane(tangent_top_at_v_next, n_next);
        let left_tangent_at_m = project_onto_plane(tangent_left_at_v_prev, n_m_k);
        let top_tangent_at_m = project_onto_plane(tangent_top_at_v_next, n_m_k);

        let spoke_at_c0 = outer[k_prev].patch_kp1_c0.tangent_u;
        let spoke_at_c2 = outer[k].patch_k_c2.tangent_v;

        let c0 = GuideCornerData {
            position: v_prev,
            tangent_u: spoke_at_c0,
            tangent_v: left_proj_at_v_prev,
            twist: Vec3::zero(),
        };
        let c2 = GuideCornerData {
            position: v_next,
            tangent_u: top_proj_at_v_next,
            tangent_v: spoke_at_c2,
            twist: Vec3::zero(),
        };
        let c3 = GuideCornerData {
            position: m_k,
            tangent_u: top_tangent_at_m,
            tangent_v: left_tangent_at_m,
            twist: Vec3::zero(),
        };

        out.push(SubPatchHexCorners {
            c0,
            c2,
            c3,
            hex_midpoint: m_k,
        });
    }
    out
}

/// Hex-edge departure direction at an outer vertex: project the
/// chord `v_next − v_prev` (or its inverse at the far end) onto
/// the face tangent plane, then remove the component along the
/// spoke tangent. Length scaled to the chord's in-plane length
/// so the cubic Hermite midpoint is close to the chord midpoint.
///
/// Falls back to the chord projection if the spoke tangent and
/// chord are parallel in the face plane (degenerate).
fn hex_edge_departure(spoke_tangent: Vec3, chord: Vec3, face_normal: Vec3) -> Vec3 {
    let chord_in_plane = project_onto_plane(chord, face_normal);
    // Remove the component along the spoke (in face plane).
    let Some(spoke_in_plane) = project_onto_plane(spoke_tangent, face_normal).normalized() else {
        return chord_in_plane;
    };
    let dot = chord_in_plane.x * spoke_in_plane.x
        + chord_in_plane.y * spoke_in_plane.y
        + chord_in_plane.z * spoke_in_plane.z;
    Vec3::new(
        chord_in_plane.x - dot * spoke_in_plane.x,
        chord_in_plane.y - dot * spoke_in_plane.y,
        chord_in_plane.z - dot * spoke_in_plane.z,
    )
}

/// Cubic Hermite midpoint given endpoints + endpoint tangents.
/// For a cubic Bezier `B(t)` with
/// `P0 = a`, `P1 = a + T_a / 3`, `P2 = b − T_b / 3`, `P3 = b`
/// the midpoint is `B(1/2) = (a + b)/2 + (T_a − T_b)/8`.
fn hex_edge_midpoint(a: Point3, t_a: Vec3, t_b: Vec3, b: Point3) -> Point3 {
    Point3::new(
        0.5 * (a.x + b.x) + (t_a.x - t_b.x) / 8.0,
        0.5 * (a.y + b.y) + (t_a.y - t_b.y) / 8.0,
        0.5 * (a.z + b.z) + (t_a.z - t_b.z) / 8.0,
    )
}

/// Average + re-normalise two unit normals. Returns the first
/// normal unchanged when the inputs cancel (should not happen
/// for well-posed inputs).
fn blend_normals(a: Vec3, b: Vec3) -> Vec3 {
    // Ensure the two normals are on the same side before summing
    // — flip `b` if it points opposite to `a`.
    let dot = a.x * b.x + a.y * b.y + a.z * b.z;
    let b_aligned = if dot < 0.0 {
        Vec3::new(-b.x, -b.y, -b.z)
    } else {
        b
    };
    let sum = Vec3::new(a.x + b_aligned.x, a.y + b_aligned.y, a.z + b_aligned.z);
    sum.normalized().unwrap_or(a)
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
#[path = "kp_guide_hex_hermite_tests.rs"]
mod tests;
