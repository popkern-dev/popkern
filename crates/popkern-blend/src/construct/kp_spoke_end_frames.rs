//! Spoke end frames for the K-P ring input (PK-547): each spoke
//! curve's end legs lie in the tangent plane the ring must share at
//! that end, so the ring is G1 across the spoke there by construction.
//!
//! At an outer corner `V_k`, the two patches sharing spoke `k` have
//! corner tangent planes `span(w, e_a)` and `span(w, e_b)`. Here `w` is
//! the spoke's end tangent, and `e_a`, `e_b` are the tangents of the two
//! hex edges meeting there. The planes agree only if `w` lies in
//! `span(e_a, e_b)`, the rim's tangent plane at `V_k`. At `t = 1` these
//! three tangents are fixed boundary data, so no interior control row
//! can repair a spoke whose end leg leaves the plane.
//!
//! Measured on cube_extreme 1:20: the seed spokes left the plane at
//! sines of 0.78, 0.69 and 0.34, and the served ring creased up to a
//! spoke normal sine of 0.325 over the last 5% of the spoke. The
//! promotion's sampled window `[0.05, 0.95]` could not see that.
//!
//! At the center `C`, every spoke's start leg must lie in the one
//! tangent plane all ring patches share there.

use popkern_core::geometry::gregory::{
    OuterCornerData, Side, Spoke, SpokeBoundaryData, TangentPlane,
};
use popkern_core::geometry::{Point3, Vec3};

/// Smallest sine between the two hex-edge tangents at an outer corner
/// for which they define the rim plane. Below it the rim runs smoothly
/// through the corner, both patch planes contain the common edge
/// direction, and any spoke leg gives G1 there.
const MIN_CORNER_SINE: f64 = 1.0e-3;

/// Unit normal of the rim tangent plane at an outer corner, spanned by
/// the two hex-edge tangents meeting there, or `None` when they are
/// (nearly) parallel.
pub(in crate::construct) fn corner_plane_normal(e_a: Vec3, e_b: Vec3) -> Option<Vec3> {
    let n = e_a.cross(&e_b);
    let scale = e_a.length() * e_b.length();
    if scale.is_nan() || scale <= 0.0 || n.length() < MIN_CORNER_SINE * scale {
        return None;
    }
    n.normalized()
}

/// Unit normal of the common center plane of a ring's spoke start
/// legs, given in cyclic ring order: Newell's average
/// `Σ d̂_k × d̂_{k+1}` over the unit legs. When the legs are coplanar it
/// is exactly their plane's normal, oriented by the ring's winding.
pub(in crate::construct) fn center_plane_normal(legs: &[Vec3]) -> Option<Vec3> {
    let units: Vec<Vec3> = legs.iter().filter_map(|d| d.normalized()).collect();
    if units.len() != legs.len() || units.len() < 3 {
        return None;
    }
    let n = units.len();
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    for k in 0..n {
        sum = sum + units[k].cross(&units[(k + 1) % n]);
    }
    sum.normalized()
}

/// The control point `handle` moved so that the leg `handle − anchor`
/// lies in the plane through `anchor` with unit normal `normal`,
/// keeping the leg's length (and so the curve's end speed). The leg is
/// projected onto the plane and rescaled. A leg along the normal has no
/// in-plane direction and is returned unchanged; the ring's gates then
/// judge it.
pub(in crate::construct) fn leg_into_plane(anchor: Point3, handle: Point3, normal: Vec3) -> Point3 {
    let leg = handle - anchor;
    let length = leg.length();
    let in_plane = leg - normal * leg.dot(&normal);
    let in_length = in_plane.length();
    if in_length.is_nan() || in_length <= 1.0e-9 * length {
        return handle;
    }
    anchor + in_plane * (length / in_length)
}

/// PK-547: rotate each spoke's end legs into the tangent planes the
/// ring must share there (`kp_spoke_end_frames`), and hand those
/// planes to the builders as the corner and center tangent planes.
/// The seed patches' own planes at `V` contain the misaligned spoke,
/// so they inherited the crease instead of preventing it. The rim plane
/// at `V_k` comes from the exact hex edges when the caller provides
/// them in the standard three-quad orientation, and otherwise from the
/// seed's cross tangents there, which are the adjacent boundaries'
/// tangents.
pub(in crate::construct) fn align_spoke_end_frames(
    spokes: &[Spoke],
    hex_boundaries: Option<&[([Point3; 10], [Point3; 10])]>,
    spoke_data: &mut [SpokeBoundaryData],
    outer_corners: &mut [OuterCornerData],
    central_vertex_tps: &mut [TangentPlane],
) {
    let oriented = |n: Vec3, like: Vec3| {
        if n.dot(&like) < 0.0 { n * -1.0 } else { n }
    };
    let legs: Vec<_> = spoke_data
        .iter()
        .map(|s| s.boundary_cps[1] - s.boundary_cps[0])
        .collect();
    let center_normal = center_plane_normal(&legs);
    for (k, spoke) in spokes.iter().enumerate() {
        let standard = spoke.a.side == Side::UOne
            && !spoke.a.reversed
            && spoke.b.side == Side::VZero
            && spoke.b.reversed;
        let data = &mut spoke_data[k];
        let (e_a, e_b) = match hex_boundaries {
            Some(hex) if standard => {
                let (a, b) = (&hex[spoke.a.patch_idx], &hex[spoke.b.patch_idx]);
                (a.1[9] - a.1[8], b.0[1] - b.0[0])
            }
            _ => (data.xtan_a_at_v, data.xtan_b_at_v),
        };
        if let Some(n_v) = corner_plane_normal(e_a, e_b) {
            data.boundary_cps[6] = leg_into_plane(data.boundary_cps[7], data.boundary_cps[6], n_v);
            let n_v = oriented(n_v, outer_corners[k].tangent_plane.normal());
            if let Some(plane) = TangentPlane::from_normal(n_v) {
                outer_corners[k].tangent_plane = plane;
            }
        }
        if let Some(n_c) = center_normal {
            data.boundary_cps[1] = leg_into_plane(data.boundary_cps[0], data.boundary_cps[1], n_c);
            let n_c = oriented(n_c, central_vertex_tps[k].normal());
            if let Some(plane) = TangentPlane::from_normal(n_c) {
                central_vertex_tps[k] = plane;
            }
        }
    }
}

#[cfg(test)]
#[path = "kp_spoke_end_frames_tests.rs"]
mod tests;
