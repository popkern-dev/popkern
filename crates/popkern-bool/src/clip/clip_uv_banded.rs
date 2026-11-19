//! Containment decided away from a face's uv boundary, against its
//! boundary sampled densely (split out of `clip_uv.rs` for the 500-line
//! ratchet).

use popkern_core::topology::{BRep, FaceHandle};

/// Containment of `q` decided away from the face's uv boundary: `None`
/// when `q` lies within `band` of the boundary (an on-boundary verdict
/// is a coin flip), else the dense inside test.
pub(super) fn banded_inside(
    brep: &BRep,
    face: FaceHandle,
    q: popkern_core::geometry::Point3,
    band: f64,
) -> Option<bool> {
    let boundary_ring = crate::clip::face_uv_boundary(brep, face);
    let mut best = f64::INFINITY;
    for (_, pc) in boundary_ring.iter().flatten() {
        let (a, b) = pc.parameter_range();
        if !(a.is_finite() && b.is_finite()) {
            continue;
        }
        for i in 0..=64 {
            let t = a + (b - a) * i as f64 / 64.0;
            let s = pc.point_at(t);
            best = best.min(((s.x - q.x).powi(2) + (s.y - q.y).powi(2)).sqrt());
        }
    }
    if best < band {
        return None;
    }
    Some(
        crate::split_inner_loop::point_inside_face_region_uv_with_density(
            brep,
            face,
            q,
            &brep.tolerance,
            64,
        ),
    )
}
