//! Where a declined result's contact may lie (design §6.2).
//!
//! A touch is located only as well as its tangency conditions it: the
//! operands' gap grows with the square of the distance from the touch,
//! so it stays under the tolerance ε over a patch √(2εL) across (L the
//! scene's scale, its operands' curvature radius), and the kernel's
//! vertex identity may place the contact's vertex anywhere it lies on
//! both operands there (a tangent crossing 1e-5 off a sphere's seam
//! joins the seam's crossing vertex, generated corpus SW-V scene 836).
//! A contact is on the scene's touch when it lies within the fixed reach
//! of it, or inside that patch on both operands' surfaces within a
//! vertex's reach (twice the linear tolerance), measured exactly on the
//! faces, inside their trims.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval as _;
use popkern_core::topology::BRep;
use popkern_core::validation::FaceUvContainment;

/// Distance from `q` to `body`'s faces, read on each face's surface at
/// the closest point inside its trim; `None` when no face's closest
/// point lies inside it.
pub fn boundary_distance(body: &BRep, q: Point3) -> Option<f64> {
    body.faces
        .iter()
        .filter_map(|(fh, f)| {
            let s = f.surface.as_ref()?;
            let (u, v) = s.geometry.closest_point(&q)?;
            match body.face_uv_trim_containment(fh, &s.geometry, u, v) {
                None | Some(FaceUvContainment::Outside) => None,
                Some(_) => Some(s.geometry.point_at(u, v).distance_to(&q)),
            }
        })
        .min_by(f64::total_cmp)
}

/// Whether contact `q` lies on the touch at `contact` (see the module
/// doc): `reach` the fixed reach, `linear` the tolerance, `scale` the
/// scene's.
pub fn on_touch(
    (a, b): (&BRep, &BRep),
    q: Point3,
    contact: Point3,
    (reach, linear, scale): (f64, f64, f64),
) -> bool {
    let d = q.distance_to(&contact);
    if d <= reach {
        return true;
    }
    let patch = (2.0 * linear * scale).sqrt();
    let on = |body: &BRep| boundary_distance(body, q).is_some_and(|g| g <= 2.0 * linear);
    d <= patch && on(a) && on(b)
}
