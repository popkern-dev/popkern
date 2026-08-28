//! Half-space split of a closed solid by a plane (the planar arm of §2.3).
//!
//! [`split_solid_by_plane`] produces up to two closed sub-solids from a
//! single input solid and an infinite cutting plane, using two opposing
//! cutting boxes and the existing boolean subtraction path. Curved
//! sheets require the connected-component approach in
//! [`crate::sheet::split_curved`].

use popkern_core::topology::BRep;

/// Compute the world-space AABB of a BRep from its vertices.
pub(super) fn brep_aabb(
    brep: &BRep,
) -> Option<(
    popkern_core::geometry::Point3,
    popkern_core::geometry::Point3,
)> {
    let mut iter = brep.vertices.values();
    let first = iter.next()?.position;
    let mut min = first;
    let mut max = first;
    for v in iter {
        let p = v.position;
        if p.x < min.x {
            min.x = p.x;
        }
        if p.y < min.y {
            min.y = p.y;
        }
        if p.z < min.z {
            min.z = p.z;
        }
        if p.x > max.x {
            max.x = p.x;
        }
        if p.y > max.y {
            max.y = p.y;
        }
        if p.z > max.z {
            max.z = p.z;
        }
    }
    Some((min, max))
}

/// Split a solid by an infinite plane into closed half-solids.
///
/// Constructs two properly-sized cutting boxes (one on each side of
/// the plane) and subtracts each from the input solid in turn,
/// producing the opposite half. Returns the pieces; returns a single
/// clone of the input when the plane misses the body.
///
/// This is the building block for multi-piece `SolidMinusSheet` when
/// the sheet is a single planar face. Unlike
/// [`crate::ops::split_slice::split_by_plane`], the cutting box is
/// sized to the body's actual extent so the intersection produces
/// correct results for arbitrarily-sized solids.
pub fn split_solid_by_plane(
    solid: &BRep,
    plane_origin: popkern_core::geometry::Point3,
    plane_normal: popkern_core::geometry::Vec3,
    tol: popkern_core::tolerance::Tolerance,
) -> Result<Vec<BRep>, crate::BoolError> {
    use popkern_core::geometry::{Point3, Transform3, Vec3};

    let normal = plane_normal
        .normalized()
        .unwrap_or(Vec3::new(0.0, 0.0, 1.0));

    let (mn, mx) = brep_aabb(solid)
        .ok_or_else(|| crate::BoolError::DegenerateInput("solid has no vertices".to_string()))?;
    let dx = mx.x - mn.x;
    let dy = mx.y - mn.y;
    let dz = mx.z - mn.z;
    let max_extent = dx.max(dy).max(dz);

    let body_centroid = Point3::new(
        0.5 * (mn.x + mx.x),
        0.5 * (mn.y + mx.y),
        0.5 * (mn.z + mx.z),
    );

    // Early-out: if the body is entirely on one side of the plane, no
    // cut happens — return the body unchanged.
    let mut min_signed = f64::INFINITY;
    let mut max_signed = f64::NEG_INFINITY;
    for v in solid.vertices.values() {
        let d = (v.position - plane_origin).dot(&normal);
        if d < min_signed {
            min_signed = d;
        }
        if d > max_signed {
            max_signed = d;
        }
    }
    let on_plane_tol = tol.linear() * 10.0;
    if min_signed > -on_plane_tol || max_signed < on_plane_tol {
        return Ok(vec![solid.clone()]);
    }

    // Keep the cutting box only ~10% larger than the body. Massively
    // oversized boxes trip a shortcut bug in the boolean engine where
    // a point-sampling no-intersection check selects a vertex that
    // lands outside the true overlap region and concludes the bodies
    // are disjoint. A small box forces the full pipeline path.
    let margin = max_extent * 0.05 + tol.linear() * 10.0;
    let box_size = max_extent + 2.0 * margin;

    let unit = BRep::unit_cube(tol).map_err(crate::BoolError::CoreError)?;

    // Build a half-space cutting box on the side of the plane indicated
    // by `side_normal`. Anchor the box center at centroid + side_normal
    // * (half_box - signed_dist_to_plane) so one edge lies on the plane
    // and the box extends `box_size` units along `side_normal`.
    let make_half_box = |side_normal: Vec3| -> Result<BRep, crate::BoolError> {
        let centroid_offset = body_centroid - plane_origin;
        let signed_dist_to_plane = centroid_offset.dot(&side_normal);

        let half_box = box_size * 0.5;
        let box_center = body_centroid + side_normal * (half_box - signed_dist_to_plane);

        let scale_t = Transform3::scale_uniform(box_size);
        let scaled_center = box_size * 0.5;
        let translate_t = Transform3::translation(
            box_center.x - scaled_center,
            box_center.y - scaled_center,
            box_center.z - scaled_center,
        );
        let combined = translate_t.compose(&scale_t);
        unit.clone()
            .transform(&combined)
            .map_err(crate::BoolError::CoreError)
    };

    let pos_box = make_half_box(normal)?;
    let neg_box = make_half_box(normal * -1.0)?;

    // Produce each half by SUBTRACTING the OPPOSITE half-box from the
    // solid. `solid - neg_box` leaves everything not covered by the
    // negative half-space — i.e. the positive half. We use subtraction
    // because `intersection_ref` hits a shortcut bug for half-space
    // intersections, while subtraction handles cube-minus-box partial
    // overlap correctly.
    let mut pieces = Vec::new();
    if let Ok(positive) = crate::subtraction_ref(solid, &neg_box, tol, None)
        && !positive.faces.is_empty()
    {
        pieces.push(positive);
    }
    if let Ok(negative) = crate::subtraction_ref(solid, &pos_box, tol, None)
        && !negative.faces.is_empty()
    {
        pieces.push(negative);
    }
    Ok(pieces)
}
