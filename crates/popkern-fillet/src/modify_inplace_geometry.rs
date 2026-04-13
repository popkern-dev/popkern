//! Pure geometry helpers for [`crate::modify_inplace`]. Split out
//! to keep each file under the 500-line structural cap.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle};

use crate::error::BlendError;
use crate::unfillet_closed_helpers::{CylinderGeometry, PlaneGeometry};

/// Sign conventions of the existing fillet, derived from the
/// support geometries and the existing contact-vertex positions.
pub(crate) struct FilletSigns {
    /// +1 if the cyl-side vertex sits on the +n_p side of the
    /// plane, -1 on the -n_p side. Determines which direction
    /// along the plane normal the rolling-ball center sits.
    pub(crate) s_plane: f64,
    /// +1 if the plane-side vertex is outside the cylinder,
    /// -1 if inside. Determines whether the fillet wraps around
    /// outside or inside the cylinder.
    pub(crate) s_radial: f64,
}

pub(crate) fn determine_fillet_signs(
    plane: &PlaneGeometry,
    cyl: &CylinderGeometry,
    v_plane_pos: &Point3,
    v_cyl_pos: &Point3,
) -> Result<FilletSigns, BlendError> {
    // s_plane: along which direction of `cyl.axis` the cyl-side
    // vertex sits relative to the plane. Computed from the AXIS
    // projection — NOT the plane normal — so it gives the correct
    // sign for `move_cyl_side_vertex`'s `plane_along + s_plane * R`
    // formula regardless of how `plane.normal` is oriented.
    let denom = cyl.axis.dot(&plane.normal);
    let plane_along = if denom.abs() < 1e-12 {
        0.0
    } else {
        let diff = Vec3::new(
            plane.origin.x - cyl.origin.x,
            plane.origin.y - cyl.origin.y,
            plane.origin.z - cyl.origin.z,
        );
        diff.dot(&plane.normal) / denom
    };
    let cyl_along = (v_cyl_pos.x - cyl.origin.x) * cyl.axis.x
        + (v_cyl_pos.y - cyl.origin.y) * cyl.axis.y
        + (v_cyl_pos.z - cyl.origin.z) * cyl.axis.z;
    let cyl_offset_along_axis = cyl_along - plane_along;
    let s_plane = if cyl_offset_along_axis.abs() < 1e-12 {
        return Err(BlendError::InvalidParameter {
            parameter: "fillet_face".to_string(),
            value: "cyl-side vertex on plane".to_string(),
            reason: "degenerate fillet (zero radius?)".to_string(),
        });
    } else if cyl_offset_along_axis > 0.0 {
        1.0
    } else {
        -1.0
    };

    let v = Vec3::new(
        v_plane_pos.x - cyl.origin.x,
        v_plane_pos.y - cyl.origin.y,
        v_plane_pos.z - cyl.origin.z,
    );
    let along = v.dot(&cyl.axis);
    let perp = Vec3::new(
        v.x - along * cyl.axis.x,
        v.y - along * cyl.axis.y,
        v.z - along * cyl.axis.z,
    );
    let dist_from_axis = perp.length();
    let s_radial = if (dist_from_axis - cyl.radius).abs() < 1e-12 {
        return Err(BlendError::InvalidParameter {
            parameter: "fillet_face".to_string(),
            value: "plane-side vertex on cylinder".to_string(),
            reason: "degenerate fillet (zero radius?)".to_string(),
        });
    } else if dist_from_axis > cyl.radius {
        1.0
    } else {
        -1.0
    };

    Ok(FilletSigns { s_plane, s_radial })
}

/// Move the plane-side vertex to its new position for `new_radius`.
/// The new vertex sits on the plane at the same angular position
/// around the cylinder axis but at radial distance
/// `cyl.radius + s_radial * new_radius`.
pub(crate) fn move_plane_side_vertex(
    old_pos: &Point3,
    plane: &PlaneGeometry,
    cyl: &CylinderGeometry,
    s_radial: f64,
    new_radius: f64,
) -> Point3 {
    let v = Vec3::new(
        old_pos.x - cyl.origin.x,
        old_pos.y - cyl.origin.y,
        old_pos.z - cyl.origin.z,
    );
    let along = v.dot(&cyl.axis);
    let perp = Vec3::new(
        v.x - along * cyl.axis.x,
        v.y - along * cyl.axis.y,
        v.z - along * cyl.axis.z,
    );
    let perp_len = perp.length();
    let new_radial = cyl.radius + s_radial * new_radius;
    let scale = if perp_len < 1e-12 {
        0.0
    } else {
        new_radial / perp_len
    };
    let lhs_const = Vec3::new(
        (cyl.origin.x - plane.origin.x) + perp.x * scale,
        (cyl.origin.y - plane.origin.y) + perp.y * scale,
        (cyl.origin.z - plane.origin.z) + perp.z * scale,
    );
    let denom = cyl.axis.dot(&plane.normal);
    let new_along = if denom.abs() < 1e-12 {
        along
    } else {
        -lhs_const.dot(&plane.normal) / denom
    };
    Point3::new(
        cyl.origin.x + new_along * cyl.axis.x + perp.x * scale,
        cyl.origin.y + new_along * cyl.axis.y + perp.y * scale,
        cyl.origin.z + new_along * cyl.axis.z + perp.z * scale,
    )
}

/// Move the cylinder-side vertex to its new position. Sits on the
/// cylinder at the same angular position, displaced from the plane
/// along the cylinder axis by `s_plane * new_radius`.
pub(crate) fn move_cyl_side_vertex(
    old_pos: &Point3,
    plane: &PlaneGeometry,
    cyl: &CylinderGeometry,
    s_plane: f64,
    new_radius: f64,
) -> Point3 {
    let v = Vec3::new(
        old_pos.x - cyl.origin.x,
        old_pos.y - cyl.origin.y,
        old_pos.z - cyl.origin.z,
    );
    let along = v.dot(&cyl.axis);
    let perp = Vec3::new(
        v.x - along * cyl.axis.x,
        v.y - along * cyl.axis.y,
        v.z - along * cyl.axis.z,
    );
    let perp_len = perp.length();
    let scale = if perp_len < 1e-12 {
        0.0
    } else {
        cyl.radius / perp_len
    };
    let new_perp = Vec3::new(perp.x * scale, perp.y * scale, perp.z * scale);
    let lhs_const = Vec3::new(
        (cyl.origin.x - plane.origin.x) + new_perp.x,
        (cyl.origin.y - plane.origin.y) + new_perp.y,
        (cyl.origin.z - plane.origin.z) + new_perp.z,
    );
    let denom = cyl.axis.dot(&plane.normal);
    let plane_along = if denom.abs() < 1e-12 {
        along
    } else {
        -lhs_const.dot(&plane.normal) / denom
    };
    let new_along = plane_along + s_plane * new_radius;
    Point3::new(
        cyl.origin.x + new_along * cyl.axis.x + new_perp.x,
        cyl.origin.y + new_along * cyl.axis.y + new_perp.y,
        cyl.origin.z + new_along * cyl.axis.z + new_perp.z,
    )
}

/// Cylinder-axis intersection with the plane — the centre of the
/// plane-side and (when offset by `s_plane * new_radius`) the
/// cyl-side contact circles.
pub(crate) fn project_to_plane_along_axis(plane: &PlaneGeometry, cyl: &CylinderGeometry) -> Point3 {
    let denom = cyl.axis.dot(&plane.normal);
    if denom.abs() < 1e-12 {
        return cyl.origin;
    }
    let diff = Vec3::new(
        plane.origin.x - cyl.origin.x,
        plane.origin.y - cyl.origin.y,
        plane.origin.z - cyl.origin.z,
    );
    let t = diff.dot(&plane.normal) / denom;
    Point3::new(
        cyl.origin.x + t * cyl.axis.x,
        cyl.origin.y + t * cyl.axis.y,
        cyl.origin.z + t * cyl.axis.z,
    )
}

pub(crate) fn project_to_cylinder_axis_at_height(
    cyl: &CylinderGeometry,
    plane: &PlaneGeometry,
    s_plane: f64,
    new_radius: f64,
) -> Point3 {
    let plane_centre = project_to_plane_along_axis(plane, cyl);
    Point3::new(
        plane_centre.x + s_plane * new_radius * cyl.axis.x,
        plane_centre.y + s_plane * new_radius * cyl.axis.y,
        plane_centre.z + s_plane * new_radius * cyl.axis.z,
    )
}

/// Replace `edge`'s curve with a new circular NURBS centred at
/// `centre` with `radius` in the plane perpendicular to `axis`.
pub(crate) fn rebuild_circle_edge(
    brep: &mut BRep,
    edge: EdgeHandle,
    centre: &Point3,
    radius: f64,
    axis: &Vec3,
) {
    let nurbs = build_circle_nurbs(centre, radius, axis);
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs),
        t_start: 0.0,
        t_end: 1.0,
    });
}

/// 9-control-point rational quadratic NURBS for a full circle in
/// the plane perpendicular to `normal`.
pub(crate) fn build_circle_nurbs(centre: &Point3, radius: f64, normal: &Vec3) -> NurbsCurve {
    let n = if normal.length() > 1e-12 {
        let s = 1.0 / normal.length();
        Vec3::new(normal.x * s, normal.y * s, normal.z * s)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let u_seed = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let dot = u_seed.dot(&n);
    let u_raw = Vec3::new(
        u_seed.x - dot * n.x,
        u_seed.y - dot * n.y,
        u_seed.z - dot * n.z,
    );
    let u_len = u_raw.length();
    let u = if u_len > 1e-12 {
        Vec3::new(u_raw.x / u_len, u_raw.y / u_len, u_raw.z / u_len)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let v = n.cross(&u);
    let r = radius;
    let sqrt2 = std::f64::consts::FRAC_1_SQRT_2;
    let pt = |a: f64, b: f64| -> Point3 {
        Point3::new(
            centre.x + a * r * u.x + b * r * v.x,
            centre.y + a * r * u.y + b * r * v.y,
            centre.z + a * r * u.z + b * r * v.z,
        )
    };
    let cps = vec![
        pt(1.0, 0.0),
        pt(1.0, 1.0),
        pt(0.0, 1.0),
        pt(-1.0, 1.0),
        pt(-1.0, 0.0),
        pt(-1.0, -1.0),
        pt(0.0, -1.0),
        pt(1.0, -1.0),
        pt(1.0, 0.0),
    ];
    let weights = vec![1.0, sqrt2, 1.0, sqrt2, 1.0, sqrt2, 1.0, sqrt2, 1.0];
    let knots = vec![
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ];
    NurbsCurve {
        degree: 2,
        knots,
        control_points: cps,
        weights,
    }
}
