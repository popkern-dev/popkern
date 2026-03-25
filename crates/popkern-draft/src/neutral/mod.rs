//! Neutral line computation for draft operations.
//!
//! The neutral line is the curve on a face that remains stationary during
//! the draft rotation. Points above and below the neutral line tilt in
//! opposite directions.

mod helpers;

use crate::error::DraftError;
use helpers::{
    build_line_curve, compute_face_centroid, compute_face_extent, compute_perpendicular,
    solve_point_on_two_planes,
};
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::types::NeutralReference;

/// A computed neutral line on a face.
///
/// The neutral line is the axis of rotation for the draft operation.
/// Points on the face rotate about this line by the draft angle.
#[derive(Debug, Clone)]
pub struct NeutralLine {
    /// 3D curve on the face representing the neutral line.
    pub curve: NurbsCurve,
    /// Unit direction along the neutral line.
    pub direction: Vec3,
    /// A point on the neutral line (rotation center).
    pub origin: Point3,
}

/// Compute the neutral line for a face given the neutral reference.
///
/// # Errors
///
/// Returns [`DraftError::InvalidNeutral`] if the neutral line cannot be
/// determined (e.g., no intersection with the face surface).
/// Returns [`DraftError::InvalidFace`] if the face has no surface geometry.
pub fn compute_neutral_line(
    brep: &BRep,
    face: FaceHandle,
    neutral_ref: &NeutralReference,
    pull_direction: &Vec3,
    tolerance: &Tolerance,
) -> Result<NeutralLine, DraftError> {
    match neutral_ref {
        NeutralReference::Plane { origin, normal } => {
            compute_plane_neutral(brep, face, origin, normal, tolerance)
        }
        NeutralReference::Edge(edge_handle) => {
            compute_edge_neutral(brep, face, *edge_handle, tolerance)
        }
        NeutralReference::Auto => compute_auto_neutral(brep, face, pull_direction, tolerance),
    }
}

/// Compute neutral line from plane intersection with the face surface.
fn compute_plane_neutral(
    brep: &BRep,
    face: FaceHandle,
    neutral_origin: &Point3,
    neutral_normal: &Vec3,
    tolerance: &Tolerance,
) -> Result<NeutralLine, DraftError> {
    let face_data = &brep.faces[face];
    let surface = face_data
        .surface
        .as_ref()
        .ok_or_else(|| DraftError::InvalidFace("face has no surface geometry".into()))?;

    match &surface.geometry {
        SurfaceGeometry::Plane {
            origin: face_origin,
            normal: face_normal,
            ..
        } => compute_planar_neutral(
            brep,
            face,
            face_origin,
            face_normal,
            neutral_origin,
            neutral_normal,
            tolerance,
        ),
        _ => compute_general_neutral(
            brep,
            face,
            surface,
            neutral_origin,
            neutral_normal,
            tolerance,
        ),
    }
}

/// Planar face + plane neutral: analytic intersection.
#[allow(clippy::too_many_arguments)]
fn compute_planar_neutral(
    brep: &BRep,
    face: FaceHandle,
    face_origin: &Point3,
    face_normal: &Vec3,
    neutral_origin: &Point3,
    neutral_normal: &Vec3,
    tolerance: &Tolerance,
) -> Result<NeutralLine, DraftError> {
    let nn = neutral_normal
        .normalized()
        .ok_or_else(|| DraftError::InvalidNeutral("neutral normal is zero".into()))?;
    let fn_norm = face_normal
        .normalized()
        .ok_or_else(|| DraftError::InvalidNeutral("face normal is zero".into()))?;

    let direction = fn_norm.cross(&nn);
    let dir_len = direction.length();
    if dir_len < tolerance.linear {
        return Err(DraftError::InvalidNeutral(
            "neutral plane parallel to face".into(),
        ));
    }
    let direction = direction * (1.0 / dir_len);

    let d1 = fn_norm.dot(&(*face_origin - Point3::new(0.0, 0.0, 0.0)));
    let d2 = nn.dot(&(*neutral_origin - Point3::new(0.0, 0.0, 0.0)));
    let point_on_line = solve_point_on_two_planes(fn_norm, d1, nn, d2, direction, tolerance)?;

    let extent = compute_face_extent(brep, face);
    let p0 = point_on_line - direction * extent;
    let p1 = point_on_line + direction * extent;
    let curve = build_line_curve(p0, p1);

    Ok(NeutralLine {
        curve,
        direction,
        origin: point_on_line,
    })
}

/// General face + plane neutral: surface-surface intersection.
fn compute_general_neutral(
    brep: &BRep,
    face: FaceHandle,
    surface: &popkern_core::geometry::SurfaceOnFace,
    neutral_origin: &Point3,
    neutral_normal: &Vec3,
    tolerance: &Tolerance,
) -> Result<NeutralLine, DraftError> {
    let neutral_plane = SurfaceGeometry::Plane {
        origin: *neutral_origin,
        normal: *neutral_normal,
        u_axis: compute_perpendicular(neutral_normal),
        v_axis: neutral_normal
            .cross(&compute_perpendicular(neutral_normal))
            .normalized()
            .unwrap_or(Vec3::new(0.0, 1.0, 0.0)),
    };

    let ((u1_min, u1_max), (v1_min, v1_max)) = surface.geometry.parameter_range();
    let range = compute_face_extent(brep, face);

    let result = popkern_core::geometry::intersection::intersect_surfaces(
        &surface.geometry,
        (u1_min, u1_max),
        (v1_min, v1_max),
        &neutral_plane,
        (-range, range),
        (-range, range),
        tolerance,
    )?;

    for isect in &result.intersections {
        if let popkern_core::geometry::intersection::SurfaceSurfaceIntersection::Curve(cos) = isect
        {
            let (t0, t1) = cos.curve_3d.parameter_range();
            let mid_t = (t0 + t1) * 0.5;
            let origin = cos.curve_3d.point_at(mid_t);
            let direction = cos
                .curve_3d
                .tangent_at(mid_t)
                .normalized()
                .unwrap_or(Vec3::new(1.0, 0.0, 0.0));

            return Ok(NeutralLine {
                curve: cos.curve_3d.clone(),
                direction,
                origin,
            });
        }
    }

    Err(DraftError::InvalidNeutral(
        "no intersection curve found between face and neutral plane".into(),
    ))
}

/// Compute neutral line from an existing edge.
fn compute_edge_neutral(
    brep: &BRep,
    face: FaceHandle,
    edge_handle: popkern_core::topology::EdgeHandle,
    _tolerance: &Tolerance,
) -> Result<NeutralLine, DraftError> {
    let edge = brep
        .edges
        .get(edge_handle)
        .ok_or_else(|| DraftError::InvalidNeutral("edge handle not found in BRep".into()))?;

    let he_handle = edge.half_edge;
    let he = &brep.half_edges[he_handle];
    let twin = &brep.half_edges[he.twin];

    let he_loop = he.loop_;
    let twin_loop = twin.loop_;

    let face_data = &brep.faces[face];
    let on_face = brep.loops[he_loop].face == face
        || brep.loops[twin_loop].face == face
        || face_data.inner_loops.contains(&he_loop)
        || face_data.inner_loops.contains(&twin_loop);

    if !on_face {
        return Err(DraftError::InvalidNeutral(
            "edge is not on the specified face".into(),
        ));
    }

    let curve_on_edge = edge
        .curve
        .as_ref()
        .ok_or_else(|| DraftError::InvalidNeutral("edge has no curve geometry".into()))?;

    let (t0, t1) = curve_on_edge.geometry.parameter_range();
    let mid_t = (t0 + t1) * 0.5;
    let origin = curve_on_edge.geometry.point_at(mid_t);
    let direction = curve_on_edge
        .geometry
        .tangent_at(mid_t)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));

    let nurbs_curve = match &curve_on_edge.geometry {
        popkern_core::geometry::curve::CurveGeometry::Nurbs(nc) => nc.clone(),
        _ => {
            let start = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
            let end = curve_on_edge.geometry.point_at(curve_on_edge.t_end);
            build_line_curve(start, end)
        }
    };

    Ok(NeutralLine {
        curve: nurbs_curve,
        direction,
        origin,
    })
}

/// Compute neutral line automatically: perpendicular to pull direction at face centroid.
fn compute_auto_neutral(
    brep: &BRep,
    face: FaceHandle,
    pull_direction: &Vec3,
    tolerance: &Tolerance,
) -> Result<NeutralLine, DraftError> {
    let centroid = compute_face_centroid(brep, face);
    let pull_norm = pull_direction
        .normalized()
        .ok_or_else(|| DraftError::InvalidNeutral("pull direction is zero".into()))?;

    let result = compute_plane_neutral(brep, face, &centroid, &pull_norm, tolerance);

    match result {
        Ok(nl) => Ok(nl),
        Err(_) => find_best_edge_neutral(brep, face, &pull_norm, tolerance),
    }
}

/// Fallback: find the face edge whose tangent is most perpendicular to pull.
fn find_best_edge_neutral(
    brep: &BRep,
    face: FaceHandle,
    pull_norm: &Vec3,
    tolerance: &Tolerance,
) -> Result<NeutralLine, DraftError> {
    let face_data = &brep.faces[face];
    let outer_loop = &brep.loops[face_data.outer_loop];

    let mut best_edge = None;
    let mut best_dot = f64::MAX;

    let start_he = outer_loop.half_edge;
    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        let edge_h = he.edge;
        let edge = &brep.edges[edge_h];

        if let Some(ref coe) = edge.curve {
            let (t0, t1) = coe.geometry.parameter_range();
            let mid_t = (t0 + t1) * 0.5;
            let tangent = coe.geometry.tangent_at(mid_t);
            if let Some(t_norm) = tangent.normalized() {
                let dot = t_norm.dot(pull_norm).abs();
                if dot < best_dot {
                    best_dot = dot;
                    best_edge = Some(edge_h);
                }
            }
        }

        current = he.next;
        if current == start_he {
            break;
        }
    }

    match best_edge {
        Some(eh) => compute_edge_neutral(brep, face, eh, tolerance),
        None => Err(DraftError::InvalidNeutral(
            "no suitable edge found for auto neutral".into(),
        )),
    }
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
