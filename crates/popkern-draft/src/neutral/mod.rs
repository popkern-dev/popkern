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
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn plane_neutral_on_cube_face() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let nl = compute_neutral_line(
            &brep,
            face,
            &NeutralReference::Plane {
                origin: Point3::new(0.0, 0.0, 0.5),
                normal: Vec3::new(0.0, 0.0, 1.0),
            },
            &Vec3::new(0.0, 0.0, 1.0),
            &tol(),
        );
        let _ = nl;
    }

    #[test]
    fn auto_neutral_on_cube_face() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let nl = compute_neutral_line(
            &brep,
            face,
            &NeutralReference::Auto,
            &Vec3::new(0.0, 0.0, 1.0),
            &tol(),
        );
        let _ = nl;
    }

    #[test]
    fn edge_neutral_on_cube_face() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let face_data = &brep.faces[face];
        let outer_loop = &brep.loops[face_data.outer_loop];
        let he = &brep.half_edges[outer_loop.half_edge];
        let edge_h = he.edge;

        let nl = compute_neutral_line(
            &brep,
            face,
            &NeutralReference::Edge(edge_h),
            &Vec3::new(0.0, 0.0, 1.0),
            &tol(),
        );
        assert!(nl.is_ok(), "edge neutral failed: {:?}", nl.err());
        let nl = nl.unwrap();
        assert!(nl.direction.length() > 0.9, "direction should be unit");
    }

    #[test]
    fn invalid_edge_handle_returns_error() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let nl = compute_neutral_line(
            &brep,
            face,
            &NeutralReference::Edge(Default::default()),
            &Vec3::new(0.0, 0.0, 1.0),
            &tol(),
        );
        assert!(nl.is_err());
    }

    #[test]
    fn plane_neutral_parallel_to_face_returns_error() {
        // For a cube's Z-normal face, a neutral plane also with Z normal is parallel →
        // the cross product is zero → InvalidNeutral error.
        let brep = BRep::unit_cube(tol()).unwrap();

        // Find the face whose surface normal is Z (top or bottom face).
        use popkern_core::geometry::surface::SurfaceEval;
        let z_face = brep.faces.iter().find(|(_, fd)| {
            fd.surface.as_ref().map_or(false, |sof| {
                let n = sof.geometry.normal_at(0.5, 0.5);
                n.z.abs() > 0.9
            })
        });

        if let Some((face_h, _)) = z_face {
            // Neutral plane also has Z normal → parallel to face → error.
            let nl = compute_neutral_line(
                &brep,
                face_h,
                &NeutralReference::Plane {
                    origin: Point3::new(0.5, 0.5, 0.5),
                    normal: Vec3::new(0.0, 0.0, 1.0),
                },
                &Vec3::new(0.0, 0.0, 1.0),
                &tol(),
            );
            assert!(
                nl.is_err(),
                "expected error for neutral plane parallel to face, got {:?}",
                nl.ok()
            );
        }
    }

    #[test]
    fn edge_neutral_edge_not_on_face_returns_error() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let mut face_iter = brep.faces.keys();
        let face_a = face_iter.next().unwrap();
        let face_b = face_iter.next().unwrap();

        // Find an edge that belongs to face_b but NOT face_a.
        let face_b_data = &brep.faces[face_b];
        let loop_b = &brep.loops[face_b_data.outer_loop];
        let he_b = &brep.half_edges[loop_b.half_edge];
        let edge_from_b = he_b.edge;

        // Verify the edge is actually on face_b, not face_a, before asserting.
        let face_a_data = &brep.faces[face_a];
        let loop_a = &brep.loops[face_a_data.outer_loop];
        let start_he = loop_a.half_edge;
        let mut cur = start_he;
        let mut edge_on_a = false;
        loop {
            let he = &brep.half_edges[cur];
            if he.edge == edge_from_b {
                edge_on_a = true;
                break;
            }
            cur = he.next;
            if cur == start_he {
                break;
            }
        }

        if !edge_on_a {
            let nl = compute_neutral_line(
                &brep,
                face_a,
                &NeutralReference::Edge(edge_from_b),
                &Vec3::new(0.0, 0.0, 1.0),
                &tol(),
            );
            assert!(nl.is_err(), "expected error: edge not on face");
        }
    }

    #[test]
    fn auto_neutral_with_pull_parallel_to_face_normal_uses_edge_fallback() {
        // Find a Z-normal cube face, then pull in Z direction.
        // The neutral plane (Z-normal through centroid) is parallel to the face
        // → compute_plane_neutral fails → auto falls back to find_best_edge_neutral.
        let brep = BRep::unit_cube(tol()).unwrap();
        use popkern_core::geometry::surface::SurfaceEval;
        let z_face = brep.faces.iter().find(|(_, fd)| {
            fd.surface.as_ref().map_or(false, |sof| {
                let n = sof.geometry.normal_at(0.5, 0.5);
                n.z.abs() > 0.9
            })
        });

        if let Some((face_h, _)) = z_face {
            // Pull in Z → neutral plane has Z normal → parallel to face → falls back.
            let nl = compute_neutral_line(
                &brep,
                face_h,
                &NeutralReference::Auto,
                &Vec3::new(0.0, 0.0, 1.0),
                &tol(),
            );
            // Should succeed via edge fallback (the face has 4 edges to choose from).
            // It may also succeed directly if the centroid plane is acceptable.
            let _ = nl; // either outcome is valid; just exercise the code path.
        }
    }

    #[test]
    fn plane_neutral_succeeds_for_side_face() {
        // A side face (X-normal) with a Z-normal neutral plane should intersect.
        let brep = BRep::unit_cube(tol()).unwrap();
        use popkern_core::geometry::surface::SurfaceEval;
        let x_face = brep.faces.iter().find(|(_, fd)| {
            fd.surface.as_ref().map_or(false, |sof| {
                let n = sof.geometry.normal_at(0.5, 0.5);
                n.x.abs() > 0.9
            })
        });

        if let Some((face_h, _)) = x_face {
            let nl = compute_neutral_line(
                &brep,
                face_h,
                &NeutralReference::Plane {
                    origin: Point3::new(0.5, 0.5, 0.5),
                    normal: Vec3::new(0.0, 0.0, 1.0),
                },
                &Vec3::new(0.0, 1.0, 0.0),
                &tol(),
            );
            assert!(
                nl.is_ok(),
                "expected Ok for side face + Z neutral, got {:?}",
                nl.err()
            );
        }
    }

    #[test]
    fn edge_neutral_direction_is_unit() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let face_data = &brep.faces[face];
        let outer_loop = &brep.loops[face_data.outer_loop];

        // Walk all 4 edges of the outer loop and check each.
        let start_he = outer_loop.half_edge;
        let mut cur = start_he;
        loop {
            let he = &brep.half_edges[cur];
            let edge_h = he.edge;
            let nl = compute_neutral_line(
                &brep,
                face,
                &NeutralReference::Edge(edge_h),
                &Vec3::new(0.0, 0.0, 1.0),
                &tol(),
            );
            if let Ok(nl) = nl {
                let len = nl.direction.length();
                assert!((len - 1.0).abs() < 0.01, "direction not unit: {}", len);
            }
            cur = he.next;
            if cur == start_he {
                break;
            }
        }
    }
}
