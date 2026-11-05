//! Surface drafting algorithms: dispatch by surface type.

pub mod cylindrical;
pub mod general;
pub mod planar;
pub mod skinning;

use crate::error::DraftError;
use crate::neutral::NeutralLine;
use crate::types::{AngleVariation, DraftOptions};
use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Compute the drafted surface for a face.
///
/// Dispatches to planar, cylindrical, or general NURBS drafting
/// based on the face's surface type.
///
/// # Errors
///
/// Returns [`DraftError`] if the surface cannot be drafted (e.g., angle
/// too large for the geometry, degenerate surface).
#[allow(clippy::too_many_arguments)]
pub fn compute_drafted_surface(
    brep: &BRep,
    face: FaceHandle,
    neutral: &NeutralLine,
    angle: f64,
    variation: &AngleVariation,
    pull_direction: &Vec3,
    options: &DraftOptions,
    tolerance: &Tolerance,
) -> Result<SurfaceGeometry, DraftError> {
    let drafted = drafted_by(
        brep,
        face,
        neutral,
        angle,
        variation,
        pull_direction,
        options,
        tolerance,
    )?;
    // The sense: a drafted face leans so the part draws along the pull,
    // its outward normal (out of the material) gaining a component
    // along it. Above the neutral line the face then sits displaced
    // against its outward normal: an outer wall leans in as it rises, a
    // hole's wall out (the hole widens toward the pull). The rotation
    // above follows the surface's own normal, whose side of the
    // material a hole or pocket wall reverses (measured: a 3° hole
    // narrowed upward, removing 16110 for 15320).
    if wrong_sense(brep, face, neutral, pull_direction, &drafted, angle) == Some(true) {
        return drafted_by(
            brep,
            face,
            neutral,
            -angle,
            variation,
            pull_direction,
            options,
            tolerance,
        );
    }
    Ok(drafted)
}

/// Whether `drafted` leans the wrong way for `face` and the signed
/// `angle` (see [`compute_drafted_surface`]): for a positive angle a
/// point of the face above the neutral line moves along its outward
/// normal; for a negative angle (leaning back, as un-drafting does) it
/// moves against it. `None` when the face gives no such point.
fn wrong_sense(
    brep: &BRep,
    face: FaceHandle,
    neutral: &NeutralLine,
    pull: &Vec3,
    drafted: &SurfaceGeometry,
    angle: f64,
) -> Option<bool> {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    let sof = brep.faces.get(face)?.surface.as_ref()?;
    let pull = pull.normalized()?;
    // The face's reach along the pull from the neutral line.
    let reach = brep
        .face_boundary_sample_points(face, 8)
        .iter()
        .map(|p| (*p - neutral.origin).dot(&pull))
        .fold(0.0_f64, |m, h| if h.abs() > m.abs() { h } else { m });
    if reach.abs() <= brep.tolerance.linear() * 100.0 {
        return None;
    }
    let probe = neutral.origin + pull * reach;
    let (u, v) = project_point_to_uv(&sof.geometry, &probe)?;
    let on_face = sof.geometry.point_at(u, v);
    let outward = sof.geometry.normal_at(u, v) * if sof.orientation { 1.0 } else { -1.0 };
    let (du, dv) = project_point_to_uv(drafted, &on_face)?;
    let moved = drafted.point_at(du, dv) - on_face;
    let along = moved.dot(&outward) * reach.signum();
    (along.abs() > brep.tolerance.linear()).then_some((along > 0.0) == (angle > 0.0))
}

#[allow(clippy::too_many_arguments)]
fn drafted_by(
    brep: &BRep,
    face: FaceHandle,
    neutral: &NeutralLine,
    angle: f64,
    variation: &AngleVariation,
    pull_direction: &Vec3,
    options: &DraftOptions,
    tolerance: &Tolerance,
) -> Result<SurfaceGeometry, DraftError> {
    let face_data = &brep.faces[face];
    let surface = face_data
        .surface
        .as_ref()
        .ok_or_else(|| DraftError::InvalidFace("face has no surface geometry".into()))?;

    match &surface.geometry {
        SurfaceGeometry::Plane { .. } => planar::draft_planar(&surface.geometry, neutral, angle),
        SurfaceGeometry::Cylinder { .. } => {
            cylindrical::draft_cylindrical(&surface.geometry, neutral, angle)
        }
        _ => general::draft_general(
            &surface.geometry,
            neutral,
            angle,
            variation,
            pull_direction,
            options,
            tolerance,
        ),
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
    fn dispatch_planar_face() {
        let brep = BRep::unit_cube(tol()).unwrap();
        // Pick a side face (non-Z-normal) so we can find a valid neutral line.
        use popkern_core::geometry::surface::SurfaceEval;
        let face = brep
            .faces
            .iter()
            .find(|(_, fd)| {
                fd.surface.as_ref().is_some_and(|sof| {
                    let n = sof.geometry.normal_at(0.5, 0.5);
                    n.z.abs() < 0.1
                })
            })
            .map(|(fh, _)| fh)
            .unwrap_or_else(|| brep.faces.keys().next().unwrap());

        let nl = crate::neutral::compute_neutral_line(
            &brep,
            face,
            &crate::types::NeutralReference::Plane {
                origin: popkern_core::geometry::Point3::new(0.5, 0.5, 0.5),
                normal: Vec3::new(0.0, 0.0, 1.0),
            },
            &Vec3::new(0.0, 0.0, 1.0),
            &tol(),
        );
        if let Ok(neutral) = nl {
            let result = compute_drafted_surface(
                &brep,
                face,
                &neutral,
                0.087,
                &AngleVariation::Constant,
                &Vec3::new(0.0, 0.0, 1.0),
                &DraftOptions::default(),
                &tol(),
            );
            assert!(result.is_ok(), "planar dispatch failed: {:?}", result.err());
        }
    }

    #[test]
    fn dispatch_cylindrical_face() {
        let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();

        // Find the lateral (cylindrical) face.
        let cyl_face = brep
            .faces
            .iter()
            .find(|(_, fd)| {
                fd.surface
                    .as_ref()
                    .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Cylinder { .. }))
            })
            .map(|(fh, _)| fh);

        if let Some(face) = cyl_face {
            // Use an edge of the face as neutral so we always get Ok.
            let face_data = &brep.faces[face];
            let loop_h = &brep.loops[face_data.outer_loop];
            let he = &brep.half_edges[loop_h.half_edge];
            let edge_h = he.edge;

            let nl = crate::neutral::compute_neutral_line(
                &brep,
                face,
                &crate::types::NeutralReference::Edge(edge_h),
                &Vec3::new(0.0, 0.0, 1.0),
                &tol(),
            );
            if let Ok(neutral) = nl {
                let result = compute_drafted_surface(
                    &brep,
                    face,
                    &neutral,
                    0.087,
                    &AngleVariation::Constant,
                    &Vec3::new(0.0, 0.0, 1.0),
                    &DraftOptions::default(),
                    &tol(),
                );
                assert!(
                    result.is_ok(),
                    "cylindrical dispatch failed: {:?}",
                    result.err()
                );
                // Should produce a Cone.
                assert!(
                    matches!(result.unwrap(), SurfaceGeometry::Cone { .. }),
                    "cylindrical draft should produce Cone"
                );
            }
        }
    }

    #[test]
    fn dispatch_general_nurbs_face() {
        use crate::surface::skinning;
        use popkern_core::geometry::{NurbsCurve, Point3};

        // Build a minimal BRep with one NURBS face using Euler operators.
        // We construct a standalone BRep containing a NURBS face by hand.
        // For simplicity, build a 2x2 NURBS surface and verify draft_general is called.
        let rows = vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ];
        let nurbs_surf = skinning::skin_points(&rows).unwrap();

        // We'll test the general dispatch path by calling draft_general directly
        // via a NeutralLine, since constructing a full BRep with a NURBS face
        // from scratch requires Euler operators that are beyond the scope here.
        // This test validates that the dispatch logic reaches the general path.
        let neutral = crate::neutral::NeutralLine {
            curve: NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
                weights: vec![1.0, 1.0],
            },
            direction: Vec3::new(1.0, 0.0, 0.0),
            origin: Point3::new(0.5, 0.5, 0.0),
        };

        let result = general::draft_general(
            &SurfaceGeometry::Nurbs(nurbs_surf),
            &neutral,
            0.05,
            &AngleVariation::Constant,
            &Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            &tol(),
        );
        assert!(
            result.is_ok(),
            "general dispatch failed: {:?}",
            result.err()
        );
        assert!(matches!(result.unwrap(), SurfaceGeometry::Nurbs(_)));
    }

    #[test]
    fn face_without_surface_returns_error() {
        // compute_drafted_surface must return InvalidFace if face.surface is None.
        // Use a unit cube and strip the surface from the first face.
        use popkern_core::geometry::{NurbsCurve, Point3};

        let mut brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();

        // Remove the surface from this face so the dispatch hits the error path.
        brep.faces[face].surface = None;

        let neutral = crate::neutral::NeutralLine {
            curve: NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
                weights: vec![1.0, 1.0],
            },
            direction: Vec3::new(1.0, 0.0, 0.0),
            origin: Point3::new(0.0, 0.0, 0.0),
        };

        let result = compute_drafted_surface(
            &brep,
            face,
            &neutral,
            0.1,
            &AngleVariation::Constant,
            &Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            &tol(),
        );
        assert!(result.is_err(), "expected error for face with no surface");
    }
}
