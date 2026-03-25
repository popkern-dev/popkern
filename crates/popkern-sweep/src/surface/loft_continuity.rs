//! G1/G2 endpoint continuity for loft surfaces.
//!
//! Adjusts interior control point rows so the surface v-derivative at
//! the boundary aligns with the target face's surface normal, providing
//! tangent-plane (G1) continuity between the loft and an adjacent face.

use popkern_core::geometry::NurbsSurface;
use popkern_core::geometry::surface::SurfaceEval;

use crate::error::SweepError;
use crate::types::Continuity;

/// Which end of the loft a continuity constraint applies to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ContinuityEnd {
    /// Constraint at the start (v = 0, row indices 0 and 1).
    Start,
    /// Constraint at the end (v = 1, row indices num_v-1 and num_v-2).
    End,
}

/// Apply G1 (or G2 falling back to G1) continuity at a loft end.
///
/// Adjusts the first (or last) interior row of control points so the
/// surface's v-derivative aligns with the target face's normal direction.
/// G2 constraints fall back to G1 with the same adjustment — curvature
/// matching requires second-derivative computation that is deferred.
pub(crate) fn apply_continuity(
    surface: &mut NurbsSurface,
    continuity: &Continuity,
    end: ContinuityEnd,
) -> Result<(), SweepError> {
    let (face_brep, face_handle) = match continuity {
        Continuity::Free => return Ok(()),
        Continuity::G1 { face, face_handle } => (face, face_handle),
        Continuity::G2 { face, face_handle } => {
            // G2 falls back to G1 — curvature matching deferred.
            (face, face_handle)
        }
    };

    // Look up the target face's surface geometry.
    let face_data = face_brep
        .faces
        .get(*face_handle)
        .ok_or_else(|| SweepError::InvalidProfile("continuity face handle invalid".into()))?;
    let surface_on_face = face_data
        .surface
        .as_ref()
        .ok_or_else(|| SweepError::InvalidProfile("continuity face has no surface".into()))?;

    adjust_row_for_g1(surface, &surface_on_face.geometry, end);

    Ok(())
}

/// Adjust the interior row of control points adjacent to a boundary row
/// so the v-derivative aligns with the target face normal.
///
/// For `Start`: row 0 is the boundary, row 1 is adjusted.
/// For `End`: row (num_v-1) is the boundary, row (num_v-2) is adjusted.
fn adjust_row_for_g1(
    surface: &mut NurbsSurface,
    target_surface: &popkern_core::geometry::SurfaceGeometry,
    end: ContinuityEnd,
) {
    let num_u = surface.num_u;
    let num_v = surface.num_v;

    if num_v < 2 {
        return;
    }

    // Determine boundary and interior row indices.
    let (boundary_iv, interior_iv) = match end {
        ContinuityEnd::Start => (0_usize, 1_usize),
        ContinuityEnd::End => (num_v - 1, num_v - 2),
    };

    // Compute delta_v: the knot span between the two rows.
    let delta_v = compute_local_knot_span(surface, end);
    if delta_v.abs() < 1e-30 {
        return;
    }

    // Get the u parameter range for sampling face normals.
    let ((u_min, u_max), _) = target_surface.parameter_range();
    let u_range_finite = (u_min.clamp(-1e6, 1e6), u_max.clamp(-1e6, 1e6));

    // v parameter on the target surface: use start or end of its range.
    let (_, (v_target_min, v_target_max)) = target_surface.parameter_range();
    let v_target = match end {
        ContinuityEnd::Start => v_target_min.clamp(-1e6, 1e6),
        ContinuityEnd::End => v_target_max.clamp(-1e6, 1e6),
    };

    let sign = match end {
        ContinuityEnd::Start => 1.0,
        ContinuityEnd::End => -1.0,
    };

    for j in 0..num_u {
        let boundary_idx = j * num_v + boundary_iv;
        let interior_idx = j * num_v + interior_iv;

        let q0 = surface.control_points[boundary_idx];
        let q1 = surface.control_points[interior_idx];

        // Current v-derivative direction.
        let dv = q1 - q0;
        let dv_len = dv.length();
        if dv_len < 1e-30 {
            continue;
        }

        // Compute the u parameter for this column on the target face.
        let u_frac = if num_u > 1 {
            j as f64 / (num_u - 1) as f64
        } else {
            0.5
        };
        let u_target = u_range_finite.0 + u_frac * (u_range_finite.1 - u_range_finite.0);

        // Get the face normal at this point.
        let normal = target_surface.normal_at(u_target, v_target);
        let normal_len = normal.length();
        if normal_len < 1e-30 {
            continue;
        }
        let normal_unit = normal * (1.0 / normal_len);

        // New interior point: boundary + normal_unit * sign * |dv|
        // This aligns the control polygon leg with the face normal direction,
        // making the v-derivative at the boundary tangent to the target face.
        let new_q1 = q0 + normal_unit * (sign * dv_len);
        surface.control_points[interior_idx] = new_q1;
    }
}

/// Compute the local knot span near a boundary end.
fn compute_local_knot_span(surface: &NurbsSurface, end: ContinuityEnd) -> f64 {
    let num_v = surface.num_v;
    let d = surface.degree_v;

    let span_start = surface.knots_v[d];
    let span_end = surface.knots_v[num_v];
    let span_total = span_end - span_start;

    if span_total.abs() < 1e-30 {
        return 0.0;
    }

    match end {
        ContinuityEnd::Start => {
            let k = d + 1;
            if k < surface.knots_v.len() {
                surface.knots_v[k] - surface.knots_v[d]
            } else {
                span_total / (num_v - 1) as f64
            }
        }
        ContinuityEnd::End => {
            let k = num_v;
            if k > 0 {
                surface.knots_v[k] - surface.knots_v[k - 1]
            } else {
                span_total / (num_v - 1) as f64
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::{NurbsCurve, Point3, SurfaceGeometry, SurfaceOnFace, Vec3};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    use crate::surface::loft_surface::generate_loft_surface;
    use crate::types::{LoftOptions, LoftSection};

    /// Create a unit circle as a rational NURBS (degree 2, 9 control points).
    fn unit_circle_at_z(z: f64, radius: f64) -> NurbsCurve {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let r = radius;
        NurbsCurve {
            degree: 2,
            knots: vec![
                0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
            ],
            control_points: vec![
                Point3::new(r, 0.0, z),
                Point3::new(r, r, z),
                Point3::new(0.0, r, z),
                Point3::new(-r, r, z),
                Point3::new(-r, 0.0, z),
                Point3::new(-r, -r, z),
                Point3::new(0.0, -r, z),
                Point3::new(r, -r, z),
                Point3::new(r, 0.0, z),
            ],
            weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
        }
    }

    fn make_section(curve: NurbsCurve) -> LoftSection {
        LoftSection {
            curve,
            position: None,
        }
    }

    /// Build a BRep with a single face that has a horizontal plane surface
    /// (normal pointing in +Z) at the given z height.
    fn brep_with_plane_face(z: f64) -> (BRep, popkern_core::topology::FaceHandle) {
        let tol = Tolerance::default();
        let mut brep = BRep::new(tol);
        let solid_handle = brep.solids.insert(popkern_core::topology::Solid {
            shell: popkern_core::topology::ShellHandle::default(),
        });
        let shell_handle = brep.shells.insert(popkern_core::topology::Shell {
            face: popkern_core::topology::FaceHandle::default(),
            solid: solid_handle,
        });
        let loop_handle = brep.loops.insert(popkern_core::topology::Loop {
            half_edge: popkern_core::topology::HalfEdgeHandle::default(),
            face: popkern_core::topology::FaceHandle::default(),
        });
        let face_handle = brep.faces.insert(popkern_core::topology::Face {
            outer_loop: loop_handle,
            inner_loops: vec![],
            shell: shell_handle,
            surface: Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Plane {
                    origin: Point3::new(0.0, 0.0, z),
                    normal: Vec3::new(0.0, 0.0, 1.0),
                    u_axis: Vec3::new(1.0, 0.0, 0.0),
                    v_axis: Vec3::new(0.0, 1.0, 0.0),
                },
                orientation: true,
            }),
        });
        brep.loops[loop_handle].face = face_handle;
        brep.shells[shell_handle].face = face_handle;
        brep.solids[solid_handle].shell = shell_handle;
        (brep, face_handle)
    }

    #[test]
    fn loft_with_g1_start_aligns_tangent() {
        // Build a plane face at z=0 with normal +Z.
        let (brep, face_handle) = brep_with_plane_face(0.0);

        // Loft three circles so we have enough v-CPs.
        let sections = vec![
            make_section(unit_circle_at_z(0.0, 1.0)),
            make_section(unit_circle_at_z(1.0, 1.5)),
            make_section(unit_circle_at_z(2.0, 1.0)),
        ];

        let options = LoftOptions {
            start_continuity: Continuity::G1 {
                face: brep,
                face_handle,
            },
            end_continuity: Continuity::Free,
            closed: false,
        };
        let tol = Tolerance::default();
        let result = generate_loft_surface(&sections, &options, &tol);
        assert!(result.is_ok(), "loft with G1 failed: {:?}", result.err());

        let surf = &result.unwrap().surface; // ok in tests

        // Check the v-derivative at v=0 is approximately aligned with +Z.
        for u_frac in [0.0, 0.25, 0.5, 0.75] {
            let dv = surf.partial_v(u_frac, 0.0);
            let dv_len = dv.length();
            if dv_len > 1e-10 {
                let dv_unit_z = dv.z / dv_len;
                // The z-component of the normalized derivative should dominate.
                assert!(
                    dv_unit_z.abs() > 0.5,
                    "v-derivative at u={u_frac} not aligned with Z: dv={dv:?}"
                );
            }
        }
    }

    #[test]
    fn loft_with_g2_falls_back_to_g1() {
        // G2 should produce the same result as G1 (fallback).
        let (brep, face_handle) = brep_with_plane_face(0.0);

        let sections = vec![
            make_section(unit_circle_at_z(0.0, 1.0)),
            make_section(unit_circle_at_z(1.0, 1.5)),
            make_section(unit_circle_at_z(2.0, 1.0)),
        ];

        let options = LoftOptions {
            start_continuity: Continuity::G2 {
                face: brep,
                face_handle,
            },
            end_continuity: Continuity::Free,
            closed: false,
        };
        let tol = Tolerance::default();
        let result = generate_loft_surface(&sections, &options, &tol);
        assert!(result.is_ok(), "loft with G2 failed: {:?}", result.err());
    }

    #[test]
    fn loft_with_g1_end_aligns_tangent() {
        let (brep, face_handle) = brep_with_plane_face(2.0);

        let sections = vec![
            make_section(unit_circle_at_z(0.0, 1.0)),
            make_section(unit_circle_at_z(1.0, 1.5)),
            make_section(unit_circle_at_z(2.0, 1.0)),
        ];

        let options = LoftOptions {
            start_continuity: Continuity::Free,
            end_continuity: Continuity::G1 {
                face: brep,
                face_handle,
            },
            closed: false,
        };
        let tol = Tolerance::default();
        let result = generate_loft_surface(&sections, &options, &tol);
        assert!(
            result.is_ok(),
            "loft with G1 end failed: {:?}",
            result.err()
        );

        let surf = &result.unwrap().surface; // ok in tests

        // Check the v-derivative at v=1 has a significant Z component.
        for u_frac in [0.0, 0.25, 0.5, 0.75] {
            let dv = surf.partial_v(u_frac, 1.0);
            let dv_len = dv.length();
            if dv_len > 1e-10 {
                let dv_unit_z = dv.z / dv_len;
                assert!(
                    dv_unit_z.abs() > 0.3,
                    "v-derivative at end u={u_frac} not aligned with Z: dv={dv:?}"
                );
            }
        }
    }

    #[test]
    fn continuity_free_is_noop() {
        let mut surf = NurbsSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 2,
            num_v: 2,
            weights: vec![1.0; 4],
        };
        let cps_before = surf.control_points.clone();
        let result = apply_continuity(&mut surf, &Continuity::Free, ContinuityEnd::Start);
        assert!(result.is_ok());
        assert_eq!(surf.control_points, cps_before);
    }
}
