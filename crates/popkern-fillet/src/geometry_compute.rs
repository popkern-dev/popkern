//! Internal geometry computation helpers for fillet and chamfer operations.
//!
//! Extracted from `lib.rs` to keep source files within the 500-line limit.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::analysis;
use crate::center_curve;
use crate::chamfer_curve;
use crate::construct::EdgeGeometry;
use crate::error::BlendFailure;
use crate::fillet_surface;
use crate::trim;
use crate::types::{ChamferKind, FilletRadius};

/// Trim curve triplet used by the overlap detector: (face_index, trim_on_face, other_trim).
pub(crate) type OverlapTrimEntry = (
    usize,
    popkern_core::geometry::NurbsCurve,
    popkern_core::geometry::NurbsCurve,
);

/// Per-face collection of trim curves for face-disappearance detection.
pub(crate) type FaceTrimMap =
    std::collections::HashMap<usize, Vec<popkern_core::geometry::NurbsCurve>>;

/// Compute fillet geometry (center curve, fillet surface, trim curves) for one edge.
pub(crate) fn compute_fillet_geometry(
    analysis: &analysis::EdgeAnalysis,
    radius: &FilletRadius,
    brep: &BRep,
    tolerance: &Tolerance,
) -> Result<EdgeGeometry, BlendFailure> {
    let edge_midpoint = get_edge_midpoint(brep, analysis)?;

    let (u1_range, v1_range) = face_param_range(brep, analysis.face_1, &analysis.surface_1);
    let (u2_range, v2_range) = face_param_range(brep, analysis.face_2, &analysis.surface_2);

    let center_result = match radius {
        FilletRadius::Constant(r) => center_curve::constant_radius_center_curve(
            &analysis.surface_1,
            u1_range,
            v1_range,
            &analysis.surface_2,
            u2_range,
            v2_range,
            *r,
            &edge_midpoint,
            tolerance,
            analysis.orientation_1,
            analysis.orientation_2,
        )
        .map_err(|e| BlendFailure::IntersectionFailed { reason: e })?,
        FilletRadius::Variable { start, end } => center_curve::variable_radius_center_curve(
            &analysis.surface_1,
            u1_range,
            v1_range,
            &analysis.surface_2,
            u2_range,
            v2_range,
            *start,
            *end,
            &edge_midpoint,
            tolerance,
            analysis.orientation_1,
            analysis.orientation_2,
        )
        .map_err(|e| BlendFailure::IntersectionFailed { reason: e })?,
        FilletRadius::ByPoints(pts) => {
            let (r_start, r_end) = if pts.len() >= 2 {
                (pts[0].1, pts[pts.len() - 1].1)
            } else if pts.len() == 1 {
                (pts[0].1, pts[0].1)
            } else {
                return Err(BlendFailure::RadiusTooLarge { max_radius: 0.0 });
            };
            center_curve::variable_radius_center_curve(
                &analysis.surface_1,
                u1_range,
                v1_range,
                &analysis.surface_2,
                u2_range,
                v2_range,
                r_start,
                r_end,
                &edge_midpoint,
                tolerance,
                analysis.orientation_1,
                analysis.orientation_2,
            )
            .map_err(|e| BlendFailure::IntersectionFailed { reason: e })?
        }
    };

    let effective_radius = match radius {
        FilletRadius::Constant(r) => *r,
        FilletRadius::Variable { start, end } => (*start + *end) * 0.5,
        FilletRadius::ByPoints(pts) => {
            if pts.is_empty() {
                return Err(BlendFailure::RadiusTooLarge { max_radius: 0.0 });
            }
            pts.iter().map(|(_, r)| *r).sum::<f64>() / pts.len() as f64
        }
    };

    let fillet_result = fillet_surface::build_fillet_surface(
        &analysis.surface_1,
        &analysis.surface_2,
        &center_result,
        effective_radius,
        8, // arc samples
    )
    .map_err(|e| BlendFailure::IntersectionFailed { reason: e })?;

    let trim_curves = trim::compute_trim_curves(
        &analysis.surface_1,
        &center_result.footprint_1,
        &analysis.surface_2,
        &center_result.footprint_2,
    )
    .map_err(|e| BlendFailure::IntersectionFailed { reason: e })?;

    Ok(EdgeGeometry {
        analysis: rebuild_analysis(analysis, brep)?,
        surface: fillet_result.surface,
        trim_3d_1: trim_curves.curve_3d_1,
        trim_3d_2: trim_curves.curve_3d_2,
        pcurve_1: trim_curves.pcurve_1,
        pcurve_2: trim_curves.pcurve_2,
        effective_radius,
    })
}

/// Compute chamfer geometry for one edge.
pub(crate) fn compute_chamfer_geometry(
    analysis: &analysis::EdgeAnalysis,
    kind: &ChamferKind,
    brep: &BRep,
    tolerance: &Tolerance,
) -> Result<EdgeGeometry, BlendFailure> {
    let chamfer_result = chamfer_curve::compute_chamfer(analysis, brep, kind, tolerance)
        .map_err(|e| BlendFailure::IntersectionFailed { reason: e })?;

    let effective_radius = match kind {
        ChamferKind::Equal(d) => *d,
        ChamferKind::TwoDistance { dist_1, dist_2, .. } => (*dist_1).max(*dist_2),
        ChamferKind::DistanceAngle { distance, .. } => *distance,
    };

    Ok(EdgeGeometry {
        analysis: rebuild_analysis(analysis, brep)?,
        surface: chamfer_result.surface,
        trim_3d_1: chamfer_result.trim_curve_1,
        trim_3d_2: chamfer_result.trim_curve_2,
        pcurve_1: chamfer_result.pcurve_1,
        pcurve_2: chamfer_result.pcurve_2,
        effective_radius,
    })
}

/// Get the 3D midpoint of an edge for center-curve selection.
pub(crate) fn get_edge_midpoint(
    brep: &BRep,
    analysis: &analysis::EdgeAnalysis,
) -> Result<popkern_core::geometry::Point3, BlendFailure> {
    let edge = brep
        .edges
        .get(analysis.edge)
        .ok_or_else(|| BlendFailure::InvalidTopology {
            reason: format!("edge {:?} not found", analysis.edge),
        })?;

    if let Some(curve) = &edge.curve {
        let t_mid = (curve.t_start + curve.t_end) * 0.5;
        Ok(curve.geometry.point_at(t_mid))
    } else {
        let v1 = brep.vertices[analysis.vertex_start].position;
        let v2 = brep.vertices[analysis.vertex_end].position;
        Ok(popkern_core::geometry::Point3::new(
            (v1.x + v2.x) * 0.5,
            (v1.y + v2.y) * 0.5,
            (v1.z + v2.z) * 0.5,
        ))
    }
}

/// Get finite parameter ranges for a surface, bounded by face vertex extents.
///
/// Uses the face's outer-loop vertices projected into UV space (via the core
/// `bounded_parameter_range` helper) instead of the old ±10 fallback.  The core
/// function applies 50 % padding (min 1.0) which gives intersection solvers
/// enough headroom while keeping the domain proportional to the actual face.
pub(crate) fn face_param_range(
    brep: &BRep,
    face_handle: popkern_core::topology::FaceHandle,
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    let positions = brep.face_outer_loop_positions(face_handle);
    popkern_core::geometry::uv_projection::bounded_parameter_range(&positions, surface)
}

/// Get finite parameter ranges for a surface (standalone, without face context).
///
/// Falls back to ±10 for infinite dimensions. Prefer [`face_param_range`] when
/// a face handle is available.
#[cfg(test)]
pub(crate) fn surface_param_range(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    use popkern_core::geometry::surface::SurfaceEval;
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let u_lo = if u_min.is_finite() { u_min } else { -10.0 };
    let u_hi = if u_max.is_finite() { u_max } else { 10.0 };
    let v_lo = if v_min.is_finite() { v_min } else { -10.0 };
    let v_hi = if v_max.is_finite() { v_max } else { 10.0 };
    ((u_lo, u_hi), (v_lo, v_hi))
}

/// Rebuild an `EdgeAnalysis` from a reference (since `EdgeAnalysis` is not `Clone`).
pub(crate) fn rebuild_analysis(
    src: &analysis::EdgeAnalysis,
    _brep: &BRep,
) -> Result<analysis::EdgeAnalysis, BlendFailure> {
    Ok(analysis::EdgeAnalysis {
        edge: src.edge,
        face_1: src.face_1,
        face_2: src.face_2,
        surface_1: src.surface_1.clone(),
        surface_2: src.surface_2.clone(),
        orientation_1: src.orientation_1,
        orientation_2: src.orientation_2,
        convexity: src.convexity,
        vertex_start: src.vertex_start,
        vertex_end: src.vertex_end,
    })
}

/// Build the inputs needed for [`crate::overlap::detect_overlaps`] and
/// [`crate::overlap::detect_face_disappearance`] from computed edge geometries.
///
/// Returns overlap trim tuples and a map of face-index to trim curves.
pub(crate) fn build_overlap_inputs(
    edge_geometries: &[EdgeGeometry],
    brep: &BRep,
) -> (Vec<OverlapTrimEntry>, FaceTrimMap) {
    let mut overlap_trims = Vec::new();
    let mut face_trim_map: std::collections::HashMap<
        usize,
        Vec<popkern_core::geometry::NurbsCurve>,
    > = std::collections::HashMap::new();

    for eg in edge_geometries {
        let f1_idx = face_handle_to_index(brep, eg.analysis.face_1);
        let f2_idx = face_handle_to_index(brep, eg.analysis.face_2);

        overlap_trims.push((f1_idx, eg.trim_3d_1.clone(), eg.trim_3d_2.clone()));
        overlap_trims.push((f2_idx, eg.trim_3d_2.clone(), eg.trim_3d_1.clone()));

        face_trim_map
            .entry(f1_idx)
            .or_default()
            .push(eg.trim_3d_1.clone());
        face_trim_map
            .entry(f2_idx)
            .or_default()
            .push(eg.trim_3d_2.clone());
    }

    (overlap_trims, face_trim_map)
}

/// Map a [`popkern_core::topology::FaceHandle`] to a numeric index for the overlap API.
pub(crate) fn face_handle_to_index(brep: &BRep, face: popkern_core::topology::FaceHandle) -> usize {
    brep.faces
        .keys()
        .enumerate()
        .find(|(_, k)| *k == face)
        .map(|(i, _)| i)
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::geometry::{Point3, Vec3};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    fn make_cube_and_analysis() -> (BRep, analysis::EdgeAnalysis) {
        let brep = BRep::unit_cube(tol()).unwrap();
        let edge = brep.edges.keys().next().unwrap();
        let mut analyses = analysis::analyze_edges(&brep, &[edge], false).unwrap();
        (brep, analyses.remove(0))
    }

    #[test]
    fn surface_param_range_plane_is_infinite_clamped() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let ((u_lo, u_hi), (v_lo, v_hi)) = surface_param_range(&plane);
        // Plane has infinite parameter range, so we clamp to ±10.
        assert_eq!(u_lo, -10.0);
        assert_eq!(u_hi, 10.0);
        assert_eq!(v_lo, -10.0);
        assert_eq!(v_hi, 10.0);
    }

    #[test]
    fn surface_param_range_cylinder_is_bounded() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let ((u_lo, u_hi), (v_lo, v_hi)) = surface_param_range(&cyl);
        // Cylinder: u (angle) in [0, 2π], v is infinite → clamped.
        assert!(u_lo.is_finite());
        assert!(u_hi.is_finite());
        // v axis is infinite for a cylinder, so we clamp to ±10.
        assert_eq!(v_lo, -10.0);
        assert_eq!(v_hi, 10.0);
    }

    #[test]
    fn rebuild_analysis_preserves_handles() {
        let (brep, analysis) = make_cube_and_analysis();
        let rebuilt = rebuild_analysis(&analysis, &brep).unwrap();
        assert_eq!(rebuilt.edge, analysis.edge);
        assert_eq!(rebuilt.face_1, analysis.face_1);
        assert_eq!(rebuilt.face_2, analysis.face_2);
        assert_eq!(rebuilt.convexity, analysis.convexity);
        assert_eq!(rebuilt.vertex_start, analysis.vertex_start);
        assert_eq!(rebuilt.vertex_end, analysis.vertex_end);
    }

    #[test]
    fn face_handle_to_index_returns_valid_index() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let faces: Vec<_> = brep.faces.keys().collect();
        for (expected_i, &face) in faces.iter().enumerate() {
            let idx = face_handle_to_index(&brep, face);
            // The index should be within the range of faces.
            assert!(idx < brep.faces.len(), "index {idx} out of range");
            // The face at that index should match.
            assert_eq!(
                brep.faces.keys().nth(idx).unwrap(),
                face,
                "index mismatch for face {expected_i}"
            );
        }
    }

    #[test]
    fn face_handle_to_index_unknown_face_returns_zero() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let bogus = popkern_core::topology::FaceHandle::default();
        let idx = face_handle_to_index(&brep, bogus);
        assert_eq!(idx, 0, "unknown face should map to index 0");
    }

    #[test]
    fn get_edge_midpoint_with_curve() {
        let (brep, analysis) = make_cube_and_analysis();
        let midpoint = get_edge_midpoint(&brep, &analysis).unwrap();
        // The midpoint should be within the bounding box of the cube [0,1]^3.
        assert!(midpoint.x >= -0.01 && midpoint.x <= 1.01);
        assert!(midpoint.y >= -0.01 && midpoint.y <= 1.01);
        assert!(midpoint.z >= -0.01 && midpoint.z <= 1.01);
    }

    #[test]
    fn build_overlap_inputs_empty_geometries() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let (overlap_trims, face_trim_map) = build_overlap_inputs(&[], &brep);
        assert!(overlap_trims.is_empty());
        assert!(face_trim_map.is_empty());
    }

    #[test]
    fn compute_fillet_geometry_variable_radius() {
        let (brep, analysis) = make_cube_and_analysis();
        let radius = crate::types::FilletRadius::Variable {
            start: 0.08,
            end: 0.12,
        };
        let result = compute_fillet_geometry(&analysis, &radius, &brep, &tol());
        assert!(
            result.is_ok(),
            "variable radius fillet geometry should compute: {:?}",
            result.err()
        );
    }

    #[test]
    fn compute_fillet_geometry_by_points_two_points_exercises_path() {
        // Two-point ByPoints exercises the ByPoints branch in compute_fillet_geometry.
        // The result may fail (degenerate curve), but the branch is reached.
        let (brep, analysis) = make_cube_and_analysis();
        let radius = crate::types::FilletRadius::ByPoints(vec![(0.0, 0.1), (1.0, 0.1)]);
        // Either Ok or Err is acceptable — we just confirm no panic and that the
        // ByPoints code branch fires (particularly r_start/r_end extraction).
        let _ = compute_fillet_geometry(&analysis, &radius, &brep, &tol());
    }

    #[test]
    fn compute_fillet_geometry_by_points_single_point_exercises_path() {
        let (brep, analysis) = make_cube_and_analysis();
        let radius = crate::types::FilletRadius::ByPoints(vec![(0.5, 0.1)]);
        // Single-point ByPoints uses pts[0].1 for both r_start and r_end.
        let _ = compute_fillet_geometry(&analysis, &radius, &brep, &tol());
    }

    #[test]
    fn compute_fillet_geometry_by_points_empty_fails() {
        let (brep, analysis) = make_cube_and_analysis();
        let radius = crate::types::FilletRadius::ByPoints(vec![]);
        let result = compute_fillet_geometry(&analysis, &radius, &brep, &tol());
        assert!(result.is_err(), "empty ByPoints should fail");
    }

    #[test]
    fn compute_chamfer_geometry_equal() {
        let (brep, analysis) = make_cube_and_analysis();
        let kind = crate::types::ChamferKind::Equal(0.1);
        let result = compute_chamfer_geometry(&analysis, &kind, &brep, &tol());
        assert!(
            result.is_ok(),
            "equal chamfer geometry should compute: {:?}",
            result.err()
        );
    }

    #[test]
    fn compute_chamfer_geometry_two_distance() {
        let (brep, analysis) = make_cube_and_analysis();
        let kind = crate::types::ChamferKind::TwoDistance {
            face_1: analysis.face_1,
            dist_1: 0.1,
            dist_2: 0.2,
        };
        let result = compute_chamfer_geometry(&analysis, &kind, &brep, &tol());
        assert!(
            result.is_ok(),
            "two-distance chamfer geometry should compute: {:?}",
            result.err()
        );
    }

    #[test]
    fn compute_chamfer_geometry_distance_angle() {
        let (brep, analysis) = make_cube_and_analysis();
        let kind = crate::types::ChamferKind::DistanceAngle {
            face: analysis.face_1,
            distance: 0.1,
            angle: std::f64::consts::FRAC_PI_4,
        };
        let result = compute_chamfer_geometry(&analysis, &kind, &brep, &tol());
        assert!(
            result.is_ok(),
            "distance-angle chamfer geometry should compute: {:?}",
            result.err()
        );
    }
}
