//! Internal geometry computation helpers for fillet and chamfer operations.
//!
//! Extracted from `lib.rs` to keep source files within the 500-line limit.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
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
pub(crate) type OverlapTrimEntry =
    (usize, popkern_core::geometry::NurbsCurve, popkern_core::geometry::NurbsCurve);

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

    let (u1_range, v1_range) = surface_param_range(&analysis.surface_1);
    let (u2_range, v2_range) = surface_param_range(&analysis.surface_2);

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

    Ok(EdgeGeometry {
        analysis: rebuild_analysis(analysis, brep)?,
        surface: chamfer_result.surface,
        trim_3d_1: chamfer_result.trim_curve_1,
        trim_3d_2: chamfer_result.trim_curve_2,
        pcurve_1: chamfer_result.pcurve_1,
        pcurve_2: chamfer_result.pcurve_2,
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

/// Get finite parameter ranges for a surface.
pub(crate) fn surface_param_range(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
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
pub(crate) fn face_handle_to_index(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
) -> usize {
    brep.faces
        .keys()
        .enumerate()
        .find(|(_, k)| *k == face)
        .map(|(i, _)| i)
        .unwrap_or(0)
}
