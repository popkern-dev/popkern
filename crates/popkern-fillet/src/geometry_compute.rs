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
    g2_cross_section: bool,
) -> Result<EdgeGeometry, BlendFailure> {
    let edge_midpoint = get_edge_midpoint(brep, analysis)?;

    let (u1_range, v1_range) = face_param_range(brep, analysis.face_1, &analysis.surface_1);
    let (u2_range, v2_range) = face_param_range(brep, analysis.face_2, &analysis.surface_2);

    // Pre-flight: the rolling-ball radius must not exceed the local
    // minimum curvature radius on either support surface, or the
    // offset surface self-intersects before the intersection solver
    // has a chance. The self-intersection checker returns the
    // largest safe radius, which we thread into
    // `BlendFailure::RadiusTooLarge.max_radius` for the remediation
    // hint. §10.4 "feasibility bounds".
    if let Some(max_safe) = self_intersection_safe_radius(
        &analysis.surface_1,
        &analysis.surface_2,
        radius,
        analysis.orientation_1,
        analysis.orientation_2,
        u1_range,
        v1_range,
        u2_range,
        v2_range,
    ) {
        let requested = largest_requested_radius(radius);
        if requested > max_safe {
            return Err(BlendFailure::RadiusTooLarge {
                requested,
                max_radius: Some(max_safe),
                location: Some(edge_midpoint),
            });
        }
    }

    // Phase 1b: when either surface is a Cone, the offset-intersect
    // path can't run (cones have no analytic offset). Route the
    // entire center-curve computation through the spine walker, which
    // builds the center curve from rolling-ball contact solutions
    // along the actual edge curve.
    let either_is_cone = matches!(
        analysis.surface_1,
        popkern_core::geometry::surface::SurfaceGeometry::Cone { .. }
    ) || matches!(
        analysis.surface_2,
        popkern_core::geometry::surface::SurfaceGeometry::Cone { .. }
    );

    let center_result = if either_is_cone {
        crate::geometry_compute_cone::compute_cone_center_via_walker(
            analysis,
            radius,
            brep,
            u1_range,
            v1_range,
            u2_range,
            v2_range,
            tolerance,
            edge_midpoint,
        )?
    } else {
        match radius {
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
            .map_err(|e| BlendFailure::IntersectionFailed {
                reason: e,
                location: Some(edge_midpoint),
            })?,
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
            .map_err(|e| BlendFailure::IntersectionFailed {
                reason: e,
                location: Some(edge_midpoint),
            })?,
            FilletRadius::ByPoints(pts) => {
                if pts.is_empty() {
                    return Err(BlendFailure::InvalidTopology {
                        reason: "ByPoints radius with empty point list".to_string(),
                    });
                }
                // Round 5l: build the center curve using a piecewise-
                // linear closure over the *full* ByPoints profile, not
                // just the endpoint radii. The intermediate points now
                // shape the center curve rather than being silently
                // dropped.
                center_curve::parametric_radius_center_curve(
                    &analysis.surface_1,
                    u1_range,
                    v1_range,
                    &analysis.surface_2,
                    u2_range,
                    v2_range,
                    |t| center_curve::bypoints_local_radius(pts, t),
                    &edge_midpoint,
                    tolerance,
                    analysis.orientation_1,
                    analysis.orientation_2,
                )
                .map_err(|e| BlendFailure::IntersectionFailed {
                    reason: e,
                    location: Some(edge_midpoint),
                })?
            }
        }
    };

    let (effective_radius, radius_at_start, radius_at_end) = match radius {
        FilletRadius::Constant(r) => (*r, *r, *r),
        FilletRadius::Variable { start, end } => ((*start + *end) * 0.5, *start, *end),
        FilletRadius::ByPoints(pts) => {
            if pts.is_empty() {
                return Err(BlendFailure::InvalidTopology {
                    reason: "ByPoints radius with empty point list".to_string(),
                });
            }
            let avg = pts.iter().map(|(_, r)| *r).sum::<f64>() / pts.len() as f64;
            // For ByPoints, look up the radii nearest the spec's
            // start (parameter 0) and end (parameter 1) points.
            let r_at = |target: f64| -> f64 {
                pts.iter()
                    .min_by(|a, b| {
                        (a.0 - target)
                            .abs()
                            .partial_cmp(&(b.0 - target).abs())
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .map(|(_, r)| *r)
                    .unwrap_or(avg)
            };
            (avg, r_at(0.0), r_at(1.0))
        }
    };

    // Round 5k: a Variable or ByPoints spec must NOT be collapsed
    // to an analytic torus, because the torus has a single minor
    // radius and would silently average the per-sample radii. The
    // `is_constant_radius` flag forces `build_fillet_surface` to
    // build a true NURBS sweep that follows the per-sample
    // footprints from `variable_radius_center_curve`.
    let is_constant_radius = matches!(radius, FilletRadius::Constant(_));

    let fillet_result = fillet_surface::build_fillet_surface(
        &analysis.surface_1,
        &analysis.surface_2,
        &center_result,
        effective_radius,
        8, // arc samples
        is_constant_radius,
        g2_cross_section,
        analysis.orientation_1,
        analysis.orientation_2,
    )
    .map_err(|e| BlendFailure::IntersectionFailed {
        reason: e,
        location: Some(edge_midpoint),
    })?;

    let trim_curves = trim::compute_trim_curves(
        &analysis.surface_1,
        &center_result.footprint_1,
        &analysis.surface_2,
        &center_result.footprint_2,
    )
    .map_err(|e| BlendFailure::IntersectionFailed {
        reason: e,
        location: Some(edge_midpoint),
    })?;

    Ok(EdgeGeometry {
        analysis: rebuild_analysis(analysis, brep)?,
        surface: fillet_result.surface,
        trim_3d_1: trim_curves.curve_3d_1,
        trim_3d_2: trim_curves.curve_3d_2,
        pcurve_1: trim_curves.pcurve_1,
        pcurve_2: trim_curves.pcurve_2,
        effective_radius,
        radius_at_start,
        radius_at_end,
    })
}

/// Compute chamfer geometry for one edge.
pub(crate) fn compute_chamfer_geometry(
    analysis: &analysis::EdgeAnalysis,
    kind: &ChamferKind,
    brep: &BRep,
    tolerance: &Tolerance,
) -> Result<EdgeGeometry, BlendFailure> {
    let chamfer_result =
        chamfer_curve::compute_chamfer(analysis, brep, kind, tolerance).map_err(|e| {
            BlendFailure::IntersectionFailed {
                reason: e,
                location: None,
            }
        })?;

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
        // Chamfers are not parameterised by per-end radii; both
        // are set to `effective_radius` so the shared construct
        // pipeline sees a constant profile.
        radius_at_start: effective_radius,
        radius_at_end: effective_radius,
    })
}

/// Largest requested radius in a radius spec. Used as the "requested"
/// value in `BlendFailure::RadiusTooLarge`, so the user sees the
/// worst-case instance of their spec, not an arbitrary midpoint.
pub(crate) fn largest_requested_radius(radius: &FilletRadius) -> f64 {
    match radius {
        FilletRadius::Constant(r) => *r,
        FilletRadius::Variable { start, end } => start.max(*end),
        FilletRadius::ByPoints(pts) => pts
            .iter()
            .map(|(_, r)| *r)
            .fold(0.0_f64, |acc, r| acc.max(r)),
    }
}

/// Pre-flight self-intersection check for a rolling-ball fillet at the
/// requested radius. Returns `Some(max_safe_radius)` when the
/// requested radius would exceed the feasibility envelope on either
/// support surface, or `None` when the radius is safe (i.e. the
/// caller should proceed with the normal geometry pipeline).
///
/// For [`FilletRadius::Constant`] we test the single value. For
/// variable/by-points specs we test the peak requested radius —
/// if the peak is unsafe, the spec is unconditionally unsafe.
#[allow(clippy::too_many_arguments)]
pub(crate) fn self_intersection_safe_radius(
    surface_1: &popkern_core::geometry::surface::SurfaceGeometry,
    surface_2: &popkern_core::geometry::surface::SurfaceGeometry,
    radius: &FilletRadius,
    orientation_1: bool,
    orientation_2: bool,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    u2_range: (f64, f64),
    v2_range: (f64, f64),
) -> Option<f64> {
    use crate::offset::self_intersect::check_offset_self_intersection;

    let peak = largest_requested_radius(radius);
    if peak <= 0.0 {
        return None;
    }

    // "inward" matches the logic in center_curve.rs — offset_surface's
    // inward flag moves AGAINST the normal, and face normals point
    // into the solid, so we need `inward = !orientation`.
    let inward_1 = !orientation_1;
    let inward_2 = !orientation_2;

    let max_1 = check_offset_self_intersection(surface_1, peak, inward_1, u1_range, v1_range);
    let max_2 = check_offset_self_intersection(surface_2, peak, inward_2, u2_range, v2_range);

    // Both surfaces reported "safe" — nothing to flag.
    match (max_1, max_2) {
        (None, None) => None,
        // One or both reported a safety bound — take the tighter
        // (smaller) bound since both surfaces have to accommodate
        // the rolling ball simultaneously.
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
    }
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
#[path = "geometry_compute_tests.rs"]
mod tests;
