//! Internal geometry computation helpers for fillet and chamfer operations.
//!
//! Extracted from `lib.rs` to keep source files within the 500-line limit.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::analysis;
use crate::chamfer_curve;
use crate::construct::EdgeGeometry;
use crate::error::BlendFailure;
use crate::fillet_surface;
use crate::geometry_compute_caps::build_edge_caps;
use crate::geometry_compute_feasibility::{
    largest_requested_radius, self_intersection_safe_radius,
};
use crate::trim;
use crate::types::{CapKind, ChamferKind, ContinuityLevel, FilletRadius};

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
    continuity: ContinuityLevel,
    cap_kind: CapKind,
) -> Result<EdgeGeometry, BlendFailure> {
    let wants_g2 = continuity == ContinuityLevel::G2;
    let edge_midpoint = get_edge_midpoint(brep, analysis)?;

    // §13.1 degeneracy pre-flight: collapsed edge, endpoint cusp,
    // degenerate support normal, near-tangent knife-edge. Fails
    // fast on clearly-degenerate input before any real geometry
    // work. Repair (subdivide / reparameterize / stitch) is a
    // separate, per-category follow-up.
    if let Err(kind) = crate::degeneracy::check_edge_degeneracy(brep, analysis, tolerance) {
        return Err(BlendFailure::InvalidTopology {
            reason: kind.describe(),
            location: Some(edge_midpoint),
        });
    }

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

    use popkern_core::geometry::surface::SurfaceGeometry;
    let either_needs_walker = matches!(
        analysis.surface_1,
        SurfaceGeometry::Cone { .. }
            | SurfaceGeometry::Torus { .. }
            | SurfaceGeometry::Cylinder { .. }
    ) || matches!(
        analysis.surface_2,
        SurfaceGeometry::Cone { .. }
            | SurfaceGeometry::Torus { .. }
            | SurfaceGeometry::Cylinder { .. }
    );
    // §7 continuity downgrade tracking — populates
    // `EdgeGeometry::continuity_downgrade` so `lib.rs` emits
    // `BlendWarning::ContinuityDowngraded`.
    let (center_result, mut downgrade_reason) =
        crate::geometry_compute_cone::compute_center_with_fallback(
            analysis,
            radius,
            brep,
            u1_range,
            v1_range,
            u2_range,
            v2_range,
            tolerance,
            edge_midpoint,
            either_needs_walker,
            wants_g2,
        )?;

    let (effective_radius, radius_at_start, radius_at_end) = match radius {
        FilletRadius::Constant(r) => (*r, *r, *r),
        FilletRadius::Variable { start, end } => ((*start + *end) * 0.5, *start, *end),
        FilletRadius::ByPoints(pts) if pts.is_empty() => {
            return Err(BlendFailure::InvalidTopology {
                reason: "ByPoints radius with empty point list".to_string(),
                location: Some(edge_midpoint),
            });
        }
        FilletRadius::ByPoints(pts) => {
            let avg = pts.iter().map(|(_, r)| *r).sum::<f64>() / pts.len() as f64;
            // Look up radius at normalized parameter `target`.
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
        continuity,
        analysis.orientation_1,
        analysis.orientation_2,
    )
    .map_err(|e| BlendFailure::IntersectionFailed {
        reason: e,
        location: Some(edge_midpoint),
    })?;

    // Surface-builder downgrade (G2 sweep errored) — walker fallback takes precedence.
    if downgrade_reason.is_none() {
        downgrade_reason = fillet_result.continuity_downgrade.clone();
    }

    // §13.2 post-construction self-intersection detection. The
    // pre-flight `self_intersection_safe_radius` check above catches
    // the "radius > support curvature radius" class of failures
    // before the output surface exists; this checks the output
    // itself for the residual cases (spindle torus, twisted NURBS
    // sweep, coincident non-neighbour UV samples). Fires on a clean
    // surface-type mismatch rather than a cryptic downstream
    // topology error.
    {
        use popkern_core::geometry::surface::SurfaceEval;
        let is_analytic = !matches!(
            fillet_result.surface,
            popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_)
        );
        let ((u_lo, u_hi), (v_lo, v_hi)) = fillet_result.surface.parameter_range();
        let fillet_u_range = (
            if u_lo.is_finite() { u_lo } else { 0.0_f64 },
            if u_hi.is_finite() { u_hi } else { 1.0_f64 },
        );
        let fillet_v_range = (
            if v_lo.is_finite() { v_lo } else { 0.0_f64 },
            if v_hi.is_finite() { v_hi } else { 1.0_f64 },
        );
        if !is_analytic
            && let Err(issue) = crate::surface_validate::validate_fillet_surface(
                &fillet_result.surface,
                fillet_u_range,
                fillet_v_range,
                tolerance,
            )
        {
            return Err(BlendFailure::SelfIntersectingBlend {
                reason: issue.describe(),
                location: Some(edge_midpoint),
            });
        }
    }

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

    // §18 telemetry: when the walker produced the center curve,
    // `center_result.contact_solutions` carries one
    // `ContactSolution` per spine sample, each recording its
    // Newton iteration count AND its line-search rejection count.
    // Sum both across the samples for per-edge totals. Legacy
    // offset-intersect paths leave `contact_solutions` as `None`,
    // which propagates as `None` in both fields.
    let (iterations_used, step_rejections_used) = center_result
        .contact_solutions
        .as_ref()
        .map(|sols| {
            let iters: usize = sols.iter().map(|s| s.iterations).sum();
            let rejs: usize = sols.iter().map(|s| s.step_rejections).sum();
            (Some(iters), Some(rejs))
        })
        .unwrap_or((None, None));

    // §8.3 caps: geometry at each blend endpoint (pure, no BRep mutation).
    let v_start_pos = brep.vertices[analysis.vertex_start].position;
    let v_end_pos = brep.vertices[analysis.vertex_end].position;
    let (start_cap, end_cap) = build_edge_caps(
        cap_kind,
        &center_result,
        radius_at_start,
        radius_at_end,
        v_start_pos,
        v_end_pos,
        &trim_curves.curve_3d_1,
        &trim_curves.curve_3d_2,
    );

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
        iterations_used,
        step_rejections_used,
        start_cap,
        end_cap,
        continuity_downgrade: crate::continuity::to_field(downgrade_reason),
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
                location: edge_midpoint_by_handle(brep, analysis.edge),
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
        // Chamfers don't run the walker/contact solver — always None.
        iterations_used: None,
        step_rejections_used: None,
        // Chamfers don't grow blend caps — the chamfer surface is
        // itself planar/ruled and terminates cleanly at the
        // perpendicular edge splits.
        start_cap: None,
        end_cap: None,
        // Chamfers are by definition G1 only — no continuity option.
        continuity_downgrade: None,
    })
}

// `largest_requested_radius` and `self_intersection_safe_radius`
// live in `geometry_compute_feasibility.rs` — extracted to keep
// this file under the 500-line structural cap. Re-exported via
// `pub(crate) use` at the top of the file.

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
            location: None,
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

/// Best-effort 3D midpoint lookup for an edge by handle. Used by
/// error-construction sites that want to populate
/// `BlendFailure::*.location` without paying the cost of analysing
/// the edge first. Returns `None` when the edge handle is stale or
/// the edge carries neither a curve nor usable endpoint vertices.
///
/// Prefer this helper over inlining the curve-midpoint-or-vertex-
/// average dance — it keeps every `location: ...` site uniform.
pub(crate) fn edge_midpoint_by_handle(
    brep: &BRep,
    edge: popkern_core::topology::EdgeHandle,
) -> Option<popkern_core::geometry::Point3> {
    let e = brep.edges.get(edge)?;
    if let Some(curve) = &e.curve {
        let t_mid = (curve.t_start + curve.t_end) * 0.5;
        return Some(curve.geometry.point_at(t_mid));
    }
    // Curve-less edges fall back to vertex midpoint via a half-edge walk.
    let he = brep.half_edges.get(e.half_edge)?;
    let origin = brep.vertices.get(he.origin)?.position;
    let next = brep.half_edges.get(he.next)?;
    let dest = brep.vertices.get(next.origin)?.position;
    Some(popkern_core::geometry::Point3::new(
        (origin.x + dest.x) * 0.5,
        (origin.y + dest.y) * 0.5,
        (origin.z + dest.z) * 0.5,
    ))
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
    popkern_core::geometry::uv_projection::face_bounded_parameter_range(brep, face_handle, surface)
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

/// Build overlap-detection inputs from computed edge geometries.
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
