#![warn(missing_docs)]
//! popkern-fillet: Fillets and chamfers for the popkern BRep kernel.
//!
//! Provides [`fillet_edges`] and [`chamfer_edges`] for applying fillets and
//! chamfers to BRep edges, plus convenience wrappers [`fillet_edge`] and
//! [`chamfer_edge`] for single-edge operations.

pub mod analysis;
pub mod blend;
pub mod center_curve;
pub mod center_curve_helpers;
pub mod chamfer_curve;
pub mod construct;
pub mod error;
pub mod fillet_surface;
pub mod offset;
pub mod overlap;
pub mod trim;
pub mod types;

pub use error::{BlendError, BlendFailure, BlendWarning};
pub use types::*;

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

use crate::construct::EdgeGeometry;

/// Fillet the specified edges of a BRep solid.
///
/// Analyzes each edge, computes the rolling-ball center curve and fillet
/// surface, derives trim curves, and applies the topology modifications.
///
/// # Errors
///
/// Returns [`BlendError::NoValidEdges`] if analysis fails for all edges, or
/// [`BlendError::PartialFailure`] if some edges fail (with a partial result).
pub fn fillet_edges(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
) -> Result<BlendResult, BlendError> {
    if specs.is_empty() {
        return Err(BlendError::NoValidEdges("no fillet specs provided".into()));
    }

    // 1. Analyze edges.
    let edges: Vec<EdgeHandle> = specs.iter().map(|s| s.edge).collect();
    let analyses = analysis::analyze_edges(brep, &edges, options.propagate_tangent)
        .map_err(BlendError::NoValidEdges)?;

    // 2. Compute geometry for each edge.
    let mut edge_geometries = Vec::new();
    let mut failures = Vec::new();

    for (i, spec) in specs.iter().enumerate() {
        if i >= analyses.len() {
            continue;
        }
        match compute_fillet_geometry(&analyses[i], &spec.radius, brep, &tolerance) {
            Ok(geom) => edge_geometries.push(geom),
            Err(failure) => failures.push((spec.edge, failure)),
        }
    }

    if edge_geometries.is_empty() && !failures.is_empty() {
        return Err(BlendError::PartialFailure {
            failed: failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures,
        });
    }

    // 3. Apply topology.
    let apply_result =
        construct::apply_blend_geometry(brep, &edge_geometries, &tolerance, &options.vertex_mode)
            .map_err(|edge_failures| BlendError::PartialFailure {
            failed: edge_failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures: edge_failures,
        })?;

    Ok(BlendResult {
        brep: apply_result.brep,
        edge_warnings: apply_result.warnings,
    })
}

/// Chamfer the specified edges of a BRep solid.
///
/// Analyzes each edge, computes offset curves on each adjacent face, builds
/// a ruled surface between them, and applies topology modifications.
///
/// # Errors
///
/// Returns [`BlendError::NoValidEdges`] if analysis fails for all edges, or
/// [`BlendError::PartialFailure`] if some edges fail (with a partial result).
pub fn chamfer_edges(
    brep: &BRep,
    specs: &[ChamferSpec],
    options: &ChamferOptions,
    tolerance: Tolerance,
) -> Result<BlendResult, BlendError> {
    if specs.is_empty() {
        return Err(BlendError::NoValidEdges("no chamfer specs provided".into()));
    }

    // 1. Analyze edges.
    let edges: Vec<EdgeHandle> = specs.iter().map(|s| s.edge).collect();
    let analyses = analysis::analyze_edges(brep, &edges, options.propagate_tangent)
        .map_err(BlendError::NoValidEdges)?;

    // 2. Compute geometry for each edge.
    let mut edge_geometries = Vec::new();
    let mut failures = Vec::new();

    for (i, spec) in specs.iter().enumerate() {
        if i >= analyses.len() {
            continue;
        }
        match compute_chamfer_geometry(&analyses[i], &spec.kind, brep, &tolerance) {
            Ok(geom) => edge_geometries.push(geom),
            Err(failure) => failures.push((spec.edge, failure)),
        }
    }

    if edge_geometries.is_empty() && !failures.is_empty() {
        return Err(BlendError::PartialFailure {
            failed: failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures,
        });
    }

    // 3. Apply topology. Chamfers always use mitred vertex mode.
    let apply_result =
        construct::apply_blend_geometry(brep, &edge_geometries, &tolerance, &VertexMode::Mitred)
            .map_err(|edge_failures| BlendError::PartialFailure {
                failed: edge_failures.len(),
                total: specs.len(),
                partial_result: Box::new(brep.clone()),
                failures: edge_failures,
            })?;

    Ok(BlendResult {
        brep: apply_result.brep,
        edge_warnings: apply_result.warnings,
    })
}

/// Convenience: fillet a single edge with a constant radius.
///
/// Equivalent to calling [`fillet_edges`] with a single [`FilletSpec`] and
/// default options.
///
/// # Errors
///
/// Returns a [`BlendError`] if the edge cannot be filleted.
pub fn fillet_edge(
    brep: &BRep,
    edge: EdgeHandle,
    radius: f64,
    tolerance: Tolerance,
) -> Result<BlendResult, BlendError> {
    fillet_edges(
        brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        }],
        &FilletOptions::default(),
        tolerance,
    )
}

/// Convenience: chamfer a single edge with equal distances.
///
/// Equivalent to calling [`chamfer_edges`] with a single [`ChamferSpec`] and
/// default options.
///
/// # Errors
///
/// Returns a [`BlendError`] if the edge cannot be chamfered.
pub fn chamfer_edge(
    brep: &BRep,
    edge: EdgeHandle,
    distance: f64,
    tolerance: Tolerance,
) -> Result<BlendResult, BlendError> {
    chamfer_edges(
        brep,
        &[ChamferSpec {
            edge,
            kind: ChamferKind::Equal(distance),
        }],
        &ChamferOptions::default(),
        tolerance,
    )
}

// ---------------------------------------------------------------------------
// Internal geometry computation
// ---------------------------------------------------------------------------

/// Compute fillet geometry (center curve, fillet surface, trim curves) for one edge.
fn compute_fillet_geometry(
    analysis: &analysis::EdgeAnalysis,
    radius: &FilletRadius,
    brep: &BRep,
    tolerance: &Tolerance,
) -> Result<EdgeGeometry, BlendFailure> {
    // Determine edge midpoint for center-curve selection.
    let edge_midpoint = get_edge_midpoint(brep, analysis)?;

    // Get surface parameter ranges.
    let (u1_range, v1_range) = surface_param_range(&analysis.surface_1);
    let (u2_range, v2_range) = surface_param_range(&analysis.surface_2);

    // Compute the center curve.
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

    // Build fillet surface.
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

    // Compute trim curves (pcurves in UV space).
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
fn compute_chamfer_geometry(
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
fn get_edge_midpoint(
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
        // Fallback: midpoint of the two endpoint vertices.
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
fn surface_param_range(
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
fn rebuild_analysis(
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
