//! Transactional materialization of face-use trim authority.

use crate::geometry::fitting::fit_curve;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};
use crate::geometry::{CurveEval, Point3, SurfaceEval};
use crate::topology::{
    BRep, BrepFaceUseTrimAuthorityFindingKind, BrepFaceUseTrimAuthorityOptions,
    BrepFaceUseTrimAuthorityReport, analyze_brep_face_use_trim_authority,
};

/// Structured face-use trim authority repair report.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepFaceUseTrimAuthorityRepairReport {
    /// Diagnostics before repair.
    pub before: BrepFaceUseTrimAuthorityReport,
    /// Diagnostics after repair.
    pub after: BrepFaceUseTrimAuthorityReport,
    /// Number of face-use pcurves materialized.
    pub materialized_face_uses: usize,
}

/// Successful face-use trim authority materialization.
#[derive(Debug, Clone)]
pub struct BrepFaceUseTrimAuthorityRepairResult {
    /// Repaired clone of the input BRep.
    pub brep: BRep,
    /// Structured repair report.
    pub report: BrepFaceUseTrimAuthorityRepairReport,
}

/// Transactional face-use trim authority repair failure.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepFaceUseTrimAuthorityRepairFailure {
    /// Structured repair report from the rolled-back clone.
    pub report: BrepFaceUseTrimAuthorityRepairReport,
}

/// Materialize explicit face-use pcurve authority on a cloned BRep.
pub fn repair_brep_face_use_trim_authority_transactionally(
    brep: &BRep,
    options: BrepFaceUseTrimAuthorityOptions,
) -> Result<BrepFaceUseTrimAuthorityRepairResult, Box<BrepFaceUseTrimAuthorityRepairFailure>> {
    let before = analyze_brep_face_use_trim_authority(brep, options);
    if !before.has_blocking_findings() {
        return Ok(BrepFaceUseTrimAuthorityRepairResult {
            brep: brep.clone(),
            report: BrepFaceUseTrimAuthorityRepairReport {
                before: before.clone(),
                after: before,
                materialized_face_uses: 0,
            },
        });
    }
    let mut repaired = brep.clone();
    let mut materialized = 0;
    for finding in &before.findings {
        if finding.kind == BrepFaceUseTrimAuthorityFindingKind::MissingPcurve {
            continue;
        }
        if repaired
            .face_use_pcurve(finding.face, finding.half_edge)
            .is_some()
        {
            continue;
        }
        let Some(pcurve) =
            synthesize_face_use_pcurve(&repaired, finding.face, finding.half_edge, options)
        else {
            continue;
        };
        if repaired.set_face_use_pcurve(finding.face, finding.half_edge, pcurve) {
            materialized += 1;
        }
    }
    let after = analyze_brep_face_use_trim_authority(&repaired, options);
    let repair_report = BrepFaceUseTrimAuthorityRepairReport {
        before,
        after,
        materialized_face_uses: materialized,
    };
    if repair_report.after.has_blocking_findings() {
        return Err(Box::new(BrepFaceUseTrimAuthorityRepairFailure {
            report: repair_report,
        }));
    }
    Ok(BrepFaceUseTrimAuthorityRepairResult {
        brep: repaired,
        report: repair_report,
    })
}

fn synthesize_face_use_pcurve(
    brep: &BRep,
    face: crate::topology::FaceHandle,
    half_edge: crate::topology::HalfEdgeHandle,
    options: BrepFaceUseTrimAuthorityOptions,
) -> Option<Pcurve> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let curve = brep.edges.get(half_edge_data.edge)?.curve.as_ref()?;
    let forward = half_edge_curve_forward(brep, half_edge, curve)?;
    let samples = options.samples.max(33);
    let (period_u, period_v) = surface.periodicity();
    let mut prev_uv: Option<(f64, f64)> = None;
    let mut uv_points = Vec::with_capacity(samples);
    let mut params = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let curve_fraction = if forward { fraction } else { 1.0 - fraction };
        let t = curve.t_start + (curve.t_end - curve.t_start) * curve_fraction;
        let point = curve.geometry.point_at(t);
        let (mut u, mut v) = if let Some((prev_u, prev_v)) = prev_uv {
            project_point_to_uv_near(surface, &point, prev_u, prev_v)?
        } else {
            project_point_to_uv(surface, &point)?
        };
        unwrap_periodic_uv(&mut u, &mut v, prev_uv, period_u, period_v);
        prev_uv = Some((u, v));
        uv_points.push(Point3::new(u, v, 0.0));
        params.push(fraction);
    }
    let pcurve = if uv_points_are_linear(&uv_points, options.tolerance) {
        Pcurve::from_nurbs(linear_nurbs(*uv_points.first()?, *uv_points.last()?))
    } else {
        let degree = 3.min(uv_points.len().saturating_sub(1));
        Pcurve::from_nurbs(fit_curve(&uv_points, &params, degree, uv_points.len()).ok()?)
    };
    let gap = sampled_gap(brep, half_edge, face, &pcurve, samples)?;
    (gap <= options.tolerance).then_some(pcurve)
}

fn sampled_gap(
    brep: &BRep,
    half_edge: crate::topology::HalfEdgeHandle,
    face: crate::topology::FaceHandle,
    pcurve: &Pcurve,
    samples: usize,
) -> Option<f64> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let curve = brep.edges.get(half_edge_data.edge)?.curve.as_ref()?;
    let forward = half_edge_curve_forward(brep, half_edge, curve)?;
    let samples = samples.max(2);
    let mut max_gap = 0.0_f64;
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let curve_fraction = if forward { fraction } else { 1.0 - fraction };
        let curve_t = curve.t_start + (curve.t_end - curve.t_start) * curve_fraction;
        let pcurve_t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * fraction;
        let uv = pcurve.point_at(pcurve_t);
        max_gap = max_gap.max(
            curve
                .geometry
                .point_at(curve_t)
                .distance_to(&surface.point_at(uv.x, uv.y)),
        );
    }
    Some(max_gap)
}

fn half_edge_curve_forward(
    brep: &BRep,
    half_edge: crate::topology::HalfEdgeHandle,
    curve: &crate::geometry::CurveOnEdge,
) -> Option<bool> {
    let origin = brep
        .vertices
        .get(brep.half_edges.get(half_edge)?.origin)?
        .position;
    Some(
        origin.distance_to(&curve.geometry.point_at(curve.t_start))
            <= origin.distance_to(&curve.geometry.point_at(curve.t_end)),
    )
}

fn unwrap_periodic_uv(
    u: &mut f64,
    v: &mut f64,
    prev_uv: Option<(f64, f64)>,
    period_u: Option<f64>,
    period_v: Option<f64>,
) {
    let Some((prev_u, prev_v)) = prev_uv else {
        return;
    };
    if let Some(period) = period_u {
        while *u - prev_u > period * 0.5 {
            *u -= period;
        }
        while prev_u - *u > period * 0.5 {
            *u += period;
        }
    }
    if let Some(period) = period_v {
        while *v - prev_v > period * 0.5 {
            *v -= period;
        }
        while prev_v - *v > period * 0.5 {
            *v += period;
        }
    }
}

fn uv_points_are_linear(points: &[Point3], tolerance: f64) -> bool {
    let Some(start) = points.first().copied() else {
        return false;
    };
    let Some(end) = points.last().copied() else {
        return false;
    };
    let direction = end - start;
    let denom = direction.dot(&direction);
    if denom <= tolerance * tolerance {
        return false;
    }
    points.iter().enumerate().all(|(index, point)| {
        let fraction = index as f64 / points.len().saturating_sub(1).max(1) as f64;
        let expected = start + direction * fraction;
        point.distance_to(&expected) <= tolerance
    })
}

fn linear_nurbs(start: Point3, end: Point3) -> crate::geometry::NurbsCurve {
    crate::geometry::NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![start, end],
        weights: vec![1.0, 1.0],
    }
}

#[cfg(test)]
#[path = "brep_face_use_trim_authority_repair_tests.rs"]
mod tests;
