use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceOnFace};
use popkern_core::tolerance::Tolerance;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::mapping::geometry::ResolvedStepPcurve;
use crate::step::parser::types::EntityId;

use crate::step::same_parameter::same_parameter_repair;

#[derive(Debug, Clone)]
pub(super) struct SelectedStepPcurve {
    pub pcurve: Pcurve,
    pub pcurve_id: EntityId,
    pub reference_curve_id: EntityId,
    pub sampled_gap: f64,
    /// True when the file's pcurve failed the same-fraction check and
    /// was rebuilt on the edge's parameterization (PK-552).
    pub same_parameter_rebuilt: bool,
}

pub(super) fn select_orient_pcurve_to_half_edge(
    pcurves: Vec<ResolvedStepPcurve>,
    curve: &CurveOnEdge,
    surface: &SurfaceOnFace,
    edge_forward: bool,
    face_id: EntityId,
    tolerance: Tolerance,
    log: &mut DiagnosticLog,
) -> Option<SelectedStepPcurve> {
    let mut best: Option<(ResolvedStepPcurve, Pcurve, f64)> = None;
    for resolved in &pcurves {
        let candidate = align_pcurve_candidate(
            resolved.pcurve.clone(),
            curve,
            surface,
            edge_forward,
            tolerance,
        );
        match &best {
            Some((_, _, best_score)) if candidate.1 >= *best_score => {}
            _ => best = Some((resolved.clone(), candidate.0, candidate.1)),
        }
    }

    let (resolved, pcurve, accepted_score) = best?;
    if accepted_score > tolerance.linear()
        && let Some((rebuilt, repair)) = pcurves
            .iter()
            .filter_map(|r| {
                same_parameter_repair(&r.pcurve, curve, surface, edge_forward, tolerance.linear())
                    .map(|repair| (r, repair))
            })
            .min_by(|a, b| a.1.max_gap.total_cmp(&b.1.max_gap))
    {
        return Some(SelectedStepPcurve {
            pcurve: repair.pcurve,
            pcurve_id: rebuilt.pcurve_id,
            reference_curve_id: rebuilt.reference_curve_id,
            sampled_gap: repair.max_gap,
            same_parameter_rebuilt: true,
        });
    }
    if accepted_score > tolerance.linear() {
        log.warn(
            Phase::MapTopology,
            Some(resolved.pcurve_id),
            format!(
                "STEP pcurve sampled disagreement {:.3e} exceeds tolerance {:.3e} on ADVANCED_FACE #{face_id}",
                accepted_score, tolerance.linear()
            ),
        );
    }
    Some(SelectedStepPcurve {
        pcurve,
        pcurve_id: resolved.pcurve_id,
        reference_curve_id: resolved.reference_curve_id,
        sampled_gap: accepted_score,
        same_parameter_rebuilt: false,
    })
}

fn align_pcurve_candidate(
    pcurve: Pcurve,
    curve: &CurveOnEdge,
    surface: &SurfaceOnFace,
    edge_forward: bool,
    tolerance: Tolerance,
) -> (Pcurve, f64) {
    let (curve_start, curve_end) = curve_endpoints_for_half_edge(curve, edge_forward);
    let direct = pcurve.clone();
    let reversed = reverse_pcurve_candidate(pcurve);
    let _ = tolerance;

    let direct_gap = pcurve_candidate_gap(
        &direct,
        curve,
        surface,
        edge_forward,
        curve_start,
        curve_end,
    );
    let reversed_gap = pcurve_candidate_gap(
        &reversed,
        curve,
        surface,
        edge_forward,
        curve_start,
        curve_end,
    );

    if reversed_gap < direct_gap {
        (reversed, reversed_gap)
    } else {
        (direct, direct_gap)
    }
}

fn reverse_pcurve_candidate(mut pcurve: Pcurve) -> Pcurve {
    match &mut pcurve.geometry {
        PcurveGeometry::Linear { start, end } => {
            std::mem::swap(start, end);
        }
        PcurveGeometry::Elliptical {
            start_angle,
            end_angle,
            ..
        } => {
            std::mem::swap(start_angle, end_angle);
        }
        PcurveGeometry::Sinusoidal { .. }
        | PcurveGeometry::Nurbs(_)
        | PcurveGeometry::CylinderArcViaRationalQuad { .. } => {
            std::mem::swap(&mut pcurve.t_start, &mut pcurve.t_end);
        }
    }
    pcurve
}

fn curve_endpoints_for_half_edge(
    curve: &CurveOnEdge,
    edge_forward: bool,
) -> (
    popkern_core::geometry::Point3,
    popkern_core::geometry::Point3,
) {
    if edge_forward {
        (
            curve.geometry.point_at(curve.t_start),
            curve.geometry.point_at(curve.t_end),
        )
    } else {
        (
            curve.geometry.point_at(curve.t_end),
            curve.geometry.point_at(curve.t_start),
        )
    }
}

fn pcurve_endpoint_images(
    pcurve: &Pcurve,
    surface: &SurfaceOnFace,
) -> (
    popkern_core::geometry::Point3,
    popkern_core::geometry::Point3,
) {
    let uv_start = pcurve.point_at(pcurve.t_start);
    let uv_end = pcurve.point_at(pcurve.t_end);
    (
        surface.geometry.point_at(uv_start.x, uv_start.y),
        surface.geometry.point_at(uv_end.x, uv_end.y),
    )
}

fn pcurve_endpoint_gap(
    pcurve: &Pcurve,
    surface: &SurfaceOnFace,
    target_start: Point3,
    target_end: Point3,
) -> (
    (
        popkern_core::geometry::Point3,
        popkern_core::geometry::Point3,
    ),
    f64,
) {
    let (pc_start, pc_end) = pcurve_endpoint_images(pcurve, surface);
    let gap = target_start
        .distance_to(&pc_start)
        .max(target_end.distance_to(&pc_end));
    ((pc_start, pc_end), gap)
}

fn pcurve_candidate_gap(
    pcurve: &Pcurve,
    curve: &CurveOnEdge,
    surface: &SurfaceOnFace,
    edge_forward: bool,
    target_start: Point3,
    target_end: Point3,
) -> f64 {
    let (_, endpoint_gap) = pcurve_endpoint_gap(pcurve, surface, target_start, target_end);
    endpoint_gap.max(pcurve_same_fraction_gap(
        pcurve,
        curve,
        surface,
        edge_forward,
        9,
    ))
}

fn pcurve_same_fraction_gap(
    pcurve: &Pcurve,
    curve: &CurveOnEdge,
    surface: &SurfaceOnFace,
    edge_forward: bool,
    samples: usize,
) -> f64 {
    let edge_start = if edge_forward {
        curve.t_start
    } else {
        curve.t_end
    };
    let edge_end = if edge_forward {
        curve.t_end
    } else {
        curve.t_start
    };
    let mut max_gap = 0.0_f64;
    for i in 0..samples.max(2) {
        let fraction = i as f64 / (samples.max(2) as f64 - 1.0);
        let edge_t = edge_start + (edge_end - edge_start) * fraction;
        let pcurve_t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * fraction;
        let edge_point = curve.geometry.point_at(edge_t);
        let uv = pcurve.point_at(pcurve_t);
        let surface_point = surface.geometry.point_at(uv.x, uv.y);
        max_gap = max_gap.max(edge_point.distance_to(&surface_point));
    }
    max_gap
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::curve::CurveGeometry;
    use popkern_core::geometry::pcurve::Pcurve;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::geometry::{Point3, Vec3};

    fn xy_plane() -> SurfaceOnFace {
        SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
            orientation: true,
        }
    }

    fn unit_x_edge() -> CurveOnEdge {
        CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: Point3::new(0.0, 0.0, 0.0),
                direction: Vec3::new(1.0, 0.0, 0.0),
            },
            t_start: 0.0,
            t_end: 1.0,
        }
    }

    fn resolved(id: EntityId, pcurve: Pcurve) -> ResolvedStepPcurve {
        ResolvedStepPcurve {
            pcurve_id: id,
            reference_curve_id: id + 1000,
            pcurve,
        }
    }

    #[test]
    fn select_orient_pcurve_aligns_endpoint_images_to_half_edge() {
        let surface = xy_plane();
        let curve = unit_x_edge();
        let wrong = Pcurve::linear((0.0, 1.0), (1.0, 1.0), 0.0, 1.0);
        let reversed = Pcurve::linear((1.0, 0.0), (0.0, 0.0), 0.0, 1.0);
        let mut log = DiagnosticLog::new();

        let selected = select_orient_pcurve_to_half_edge(
            vec![resolved(100, wrong), resolved(101, reversed)],
            &curve,
            &surface,
            true,
            42,
            Tolerance::default(),
            &mut log,
        )
        .expect("one pcurve candidate should be selected");

        assert_eq!(selected.pcurve_id, 101);
        let (start, end) = pcurve_endpoint_images(&selected.pcurve, &surface);
        assert!(start.distance_to(&Point3::new(0.0, 0.0, 0.0)) <= Tolerance::default().linear());
        assert!(end.distance_to(&Point3::new(1.0, 0.0, 0.0)) <= Tolerance::default().linear());
    }

    #[test]
    fn select_orient_pcurve_prefers_same_fraction_trace_over_endpoint_tie() {
        let surface = xy_plane();
        let curve = unit_x_edge();
        let wrong = Pcurve::from_nurbs(popkern_core::geometry::NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        });
        let correct = Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0);
        let mut log = DiagnosticLog::new();

        let selected = select_orient_pcurve_to_half_edge(
            vec![resolved(100, wrong), resolved(101, correct)],
            &curve,
            &surface,
            true,
            42,
            Tolerance::default(),
            &mut log,
        )
        .expect("one pcurve candidate should be selected");

        assert!(
            pcurve_same_fraction_gap(&selected.pcurve, &curve, &surface, true, 9)
                <= Tolerance::default().linear()
        );
        assert_eq!(selected.pcurve_id, 101);
    }
}
