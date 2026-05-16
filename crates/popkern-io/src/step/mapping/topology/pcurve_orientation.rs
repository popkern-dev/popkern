use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceOnFace};
use popkern_core::tolerance::Tolerance;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::types::EntityId;

pub(super) fn select_orient_pcurve_to_half_edge(
    pcurves: Vec<Pcurve>,
    curve: &CurveOnEdge,
    surface: &SurfaceOnFace,
    edge_forward: bool,
    face_id: EntityId,
    tolerance: Tolerance,
    log: &mut DiagnosticLog,
) -> Option<Pcurve> {
    let mut best: Option<(Pcurve, f64)> = None;
    for pcurve in pcurves {
        let candidate = align_pcurve_candidate(pcurve, curve, surface, edge_forward, tolerance);
        match &best {
            Some((_, best_gap)) if candidate.1 >= *best_gap => {}
            _ => best = Some(candidate),
        }
    }

    let (pcurve, accepted_gap) = best?;
    if accepted_gap > tolerance.linear() {
        log.warn(
            Phase::MapTopology,
            Some(face_id),
            format!(
                "STEP pcurve endpoint disagreement {:.3e} exceeds tolerance {:.3e}",
                accepted_gap,
                tolerance.linear()
            ),
        );
    }
    Some(pcurve)
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

    let (_, direct_gap) = pcurve_endpoint_gap(&direct, surface, curve_start, curve_end);
    let (_, reversed_gap) = pcurve_endpoint_gap(&reversed, surface, curve_start, curve_end);

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

    #[test]
    fn select_orient_pcurve_aligns_endpoint_images_to_half_edge() {
        let surface = xy_plane();
        let curve = unit_x_edge();
        let wrong = Pcurve::linear((0.0, 1.0), (1.0, 1.0), 0.0, 1.0);
        let reversed = Pcurve::linear((1.0, 0.0), (0.0, 0.0), 0.0, 1.0);
        let mut log = DiagnosticLog::new();

        let selected = select_orient_pcurve_to_half_edge(
            vec![wrong, reversed],
            &curve,
            &surface,
            true,
            42,
            Tolerance::default(),
            &mut log,
        )
        .expect("one pcurve candidate should be selected");

        let (start, end) = pcurve_endpoint_images(&selected, &surface);
        assert!(start.distance_to(&Point3::new(0.0, 0.0, 0.0)) <= Tolerance::default().linear());
        assert!(end.distance_to(&Point3::new(1.0, 0.0, 0.0)) <= Tolerance::default().linear());
    }
}
