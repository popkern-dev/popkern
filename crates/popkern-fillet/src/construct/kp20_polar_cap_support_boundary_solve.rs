//! External cap/support boundary-row solve for committed KP20 cap candidates.
#![allow(clippy::duplicate_mod)]

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::fitting::fit_curve;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

#[path = "kp20_polar_support_cell_boundary_trim.rs"]
#[allow(dead_code)]
mod boundary_trim;
pub(super) use boundary_trim::BoundarySide;
#[path = "kp20_polar_support_cell_row_fit.rs"]
#[allow(dead_code)]
mod row_fit;
#[path = "kp20_polar_support_cell_surface_access.rs"]
#[allow(dead_code)]
mod surface_access;

const SAMPLE_COUNT: usize = 12;

#[derive(Debug, Clone, Copy, Default)]
pub(super) struct CapSupportBoundarySolveReport {
    pub(super) eligible_edges: usize,
    pub(super) sampled_edges: usize,
    pub(super) fit_rows: usize,
    pub(super) reauthored_edge_curves: usize,
    pub(super) missing_cap_pcurve_edges: usize,
    pub(super) missing_support_pcurve_edges: usize,
    pub(super) unrecognized_cap_pcurve_edges: usize,
    pub(super) missing_support_surface_edges: usize,
    pub(super) non_rectangular_cap_edges: usize,
    pub(super) max_before_gap: f64,
    pub(super) max_after_gap: f64,
    pub(super) max_boundary_displacement: f64,
}

pub(super) fn solve_external_cap_support_boundaries(
    brep: &mut BRep,
    cap_faces: &BTreeSet<FaceHandle>,
) -> Result<CapSupportBoundarySolveReport, CoreError> {
    let mut report = CapSupportBoundarySolveReport::default();
    let mut row_samples: BTreeMap<(FaceHandle, u8), BoundaryRowSamples> = BTreeMap::new();
    let mut edge_samples = Vec::new();

    for &cap_face in cap_faces {
        let Ok(half_edges) = boundary_trim::face_outer_half_edges(brep, cap_face) else {
            continue;
        };
        for half_edge in half_edges {
            collect_external_cap_support_edge(
                brep,
                cap_faces,
                cap_face,
                half_edge,
                &mut row_samples,
                &mut edge_samples,
                &mut report,
            );
        }
    }

    for samples in row_samples.values() {
        let before = surface_access::rectangular_binonic_on_face(brep, samples.face)?;
        let mut after = before.clone();
        row_fit::fit_boundary_samples(&mut after, samples.side, &samples.samples);
        report.max_boundary_displacement =
            report
                .max_boundary_displacement
                .max(row_fit::max_boundary_displacement(
                    &before,
                    &after,
                    samples.side,
                ));
        surface_access::set_rectangular_binonic_on_face(brep, samples.face, after)?;
        report.fit_rows += 1;
    }

    for edge_sample in &edge_samples {
        reauthor_edge_curve_from_samples(brep, edge_sample.cap_half_edge, &edge_sample.targets)?;
        report.reauthored_edge_curves += 1;
    }

    for edge_sample in &edge_samples {
        if let Some(after_gap) = external_cap_support_edge_gap(brep, edge_sample.cap_half_edge) {
            report.max_after_gap = report.max_after_gap.max(after_gap);
        }
    }

    Ok(report)
}

#[derive(Debug, Clone)]
struct BoundaryRowSamples {
    face: FaceHandle,
    side: BoundarySide,
    samples: Vec<row_fit::BoundaryRowFitSample>,
}

#[derive(Debug, Clone)]
struct EdgeSolveSamples {
    cap_half_edge: HalfEdgeHandle,
    targets: Vec<Point3>,
}

fn collect_external_cap_support_edge(
    brep: &BRep,
    cap_faces: &BTreeSet<FaceHandle>,
    cap_face: FaceHandle,
    half_edge: HalfEdgeHandle,
    row_samples: &mut BTreeMap<(FaceHandle, u8), BoundaryRowSamples>,
    edge_samples: &mut Vec<EdgeSolveSamples>,
    report: &mut CapSupportBoundarySolveReport,
) {
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return;
    };
    let twin = half_edge_data.twin;
    let Some(twin_data) = brep.half_edges.get(twin) else {
        return;
    };
    let Some(support_face) = brep.loops.get(twin_data.loop_).map(|loop_| loop_.face) else {
        return;
    };
    if cap_faces.contains(&support_face) {
        return;
    }
    report.eligible_edges += 1;

    let Some(cap_pcurve) = brep.half_edges[half_edge].pcurve.as_ref() else {
        report.missing_cap_pcurve_edges += 1;
        return;
    };
    let Some(support_pcurve) = brep.half_edges[twin].pcurve.as_ref() else {
        report.missing_support_pcurve_edges += 1;
        return;
    };
    let Some((cap_side, cap_span_start, cap_span_end)) = boundary_trim::boundary_span(cap_pcurve)
    else {
        report.unrecognized_cap_pcurve_edges += 1;
        return;
    };
    let Ok(cap_patch) = surface_access::rectangular_binonic_on_face(brep, cap_face) else {
        report.non_rectangular_cap_edges += 1;
        return;
    };
    let Some(support_surface) = brep
        .faces
        .get(support_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| surface.geometry.clone())
    else {
        report.missing_support_surface_edges += 1;
        return;
    };

    let cap_start = brep.vertices[half_edge_data.origin].position;
    let cap_end = brep.vertices[twin_data.origin].position;
    let cap_forward = cap_boundary_forward(
        &cap_patch,
        cap_side,
        cap_span_start,
        cap_span_end,
        cap_start,
        cap_end,
    );
    let support_forward =
        support_pcurve_forward(&support_surface, support_pcurve, cap_start, cap_end);
    let (support_t0, support_t1) = support_pcurve.parameter_range();

    let mut samples = Vec::with_capacity(SAMPLE_COUNT);
    let mut targets = Vec::with_capacity(SAMPLE_COUNT);
    let mut before_gap = 0.0_f64;
    for index in 0..SAMPLE_COUNT {
        let fraction = index as f64 / (SAMPLE_COUNT - 1) as f64;
        let cap_fraction = if cap_forward {
            fraction
        } else {
            1.0 - fraction
        };
        let support_fraction = if support_forward {
            fraction
        } else {
            1.0 - fraction
        };
        let cap_parameter = lerp_scalar(cap_span_start, cap_span_end, cap_fraction);
        let support_parameter = lerp_scalar(support_t0, support_t1, support_fraction);
        let mut target = support_surface_point(&support_surface, support_pcurve, support_parameter);
        if index == 0 {
            target = cap_start;
        } else if index + 1 == SAMPLE_COUNT {
            target = cap_end;
        }
        before_gap = before_gap
            .max(row_fit::eval_boundary(&cap_patch, cap_side, cap_parameter).distance_to(&target));
        samples.push(row_fit::BoundaryRowFitSample {
            parameter: cap_parameter,
            target,
        });
        targets.push(target);
    }

    let key = (cap_face, boundary_side_key(cap_side));
    row_samples
        .entry(key)
        .or_insert_with(|| BoundaryRowSamples {
            face: cap_face,
            side: cap_side,
            samples: Vec::new(),
        })
        .samples
        .extend(samples);
    edge_samples.push(EdgeSolveSamples {
        cap_half_edge: half_edge,
        targets,
    });
    report.sampled_edges += 1;
    report.max_before_gap = report.max_before_gap.max(before_gap);
}

fn external_cap_support_edge_gap(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<f64> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let twin = half_edge_data.twin;
    let twin_data = brep.half_edges.get(twin)?;
    let cap_face = brep.loops.get(half_edge_data.loop_)?.face;
    let support_face = brep.loops.get(twin_data.loop_)?.face;
    let cap_pcurve = brep.half_edges[half_edge].pcurve.as_ref()?;
    let support_pcurve = brep.half_edges[twin].pcurve.as_ref()?;
    let (cap_side, cap_span_start, cap_span_end) = boundary_trim::boundary_span(cap_pcurve)?;
    let cap_patch = surface_access::rectangular_binonic_on_face(brep, cap_face).ok()?;
    let support_surface = brep
        .faces
        .get(support_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| surface.geometry.clone())?;
    let cap_start = brep.vertices[half_edge_data.origin].position;
    let cap_end = brep.vertices[twin_data.origin].position;
    let cap_forward = cap_boundary_forward(
        &cap_patch,
        cap_side,
        cap_span_start,
        cap_span_end,
        cap_start,
        cap_end,
    );
    let support_forward =
        support_pcurve_forward(&support_surface, support_pcurve, cap_start, cap_end);
    let (support_t0, support_t1) = support_pcurve.parameter_range();
    let mut gap = 0.0_f64;
    for index in 0..SAMPLE_COUNT {
        let fraction = index as f64 / (SAMPLE_COUNT - 1) as f64;
        let cap_fraction = if cap_forward {
            fraction
        } else {
            1.0 - fraction
        };
        let support_fraction = if support_forward {
            fraction
        } else {
            1.0 - fraction
        };
        let cap_parameter = lerp_scalar(cap_span_start, cap_span_end, cap_fraction);
        let support_parameter = lerp_scalar(support_t0, support_t1, support_fraction);
        gap = gap.max(
            row_fit::eval_boundary(&cap_patch, cap_side, cap_parameter).distance_to(
                &support_surface_point(&support_surface, support_pcurve, support_parameter),
            ),
        );
    }
    Some(gap)
}

fn cap_boundary_forward(
    cap_patch: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    side: BoundarySide,
    start: f64,
    end: f64,
    half_edge_start: Point3,
    half_edge_end: Point3,
) -> bool {
    let p0 = row_fit::eval_boundary(cap_patch, side, start);
    let p1 = row_fit::eval_boundary(cap_patch, side, end);
    p0.distance_to(&half_edge_start) + p1.distance_to(&half_edge_end)
        <= p1.distance_to(&half_edge_start) + p0.distance_to(&half_edge_end)
}

fn support_pcurve_forward(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    half_edge_start: Point3,
    half_edge_end: Point3,
) -> bool {
    let (t0, t1) = pcurve.parameter_range();
    let p0 = support_surface_point(surface, pcurve, t0);
    let p1 = support_surface_point(surface, pcurve, t1);
    p0.distance_to(&half_edge_start) + p1.distance_to(&half_edge_end)
        <= p1.distance_to(&half_edge_start) + p0.distance_to(&half_edge_end)
}

fn support_surface_point(surface: &SurfaceGeometry, pcurve: &Pcurve, t: f64) -> Point3 {
    let uv = pcurve.point_at(t);
    surface.point_at(uv.x, uv.y)
}

fn reauthor_edge_curve_from_samples(
    brep: &mut BRep,
    cap_half_edge: HalfEdgeHandle,
    samples: &[Point3],
) -> Result<(), CoreError> {
    let Some(half_edge) = brep.half_edges.get(cap_half_edge) else {
        return Err(CoreError::InvalidHandle(
            "cap/support boundary solve half-edge invalid".into(),
        ));
    };
    let edge = half_edge.edge;
    let mut edge_samples = samples.to_vec();
    if brep.edges[edge].half_edge != cap_half_edge {
        edge_samples.reverse();
    }
    let curve = fit_nurbs_curve(&edge_samples)
        .map(|geometry| CurveOnEdge {
            geometry: CurveGeometry::Nurbs(geometry),
            t_start: 0.0,
            t_end: 1.0,
        })
        .or_else(|| line_curve_from_endpoints(brep, &edge_samples))
        .ok_or_else(|| {
            CoreError::PreconditionFailed(
                "cap/support boundary solve found degenerate edge samples".into(),
            )
        })?;
    brep.edges[edge].curve = Some(curve);
    Ok(())
}

fn fit_nurbs_curve(points: &[Point3]) -> Option<NurbsCurve> {
    if points.len() < 4 {
        return None;
    }
    let params = (0..points.len())
        .map(|index| index as f64 / (points.len() - 1) as f64)
        .collect::<Vec<_>>();
    fit_curve(points, &params, 3, points.len()).ok()
}

fn line_curve_from_endpoints(brep: &BRep, points: &[Point3]) -> Option<CurveOnEdge> {
    let start = *points.first()?;
    let end = *points.last()?;
    CurveOnEdge::line_segment(start, end, brep.tolerance.linear().max(1.0e-12))
}

fn boundary_side_key(side: BoundarySide) -> u8 {
    match side {
        BoundarySide::UMin => 0,
        BoundarySide::UMax => 1,
        BoundarySide::VMin => 2,
        BoundarySide::VMax => 3,
    }
}

fn lerp_scalar(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}
