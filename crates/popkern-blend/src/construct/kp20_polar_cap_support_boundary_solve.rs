//! External cap/support boundary-row solve for committed KP20 cap candidates.

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use crate::construct::kp20_polar_support_cell_boundary_trim as boundary_trim;
pub(super) use boundary_trim::BoundarySide;
#[path = "kp20_polar_cap_support_boundary_fit.rs"]
mod boundary_fit;
#[path = "kp20_polar_cap_support_boundary_constraints.rs"]
mod constraints;
#[path = "kp20_polar_cap_support_edge_curve.rs"]
mod edge_curve;
#[path = "kp20_polar_cap_support_boundary_geometry.rs"]
mod geometry;
#[path = "kp20_polar_cap_support_boundary_solve_report.rs"]
mod report;
use crate::construct::kp20_polar_support_cell_row_fit as row_fit;
use crate::construct::kp20_polar_support_cell_surface_access as surface_access;

pub(super) use report::CapSupportBoundarySolveReport;

const SAMPLE_COUNT: usize = 12;
pub(super) const MAX_BOUNDARY_ROW_DISPLACEMENT: f64 = boundary_fit::MAX_BOUNDARY_ROW_DISPLACEMENT;

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
        report.max_samples_per_row = report.max_samples_per_row.max(samples.samples.len());
        if samples.samples.len() > SAMPLE_COUNT {
            report.multi_edge_rows += 1;
        }
        report.max_row_parameter_span = report
            .max_row_parameter_span
            .max(row_parameter_span(&samples.samples));
        report.max_row_before_sample_gap =
            report
                .max_row_before_sample_gap
                .max(boundary_fit::row_sample_gap(
                    &before,
                    samples.side,
                    &samples.samples,
                ));
        let hard_constraints = constraints::hard_constraints_for_row(samples, &mut report);
        let after = boundary_fit::best_boundary_fit_candidate(
            &before,
            samples.side,
            &samples.samples,
            &hard_constraints,
        );
        report.row_fit_candidates += after.candidate_count;
        report.constrained_row_fit_candidates += after.constrained_candidate_count;
        report.constrained_row_singular_candidates += after.constrained_singular_candidates;
        if after
            .min_hard_satisfying_constrained_displacement
            .is_finite()
        {
            let d = after.min_hard_satisfying_constrained_displacement;
            if d > report.max_hard_satisfying_constrained_displacement {
                report.max_hard_satisfying_constrained_displacement = d;
                report.worst_hard_satisfying_face = Some(samples.face);
                report.worst_hard_satisfying_side = Some(samples.side);
            }
            report.min_hard_satisfying_constrained_displacement =
                report.min_hard_satisfying_constrained_displacement.min(d);
            if d > MAX_BOUNDARY_ROW_DISPLACEMENT {
                report.hard_satisfying_rows_over_cap += 1;
            }
        }
        if after.best_unscaled_constrained_residual.is_finite() {
            report.max_best_unscaled_constrained_residual = report
                .max_best_unscaled_constrained_residual
                .max(after.best_unscaled_constrained_residual);
        }
        if after.hard_constraint_count > 0 {
            report.hard_constraint_rows += 1;
            report.hard_constraint_samples += after.hard_constraint_count;
        }
        report.max_hard_constraint_residual = report
            .max_hard_constraint_residual
            .max(after.hard_constraint_residual);
        report.displacement_scaled_candidates += after.displacement_scaled_candidates;
        report.displacement_rejected_candidates += after.displacement_rejected_candidates;
        if after.within_displacement {
            report.displacement_bounded_rows += 1;
        } else {
            report.displacement_unbounded_rows += 1;
        }
        if after.regularization == 0.0 {
            report.zero_regularization_rows += 1;
        }
        report.max_boundary_displacement = report.max_boundary_displacement.max(after.displacement);
        if after.sample_gap > report.max_row_after_sample_gap {
            report.max_row_after_sample_gap = after.sample_gap;
            report.max_row_after_sample_face = Some(samples.face);
            report.max_row_after_sample_side = Some(samples.side);
            report.max_row_after_sample_parameter = after.sample_parameter;
            if let Some(source_sample) =
                source_sample_near_parameter(&samples.source_samples, after.sample_parameter)
            {
                report.max_row_after_sample_edge = Some(source_sample.edge);
                report.max_row_after_sample_cap_half_edge = Some(source_sample.cap_half_edge);
                report.max_row_after_sample_support_face = Some(source_sample.support_face);
                report.max_row_after_sample_support_parameter = source_sample.support_parameter;
                report.max_row_after_sample_target_source = Some(source_sample.target_source);
                report.max_row_after_sample_before_gap = source_sample.before_gap;
                report.max_row_after_sample_target_edge_curve_gap =
                    source_sample.target_edge_curve_gap;
            }
        }
        for edge_curve_sample in &samples.edge_curve_samples {
            let after_gap =
                row_fit::eval_boundary(&after.patch, samples.side, edge_curve_sample.parameter)
                    .distance_to(&edge_curve_sample.point);
            if after_gap > report.max_row_after_edge_curve_gap {
                report.max_row_after_edge_curve_gap = after_gap;
                report.max_row_after_edge_curve_edge = Some(edge_curve_sample.edge);
                report.max_row_after_edge_curve_face = Some(samples.face);
                report.max_row_after_edge_curve_side = Some(samples.side);
                report.max_row_after_edge_curve_parameter = edge_curve_sample.parameter;
            }
        }
        surface_access::set_rectangular_binonic_on_face(brep, samples.face, after.patch)?;
        report.fit_rows += 1;
    }

    for edge_sample in &edge_samples {
        edge_curve::reauthor_from_samples(brep, edge_sample.cap_half_edge, &edge_sample.targets)?;
        report.reauthored_edge_curves += 1;
    }

    for edge_sample in &edge_samples {
        if let Some(after_gap) = external_cap_support_edge_gap(brep, edge_sample) {
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
    source_samples: Vec<BoundaryRowSourceSample>,
    edge_curve_samples: Vec<BoundaryRowEdgeCurveSample>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum BoundaryRowTargetSource {
    StartVertex,
    SupportPcurve,
    EndVertex,
}

#[derive(Debug, Clone, Copy)]
struct BoundaryRowSourceSample {
    parameter: f64,
    edge: EdgeHandle,
    cap_half_edge: HalfEdgeHandle,
    support_face: FaceHandle,
    support_parameter: f64,
    target_source: BoundaryRowTargetSource,
    target: Point3,
    before_gap: f64,
    target_edge_curve_gap: f64,
}

#[derive(Debug, Clone, Copy)]
struct BoundaryRowEdgeCurveSample {
    parameter: f64,
    point: Point3,
    edge: EdgeHandle,
}

#[derive(Debug, Clone)]
struct EdgeSolveSamples {
    cap_half_edge: HalfEdgeHandle,
    cap_side: BoundarySide,
    cap_span_start: f64,
    cap_span_end: f64,
    cap_forward: bool,
    targets: Vec<Point3>,
}

fn row_parameter_span(samples: &[row_fit::BoundaryRowFitSample]) -> f64 {
    let Some(first) = samples.first() else {
        return 0.0;
    };
    let mut min_parameter = first.parameter;
    let mut max_parameter = first.parameter;
    for sample in samples.iter().skip(1) {
        min_parameter = min_parameter.min(sample.parameter);
        max_parameter = max_parameter.max(sample.parameter);
    }
    max_parameter - min_parameter
}

fn source_sample_near_parameter(
    samples: &[BoundaryRowSourceSample],
    parameter: f64,
) -> Option<&BoundaryRowSourceSample> {
    samples.iter().min_by(|a, b| {
        (a.parameter - parameter)
            .abs()
            .partial_cmp(&(b.parameter - parameter).abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    })
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
    let cap_forward = geometry::cap_boundary_forward(
        &cap_patch,
        cap_side,
        cap_span_start,
        cap_span_end,
        cap_start,
        cap_end,
    );
    let support_forward =
        geometry::support_pcurve_forward(&support_surface, support_pcurve, cap_start, cap_end);
    let (support_t0, support_t1) = support_pcurve.parameter_range();
    let edge_curve = brep
        .edges
        .get(half_edge_data.edge)
        .and_then(|edge| edge.curve.as_ref());
    let edge_curve_forward = edge_curve
        .map(|curve| geometry::edge_curve_forward(curve, cap_start, cap_end))
        .unwrap_or(true);

    let mut samples = Vec::with_capacity(SAMPLE_COUNT);
    let mut source_samples = Vec::with_capacity(SAMPLE_COUNT);
    let mut edge_curve_samples = Vec::with_capacity(SAMPLE_COUNT);
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
        let cap_parameter = geometry::lerp_scalar(cap_span_start, cap_span_end, cap_fraction);
        let support_parameter = geometry::lerp_scalar(support_t0, support_t1, support_fraction);
        let mut target =
            geometry::support_surface_point(&support_surface, support_pcurve, support_parameter);
        let mut target_source = BoundaryRowTargetSource::SupportPcurve;
        if index == 0 {
            target = cap_start;
            target_source = BoundaryRowTargetSource::StartVertex;
        } else if index + 1 == SAMPLE_COUNT {
            target = cap_end;
            target_source = BoundaryRowTargetSource::EndVertex;
        }
        let mut target_edge_curve_gap = 0.0_f64;
        if let Some(edge_curve) = edge_curve {
            let edge_fraction = if edge_curve_forward {
                fraction
            } else {
                1.0 - fraction
            };
            let edge_parameter =
                geometry::lerp_scalar(edge_curve.t_start, edge_curve.t_end, edge_fraction);
            let edge_point = edge_curve.geometry.point_at(edge_parameter);
            target_edge_curve_gap = target.distance_to(&edge_point);
            if target_edge_curve_gap > report.max_row_target_edge_curve_gap {
                report.max_row_target_edge_curve_gap = target_edge_curve_gap;
                report.max_row_target_edge_curve_edge = Some(half_edge_data.edge);
                report.max_row_target_edge_curve_face = Some(cap_face);
                report.max_row_target_edge_curve_side = Some(cap_side);
                report.max_row_target_edge_curve_parameter = cap_parameter;
            }
            edge_curve_samples.push(BoundaryRowEdgeCurveSample {
                parameter: cap_parameter,
                point: edge_point,
                edge: half_edge_data.edge,
            });
        }
        let before_sample_gap =
            row_fit::eval_boundary(&cap_patch, cap_side, cap_parameter).distance_to(&target);
        before_gap = before_gap.max(before_sample_gap);
        samples.push(row_fit::BoundaryRowFitSample {
            parameter: cap_parameter,
            target,
        });
        source_samples.push(BoundaryRowSourceSample {
            parameter: cap_parameter,
            edge: half_edge_data.edge,
            cap_half_edge: half_edge,
            support_face,
            support_parameter,
            target_source,
            target,
            before_gap: before_sample_gap,
            target_edge_curve_gap,
        });
        targets.push(target);
    }

    let key = (cap_face, boundary_side_key(cap_side));
    let row = row_samples
        .entry(key)
        .or_insert_with(|| BoundaryRowSamples {
            face: cap_face,
            side: cap_side,
            samples: Vec::new(),
            source_samples: Vec::new(),
            edge_curve_samples: Vec::new(),
        });
    row.samples.extend(samples);
    row.source_samples.extend(source_samples);
    row.edge_curve_samples.extend(edge_curve_samples);
    edge_samples.push(EdgeSolveSamples {
        cap_half_edge: half_edge,
        cap_side,
        cap_span_start,
        cap_span_end,
        cap_forward,
        targets,
    });
    report.sampled_edges += 1;
    report.max_before_gap = report.max_before_gap.max(before_gap);
}

fn external_cap_support_edge_gap(brep: &BRep, edge_sample: &EdgeSolveSamples) -> Option<f64> {
    let half_edge_data = brep.half_edges.get(edge_sample.cap_half_edge)?;
    let cap_face = brep.loops.get(half_edge_data.loop_)?.face;
    let cap_patch = surface_access::rectangular_binonic_on_face(brep, cap_face).ok()?;
    let mut gap = 0.0_f64;
    for index in 0..SAMPLE_COUNT {
        let fraction = index as f64 / (SAMPLE_COUNT - 1) as f64;
        let cap_fraction = if edge_sample.cap_forward {
            fraction
        } else {
            1.0 - fraction
        };
        let cap_parameter = geometry::lerp_scalar(
            edge_sample.cap_span_start,
            edge_sample.cap_span_end,
            cap_fraction,
        );
        let target = edge_sample.targets.get(index)?;
        gap = gap.max(
            row_fit::eval_boundary(&cap_patch, edge_sample.cap_side, cap_parameter)
                .distance_to(target),
        );
    }
    Some(gap)
}

fn boundary_side_key(side: BoundarySide) -> u8 {
    match side {
        BoundarySide::UMin => 0,
        BoundarySide::UMax => 1,
        BoundarySide::VMin => 2,
        BoundarySide::VMax => 3,
    }
}
