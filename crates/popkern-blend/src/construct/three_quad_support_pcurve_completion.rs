//! Cap/support seam completion for Gregory faces after surface assignment.
use super::helpers;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};
use std::collections::BTreeSet;
#[path = "three_quad_support_pcurve_sampling.rs"]
mod support_pcurve_sampling;
use support_pcurve_sampling::{
    chord_samples_projected_to_support, samples_from_curve_for_half_edge,
    samples_from_support_pcurve, support_pcurve_endpoints_match, support_pcurve_matches_edge_curve,
};
#[path = "three_quad_support_pcurve_completion_quality.rs"]
mod support_pcurve_completion_quality;
use support_pcurve_completion_quality::{
    oriented_pcurve_sample_gap, support_pcurve_candidate_improves_samples,
};
#[path = "three_quad_support_pcurve_completion_fit.rs"]
mod support_pcurve_completion_fit;
use support_pcurve_completion_fit::{
    best_pcurve_for_surface_samples, fit_nurbs_curve, line_geometry_from_endpoints,
};
#[path = "three_quad_support_pcurve_seam_census.rs"]
mod seam_census;
#[path = "three_quad_support_pcurve_completion_cap_hint.rs"]
mod support_pcurve_completion_cap_hint;
use support_pcurve_completion_cap_hint::{
    CapBoundaryPcurveSource, cap_boundary_pcurve_candidate,
    cap_boundary_pcurve_candidate_is_acceptable, existing_cap_boundary_pcurve_is_acceptable,
};

const SAMPLE_COUNT: usize = 12;
const SUPPORT_PCURVE_ENDPOINT_AUTHORITY_TOL: f64 = 1.0e-5;

#[derive(Clone, Copy)]
struct CapSupportCompletionPolicy(bool, bool, bool);

const REPAIR_STALE_CAP_SUPPORT_SEAMS: CapSupportCompletionPolicy =
    CapSupportCompletionPolicy(true, true, true);

pub(crate) fn complete_missing_cap_support_seams(brep: &mut BRep, cap_face: FaceHandle) {
    let _ = complete_missing_cap_support_seams_when(
        brep,
        cap_face,
        |_| true,
        REPAIR_STALE_CAP_SUPPORT_SEAMS,
    );
}

pub(crate) fn complete_missing_external_cap_support_pcurves(
    brep: &mut BRep,
    cap_face: FaceHandle,
    cap_faces: &BTreeSet<FaceHandle>,
) -> CapSupportPcurveCompletionReport {
    seam_census::trace_cap_support_seam_census(brep, cap_face, cap_faces);
    complete_missing_cap_support_seams_when(
        brep,
        cap_face,
        |support_face| !cap_faces.contains(&support_face),
        REPAIR_STALE_CAP_SUPPORT_SEAMS,
    )
}

#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct CapSupportPcurveCompletionReport {
    pub(crate) eligible_edges: usize,
    pub(crate) authored_cap_boundary_pcurves: usize,
    pub(crate) repaired_cap_boundary_pcurves: usize,
    pub(crate) inferred_cap_boundary_pcurves: usize,
    pub(crate) declined_cap_boundary_pcurves: usize,
    pub(crate) authored_support_pcurves: usize,
    pub(crate) repaired_support_pcurves: usize,
}

#[derive(Debug, Clone, Copy, Default)]
struct CapSupportHalfEdgeCompletion {
    authored_cap_boundary_pcurve: bool,
    repaired_cap_boundary_pcurve: bool,
    inferred_cap_boundary_pcurve: bool,
    declined_cap_boundary_pcurve: bool,
    authored_support_pcurve: bool,
    repaired_support_pcurve: bool,
}

fn complete_missing_cap_support_seams_when<F>(
    brep: &mut BRep,
    cap_face: FaceHandle,
    include_support_face: F,
    policy: CapSupportCompletionPolicy,
) -> CapSupportPcurveCompletionReport
where
    F: Fn(FaceHandle) -> bool + Copy,
{
    let Some(cap_surface) = brep
        .faces
        .get(cap_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| surface.geometry.clone())
    else {
        return CapSupportPcurveCompletionReport::default();
    };
    if policy.0
        && !matches!(
            cap_surface,
            SurfaceGeometry::Gregory(_) | SurfaceGeometry::Nurbs(_)
        )
    {
        return CapSupportPcurveCompletionReport::default();
    }
    let mut report = CapSupportPcurveCompletionReport::default();
    for half_edge in cap_face_half_edges(brep, cap_face) {
        let Some(support_face) = twin_face_for_half_edge(brep, half_edge) else {
            continue;
        };
        if include_support_face(support_face) {
            report.eligible_edges += 1;
            let completion = complete_cap_support_half_edge(brep, &cap_surface, half_edge, policy);
            report.authored_cap_boundary_pcurves +=
                completion.authored_cap_boundary_pcurve as usize;
            report.repaired_cap_boundary_pcurves +=
                completion.repaired_cap_boundary_pcurve as usize;
            report.inferred_cap_boundary_pcurves +=
                completion.inferred_cap_boundary_pcurve as usize;
            report.declined_cap_boundary_pcurves +=
                completion.declined_cap_boundary_pcurve as usize;
            report.authored_support_pcurves += completion.authored_support_pcurve as usize;
            report.repaired_support_pcurves += completion.repaired_support_pcurve as usize;
        }
    }
    report
}

fn complete_cap_support_half_edge(
    brep: &mut BRep,
    cap_surface: &SurfaceGeometry,
    half_edge: HalfEdgeHandle,
    policy: CapSupportCompletionPolicy,
) -> CapSupportHalfEdgeCompletion {
    let Some(half_edge_data) = brep.half_edges.get(half_edge).cloned() else {
        return CapSupportHalfEdgeCompletion::default();
    };
    let twin = half_edge_data.twin;
    let Some(twin_data) = brep.half_edges.get(twin).cloned() else {
        return CapSupportHalfEdgeCompletion::default();
    };
    let Some(support_face) = brep
        .loops
        .get(twin_data.loop_)
        .map(|loop_data| loop_data.face)
    else {
        return CapSupportHalfEdgeCompletion::default();
    };
    let Some(support_surface) = brep
        .faces
        .get(support_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| surface.geometry.clone())
    else {
        return CapSupportHalfEdgeCompletion::default();
    };
    let edge = half_edge_data.edge;
    let edge_curve = brep
        .edges
        .get(edge)
        .and_then(|edge| edge.curve.as_ref())
        .cloned();
    let support_pcurve = brep.half_edges.get(twin).and_then(|he| he.pcurve.as_ref());
    let support_owns_seam = matches!(support_surface, SurfaceGeometry::Cylinder { .. })
        && support_pcurve.is_some_and(|pcurve| {
            support_pcurve_endpoints_match(
                brep,
                half_edge,
                twin,
                &support_surface,
                pcurve,
                SUPPORT_PCURVE_ENDPOINT_AUTHORITY_TOL,
            ) && edge_curve.as_ref().is_none_or(|curve| {
                support_pcurve_matches_edge_curve(
                    brep,
                    half_edge,
                    twin,
                    &support_surface,
                    pcurve,
                    curve,
                    SUPPORT_PCURVE_ENDPOINT_AUTHORITY_TOL,
                )
            })
        });
    let Some(mut samples) = cap_oriented_samples(
        brep,
        half_edge,
        twin,
        cap_surface,
        &support_surface,
        support_owns_seam,
    ) else {
        return CapSupportHalfEdgeCompletion::default();
    };
    if !samples.is_empty() {
        samples[0] = brep.vertices[brep.half_edges[half_edge].origin].position;
        let last = samples.len() - 1;
        samples[last] = brep.vertices[brep.half_edges[twin].origin].position;
    }
    let mut completion = CapSupportHalfEdgeCompletion::default();
    let mut support_samples = samples.clone();
    support_samples.reverse();
    let needs_edge_curve = edge_curve.is_none()
        || (policy.1
            && edge_curve.as_ref().is_some_and(|curve| {
                support_owns_seam
                    && edge_curve_sample_gap(brep, half_edge, curve, &samples) > 1.0e-5
            }));
    let cap_uv_samples = samples
        .iter()
        .map(|point| cap_uv_point(cap_surface, point))
        .collect::<Vec<_>>();
    let cap_boundary_hint = brep.surface_boundary_hint(half_edge);
    let cap_boundary_candidate = cap_boundary_pcurve_candidate(
        cap_boundary_hint,
        brep.half_edges[half_edge].pcurve.as_ref(),
        &cap_uv_samples,
        cap_surface,
        &samples,
    );
    let accepted_cap_boundary_candidate = if let Some(candidate) = cap_boundary_candidate {
        let acceptable = cap_boundary_pcurve_candidate_is_acceptable(
            brep,
            half_edge,
            cap_surface,
            &candidate,
            &samples,
        );
        if crate::direct_apex_trace::enabled() {
            let candidate_gap = support_pcurve_completion_quality::pcurve_surface_sample_gap(
                cap_surface,
                &candidate.pcurve,
                &samples,
            );
            eprintln!(
                "[cap-hint-trace] he {half_edge:?} boundary={:?} source={:?} acceptable={acceptable} candidate_gap={candidate_gap:.4e} samples={}",
                candidate.boundary,
                candidate.source,
                samples.len(),
            );
        }
        completion.declined_cap_boundary_pcurve = !acceptable;
        acceptable.then_some(candidate)
    } else if cap_boundary_hint.is_some() {
        completion.declined_cap_boundary_pcurve = true;
        None
    } else {
        None
    };
    let needs_cap_pcurve = if accepted_cap_boundary_candidate.is_some() {
        true
    } else {
        match brep.half_edges[half_edge].pcurve.as_ref() {
            None => true,
            Some(pcurve)
                if cap_boundary_hint.is_some_and(|boundary| {
                    existing_cap_boundary_pcurve_is_acceptable(
                        brep,
                        half_edge,
                        cap_surface,
                        boundary,
                        pcurve,
                        &samples,
                    )
                }) =>
            {
                false
            }
            Some(pcurve) if policy.2 => {
                oriented_pcurve_sample_gap(brep, half_edge, cap_surface, pcurve, &samples) > 1.0e-5
            }
            Some(_) => false,
        }
    };
    let needs_support_pcurve = match brep.half_edges[twin].pcurve.as_ref() {
        None => true,
        Some(_) if !policy.2 => false,
        Some(pcurve) => {
            oriented_pcurve_sample_gap(brep, twin, &support_surface, pcurve, &support_samples)
                > 1.0e-5
        }
    };
    if !needs_cap_pcurve && !needs_support_pcurve && !needs_edge_curve {
        return CapSupportHalfEdgeCompletion::default();
    }
    if needs_edge_curve {
        let mut edge_samples = samples.clone();
        if brep
            .edges
            .get(edge)
            .is_some_and(|edge_data| edge_data.half_edge != half_edge)
        {
            edge_samples.reverse();
        }
        let geometry = fit_nurbs_curve(&edge_samples)
            .map(CurveGeometry::Nurbs)
            .unwrap_or_else(|| line_geometry_from_endpoints(&edge_samples));
        brep.edges[edge].curve = Some(CurveOnEdge {
            geometry,
            t_start: 0.0,
            t_end: 1.0,
        });
    }
    if needs_cap_pcurve {
        let authored_cap_pcurve = brep.half_edges[half_edge].pcurve.is_none();
        brep.half_edges[half_edge].pcurve = if let Some(candidate) = accepted_cap_boundary_candidate
        {
            completion.authored_cap_boundary_pcurve = authored_cap_pcurve;
            completion.repaired_cap_boundary_pcurve = !authored_cap_pcurve;
            completion.inferred_cap_boundary_pcurve =
                candidate.source == CapBoundaryPcurveSource::Inferred;
            Some(candidate.pcurve)
        } else {
            completion.authored_cap_boundary_pcurve = authored_cap_pcurve;
            completion.repaired_cap_boundary_pcurve = !authored_cap_pcurve;
            Some(best_pcurve_for_surface_samples(
                cap_surface,
                &cap_uv_samples,
                &samples,
            ))
        };
    }
    if needs_support_pcurve {
        let uv_samples = support_samples
            .iter()
            .map(|point| unclamped_uv_point(&support_surface, point))
            .collect::<Vec<_>>();
        let candidate =
            best_pcurve_for_surface_samples(&support_surface, &uv_samples, &support_samples);
        if support_pcurve_candidate_improves_samples(
            brep,
            twin,
            &support_surface,
            &candidate,
            &support_samples,
        ) {
            let authored_support_pcurve = brep.half_edges[twin].pcurve.is_none();
            brep.half_edges[twin].pcurve = Some(candidate);
            completion.authored_support_pcurve = authored_support_pcurve;
            completion.repaired_support_pcurve = !authored_support_pcurve;
            return CapSupportHalfEdgeCompletion {
                authored_support_pcurve: completion.authored_support_pcurve,
                repaired_support_pcurve: completion.repaired_support_pcurve,
                ..completion
            };
        }
    }
    completion
}

fn cap_oriented_samples(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    twin: HalfEdgeHandle,
    cap_surface: &SurfaceGeometry,
    support_surface: &SurfaceGeometry,
    prefer_support_pcurve: bool,
) -> Option<Vec<Point3>> {
    if prefer_support_pcurve
        && let Some(support_pcurve) = brep.half_edges.get(twin)?.pcurve.as_ref()
    {
        return Some(samples_from_support_pcurve(
            brep,
            half_edge,
            twin,
            support_surface,
            support_pcurve,
        ));
    }
    if let Some(curve) = brep
        .half_edges
        .get(half_edge)
        .and_then(|half_edge| brep.edges.get(half_edge.edge))
        .and_then(|edge| edge.curve.as_ref())
    {
        return Some(samples_from_curve_for_half_edge(brep, half_edge, curve));
    }
    if let Some(cap_pcurve) = brep.half_edges.get(half_edge)?.pcurve.as_ref() {
        return Some(samples_from_support_pcurve(
            brep,
            half_edge,
            twin,
            cap_surface,
            cap_pcurve,
        ));
    }
    if let Some(support_pcurve) = brep.half_edges.get(twin)?.pcurve.as_ref() {
        return Some(samples_from_support_pcurve(
            brep,
            half_edge,
            twin,
            support_surface,
            support_pcurve,
        ));
    }
    let chord_samples = chord_samples_projected_to_support(brep, half_edge, twin, support_surface);
    if chord_samples.is_some() && crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace] cap seam samples from support-projected chord: {half_edge:?}"
        );
    }
    chord_samples
}

fn edge_curve_sample_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
    target_samples: &[Point3],
) -> f64 {
    samples_from_curve_for_half_edge(brep, half_edge, curve)
        .iter()
        .zip(target_samples)
        .map(|(sample, target)| sample.distance_to(target))
        .fold(0.0, f64::max)
}

fn cap_uv_point(surface: &SurfaceGeometry, point: &Point3) -> Point3 {
    let (mut u, mut v) = helpers::project_to_uv(surface, point);
    let ((u_lo, u_hi), (v_lo, v_hi)) = surface.parameter_range();
    if u_lo.is_finite() && u_hi.is_finite() {
        u = u.clamp(u_lo, u_hi);
    }
    if v_lo.is_finite() && v_hi.is_finite() {
        v = v.clamp(v_lo, v_hi);
    }
    Point3::new(u, v, 0.0)
}

fn unclamped_uv_point(surface: &SurfaceGeometry, point: &Point3) -> Point3 {
    let (u, v) = helpers::project_to_uv(surface, point);
    Point3::new(u, v, 0.0)
}

fn twin_face_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<FaceHandle> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let twin_data = brep.half_edges.get(half_edge_data.twin)?;
    brep.loops
        .get(twin_data.loop_)
        .map(|loop_data| loop_data.face)
}

fn cap_face_half_edges(brep: &BRep, cap_face: FaceHandle) -> Vec<HalfEdgeHandle> {
    let mut out = Vec::new();
    for (_, edge_data) in &brep.edges {
        let first = edge_data.half_edge;
        let Some(first_data) = brep.half_edges.get(first) else {
            continue;
        };
        let second = first_data.twin;
        let Some(second_data) = brep.half_edges.get(second) else {
            continue;
        };
        let first_face = brep
            .loops
            .get(first_data.loop_)
            .map(|loop_data| loop_data.face);
        let second_face = brep
            .loops
            .get(second_data.loop_)
            .map(|loop_data| loop_data.face);
        match (first_face, second_face) {
            (Some(face), Some(other)) if face == cap_face && other != cap_face => out.push(first),
            (Some(other), Some(face)) if face == cap_face && other != cap_face => out.push(second),
            _ => {}
        };
    }
    out
}

#[cfg(test)]
#[path = "three_quad_support_pcurve_completion_tests.rs"]
mod tests;
