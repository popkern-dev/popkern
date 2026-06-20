//! Canonical cap-side pcurves from cap boundary evidence.

use popkern_core::geometry::Point3;
use popkern_core::geometry::edge_triple::BoundarySide;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, HalfEdgeHandle};

use super::support_pcurve_completion_quality::{
    oriented_pcurve_sample_gap, pcurve_surface_sample_gap,
};

pub(super) const MAX_CAP_BOUNDARY_SAMPLE_GAP: f64 = 1.0e-1;

const MAX_INFERRED_CAP_BOUNDARY_UV_DEVIATION: f64 = 5.0e-2;
const PARAMETER_OVERSHOOT_TOL: f64 = 1.0e-3;
const BOUNDARY_COORDINATE_TOL: f64 = 1.0e-9;
const MIN_BOUNDARY_SPAN: f64 = 1.0e-9;

#[derive(Debug, Clone)]
pub(super) struct CapBoundaryPcurveCandidate {
    pub(super) boundary: BoundarySide,
    pub(super) pcurve: Pcurve,
    pub(super) source: CapBoundaryPcurveSource,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum CapBoundaryPcurveSource {
    Hint,
    Inferred,
}

pub(super) fn cap_boundary_pcurve_candidate(
    boundary_hint: Option<BoundarySide>,
    existing_cap_pcurve: Option<&Pcurve>,
    uv_samples: &[Point3],
) -> Option<CapBoundaryPcurveCandidate> {
    if let Some(boundary) = boundary_hint {
        return hinted_cap_boundary_pcurve(boundary, uv_samples).map(|pcurve| {
            CapBoundaryPcurveCandidate {
                boundary,
                pcurve,
                source: CapBoundaryPcurveSource::Hint,
            }
        });
    }
    if let Some((boundary, pcurve)) =
        existing_cap_pcurve.and_then(inferred_cap_boundary_pcurve_from_existing)
    {
        return Some(CapBoundaryPcurveCandidate {
            boundary,
            pcurve,
            source: CapBoundaryPcurveSource::Inferred,
        });
    }
    inferred_cap_boundary_pcurve(uv_samples).map(|(boundary, pcurve)| CapBoundaryPcurveCandidate {
        boundary,
        pcurve,
        source: CapBoundaryPcurveSource::Inferred,
    })
}

pub(super) fn cap_boundary_pcurve_candidate_is_acceptable(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    cap_surface: &SurfaceGeometry,
    candidate: &CapBoundaryPcurveCandidate,
    samples: &[Point3],
) -> bool {
    let candidate_gap = pcurve_surface_sample_gap(cap_surface, &candidate.pcurve, samples);
    if !candidate_gap.is_finite() || candidate_gap > MAX_CAP_BOUNDARY_SAMPLE_GAP {
        return false;
    }
    let Some(existing) = brep
        .half_edges
        .get(half_edge)
        .and_then(|he| he.pcurve.as_ref())
    else {
        return true;
    };
    let existing_gap = oriented_pcurve_sample_gap(brep, half_edge, cap_surface, existing, samples);
    if pcurve_matches_boundary(existing, candidate.boundary) {
        candidate_gap + 1.0e-12 < existing_gap
    } else {
        candidate_gap <= existing_gap + MAX_CAP_BOUNDARY_SAMPLE_GAP
    }
}

pub(super) fn existing_cap_boundary_pcurve_is_acceptable(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    cap_surface: &SurfaceGeometry,
    boundary: BoundarySide,
    pcurve: &Pcurve,
    samples: &[Point3],
) -> bool {
    pcurve_matches_boundary(pcurve, boundary)
        && oriented_pcurve_sample_gap(brep, half_edge, cap_surface, pcurve, samples)
            <= MAX_CAP_BOUNDARY_SAMPLE_GAP
}

fn hinted_cap_boundary_pcurve(boundary: BoundarySide, uv_samples: &[Point3]) -> Option<Pcurve> {
    let start = uv_samples.first()?;
    let end = uv_samples.last()?;
    let start_t = unit_boundary_parameter(boundary, *start)?;
    let end_t = unit_boundary_parameter(boundary, *end)?;
    if (start_t - end_t).abs() <= MIN_BOUNDARY_SPAN {
        return None;
    }
    let start_t = start_t.clamp(0.0, 1.0);
    let end_t = end_t.clamp(0.0, 1.0);
    match boundary {
        BoundarySide::UMin => Some(Pcurve::linear((0.0, start_t), (0.0, end_t), 0.0, 1.0)),
        BoundarySide::UMax => Some(Pcurve::linear((1.0, start_t), (1.0, end_t), 0.0, 1.0)),
        BoundarySide::VMin => Some(Pcurve::linear((start_t, 0.0), (end_t, 0.0), 0.0, 1.0)),
        BoundarySide::VMax => Some(Pcurve::linear((start_t, 1.0), (end_t, 1.0), 0.0, 1.0)),
        BoundarySide::Indexed(_) => None,
    }
}

fn inferred_cap_boundary_pcurve(uv_samples: &[Point3]) -> Option<(BoundarySide, Pcurve)> {
    let boundary = [
        BoundarySide::UMin,
        BoundarySide::UMax,
        BoundarySide::VMin,
        BoundarySide::VMax,
    ]
    .into_iter()
    .filter_map(|boundary| {
        boundary_uv_deviation(boundary, uv_samples)
            .filter(|deviation| *deviation <= MAX_INFERRED_CAP_BOUNDARY_UV_DEVIATION)
            .map(|deviation| (boundary, deviation))
    })
    .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))?
    .0;
    hinted_cap_boundary_pcurve(boundary, uv_samples).map(|pcurve| (boundary, pcurve))
}

fn inferred_cap_boundary_pcurve_from_existing(pcurve: &Pcurve) -> Option<(BoundarySide, Pcurve)> {
    let (t0, t1) = pcurve.parameter_range();
    let uv_samples = (0..=4)
        .map(|index| {
            let fraction = index as f64 / 4.0;
            pcurve.point_at(t0 + fraction * (t1 - t0))
        })
        .collect::<Vec<_>>();
    inferred_cap_boundary_pcurve(&uv_samples)
}

fn pcurve_matches_boundary(pcurve: &Pcurve, boundary: BoundarySide) -> bool {
    pcurve_boundary_side(pcurve) == Some(boundary)
}

fn boundary_uv_deviation(boundary: BoundarySide, uv_samples: &[Point3]) -> Option<f64> {
    uv_samples
        .iter()
        .map(|uv| match boundary {
            BoundarySide::UMin => Some(uv.x.abs()),
            BoundarySide::UMax => Some((uv.x - 1.0).abs()),
            BoundarySide::VMin => Some(uv.y.abs()),
            BoundarySide::VMax => Some((uv.y - 1.0).abs()),
            BoundarySide::Indexed(_) => None,
        })
        .try_fold(0.0_f64, |max_deviation, deviation| {
            let deviation = deviation?;
            deviation
                .is_finite()
                .then_some(max_deviation.max(deviation))
        })
}

fn unit_boundary_parameter(boundary: BoundarySide, uv: Point3) -> Option<f64> {
    let parameter = match boundary {
        BoundarySide::UMin | BoundarySide::UMax => uv.y,
        BoundarySide::VMin | BoundarySide::VMax => uv.x,
        BoundarySide::Indexed(_) => return None,
    };
    (parameter.is_finite()
        && (-PARAMETER_OVERSHOOT_TOL..=1.0 + PARAMETER_OVERSHOOT_TOL).contains(&parameter))
    .then_some(parameter)
}

fn pcurve_boundary_side(pcurve: &Pcurve) -> Option<BoundarySide> {
    let (t0, t1) = pcurve.parameter_range();
    let start = pcurve.point_at(t0);
    let end = pcurve.point_at(t1);
    if close(start.x, 0.0) && close(end.x, 0.0) && usable_span(start.y, end.y) {
        Some(BoundarySide::UMin)
    } else if close(start.x, 1.0) && close(end.x, 1.0) && usable_span(start.y, end.y) {
        Some(BoundarySide::UMax)
    } else if close(start.y, 0.0) && close(end.y, 0.0) && usable_span(start.x, end.x) {
        Some(BoundarySide::VMin)
    } else if close(start.y, 1.0) && close(end.y, 1.0) && usable_span(start.x, end.x) {
        Some(BoundarySide::VMax)
    } else {
        None
    }
}

fn usable_span(a: f64, b: f64) -> bool {
    (a - b).abs() > MIN_BOUNDARY_SPAN && (0.0..=1.0).contains(&a) && (0.0..=1.0).contains(&b)
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= BOUNDARY_COORDINATE_TOL
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hinted_cap_boundary_pcurve_authors_exact_unit_side() {
        let pcurve = cap_boundary_pcurve_candidate(
            Some(BoundarySide::UMax),
            None,
            &[Point3::new(0.97, 0.2, 0.0), Point3::new(0.94, 0.8, 0.0)],
        )
        .unwrap()
        .pcurve;
        assert!(pcurve_matches_boundary(&pcurve, BoundarySide::UMax));
        assert_eq!(pcurve.point_at(0.0).x, 1.0);
        assert_eq!(pcurve.point_at(1.0).x, 1.0);
    }

    #[test]
    fn inferred_cap_boundary_pcurve_snaps_near_unit_side() {
        let candidate = cap_boundary_pcurve_candidate(
            None,
            None,
            &[Point3::new(0.98, 0.2, 0.0), Point3::new(0.97, 0.8, 0.0)],
        )
        .unwrap();
        assert_eq!(candidate.source, CapBoundaryPcurveSource::Inferred);
        assert_eq!(candidate.boundary, BoundarySide::UMax);
        assert!(pcurve_matches_boundary(
            &candidate.pcurve,
            BoundarySide::UMax
        ));
    }

    #[test]
    fn inferred_cap_boundary_pcurve_rejects_interior_trace() {
        assert!(
            cap_boundary_pcurve_candidate(
                None,
                None,
                &[Point3::new(0.75, 0.2, 0.0), Point3::new(0.72, 0.8, 0.0)],
            )
            .is_none()
        );
    }

    #[test]
    fn inferred_cap_boundary_pcurve_snaps_existing_trace() {
        let existing = Pcurve::linear((0.98, 0.2), (0.97, 0.8), 0.0, 1.0);
        let candidate = cap_boundary_pcurve_candidate(
            None,
            Some(&existing),
            &[Point3::new(0.4, 0.4, 0.0), Point3::new(0.6, 0.6, 0.0)],
        )
        .unwrap();
        assert_eq!(candidate.source, CapBoundaryPcurveSource::Inferred);
        assert_eq!(candidate.boundary, BoundarySide::UMax);
        assert!(pcurve_matches_boundary(
            &candidate.pcurve,
            BoundarySide::UMax
        ));
    }
}
