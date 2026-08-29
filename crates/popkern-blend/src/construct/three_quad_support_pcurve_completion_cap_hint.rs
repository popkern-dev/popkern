//! Canonical cap-side pcurves from cap boundary evidence.

use popkern_core::geometry::edge_triple::BoundarySide;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, SurfaceEval};
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
    cap_surface: &SurfaceGeometry,
    samples: &[Point3],
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
    if let Some((boundary, pcurve)) = inferred_cap_boundary_pcurve(uv_samples) {
        return Some(CapBoundaryPcurveCandidate {
            boundary,
            pcurve,
            source: CapBoundaryPcurveSource::Inferred,
        });
    }
    inferred_cap_boundary_pcurve_from_surface_fit(cap_surface, samples).map(|(boundary, pcurve)| {
        CapBoundaryPcurveCandidate {
            boundary,
            pcurve,
            source: CapBoundaryPcurveSource::Inferred,
        }
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

/// Border inference from 3D evidence when UV projections are too noisy
/// to recognize the border. On small caps a within-gate 3D boundary
/// offset can exceed any fixed UV threshold (a committed KP15 cap at
/// scale ~2e-2 with a 2.6e-3 boundary residual projects ~0.12 charts
/// off), and nearest-point projections can jump borders near chart
/// corners, leaving diagonal off-border pcurves. This path derives
/// each candidate border's span by 1-D projection of the seam's
/// endpoint samples onto the surface border curve, then selects the
/// border by the 3D sample gap of the resulting exact Linear border
/// pcurve — the same metric the downstream acceptance enforces, so
/// the inference is scale-honest and projection-independent.
fn inferred_cap_boundary_pcurve_from_surface_fit(
    cap_surface: &SurfaceGeometry,
    samples: &[Point3],
) -> Option<(BoundarySide, Pcurve)> {
    let first = samples.first()?;
    let last = samples.last()?;
    let mut best: Option<(f64, BoundarySide, Pcurve)> = None;
    for boundary in [
        BoundarySide::UMin,
        BoundarySide::UMax,
        BoundarySide::VMin,
        BoundarySide::VMax,
    ] {
        let start_t = border_parameter_of_nearest_point(cap_surface, boundary, first);
        let end_t = border_parameter_of_nearest_point(cap_surface, boundary, last);
        if (start_t - end_t).abs() <= MIN_BOUNDARY_SPAN {
            continue;
        }
        let Some(pcurve) = border_linear_pcurve(boundary, start_t, end_t) else {
            continue;
        };
        let gap = pcurve_surface_sample_gap(cap_surface, &pcurve, samples);
        if !gap.is_finite() || gap > MAX_CAP_BOUNDARY_SAMPLE_GAP {
            continue;
        }
        if best.as_ref().is_none_or(|(best_gap, _, _)| gap < *best_gap) {
            best = Some((gap, boundary, pcurve));
        }
    }
    best.map(|(_, boundary, pcurve)| (boundary, pcurve))
}

/// Parameter along the unit border curve of `surface` whose 3D point
/// is nearest to `point` — coarse 17-point grid refined twice around
/// the running best (deterministic, derivative-free).
fn border_parameter_of_nearest_point(
    surface: &SurfaceGeometry,
    boundary: BoundarySide,
    point: &Point3,
) -> f64 {
    let eval = |t: f64| -> f64 {
        let (u, v) = match boundary {
            BoundarySide::UMin => (0.0, t),
            BoundarySide::UMax => (1.0, t),
            BoundarySide::VMin => (t, 0.0),
            BoundarySide::VMax => (t, 1.0),
            BoundarySide::Indexed(_) => return f64::INFINITY,
        };
        surface.point_at(u, v).distance_to(point)
    };
    let mut lo = 0.0_f64;
    let mut hi = 1.0_f64;
    let mut best_t = 0.0_f64;
    for _ in 0..3 {
        let mut best_d = f64::INFINITY;
        for k in 0..=16 {
            let t = lo + (hi - lo) * k as f64 / 16.0;
            let d = eval(t);
            if d < best_d {
                best_d = d;
                best_t = t;
            }
        }
        let span = (hi - lo) / 16.0;
        lo = (best_t - span).max(0.0);
        hi = (best_t + span).min(1.0);
    }
    best_t
}

/// Exact Linear pcurve along one canonical unit-chart border.
fn border_linear_pcurve(boundary: BoundarySide, start_t: f64, end_t: f64) -> Option<Pcurve> {
    match boundary {
        BoundarySide::UMin => Some(Pcurve::linear((0.0, start_t), (0.0, end_t), 0.0, 1.0)),
        BoundarySide::UMax => Some(Pcurve::linear((1.0, start_t), (1.0, end_t), 0.0, 1.0)),
        BoundarySide::VMin => Some(Pcurve::linear((start_t, 0.0), (end_t, 0.0), 0.0, 1.0)),
        BoundarySide::VMax => Some(Pcurve::linear((start_t, 1.0), (end_t, 1.0), 0.0, 1.0)),
        BoundarySide::Indexed(_) => None,
    }
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
    use popkern_core::geometry::Vec3;

    fn unit_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    #[test]
    fn hinted_cap_boundary_pcurve_authors_exact_unit_side() {
        let pcurve = cap_boundary_pcurve_candidate(
            Some(BoundarySide::UMax),
            None,
            &[Point3::new(0.97, 0.2, 0.0), Point3::new(0.94, 0.8, 0.0)],
            &unit_plane(),
            &[],
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
            &unit_plane(),
            &[],
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
                &unit_plane(),
                &[
                    Point3::new(0.75, 0.2, 0.0),
                    Point3::new(0.735, 0.5, 0.0),
                    Point3::new(0.72, 0.8, 0.0)
                ],
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
            &unit_plane(),
            &[],
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
    fn surface_fit_recovers_border_when_uv_projections_are_noisy() {
        // 3D samples lie exactly on the plane's u=0 border between
        // v=0.2 and v=0.9, but the supplied UV projections are
        // corner-jumped diagonal garbage (the committed-cap failure
        // mode) — the surface-fit path must still author the exact
        // UMin border span.
        let samples = [
            Point3::new(0.0, 0.2, 0.0),
            Point3::new(0.0, 0.55, 0.0),
            Point3::new(0.0, 0.9, 0.0),
        ];
        let noisy_uv = [
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.3, 0.5, 0.0),
            Point3::new(0.557, 0.0, 0.0),
        ];
        let candidate =
            cap_boundary_pcurve_candidate(None, None, &noisy_uv, &unit_plane(), &samples).unwrap();
        assert_eq!(candidate.source, CapBoundaryPcurveSource::Inferred);
        assert_eq!(candidate.boundary, BoundarySide::UMin);
        let start = candidate.pcurve.point_at(0.0);
        let end = candidate.pcurve.point_at(1.0);
        assert_eq!(start.x, 0.0);
        assert_eq!(end.x, 0.0);
        assert!((start.y - 0.2).abs() < 1e-2, "start span {start:?}");
        assert!((end.y - 0.9).abs() < 1e-2, "end span {end:?}");
    }
}
