//! Sample-quality gates for cap/support pcurve completion.

use popkern_core::geometry::Point3;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, HalfEdgeHandle};

pub(super) fn pcurve_surface_sample_gap(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    samples: &[Point3],
) -> f64 {
    if samples.len() < 2 {
        return f64::INFINITY;
    }
    let (t0, t1) = pcurve.parameter_range();
    samples
        .iter()
        .enumerate()
        .map(|(index, sample)| {
            let s = index as f64 / (samples.len() - 1) as f64;
            let point = pcurve_surface_point(surface, pcurve, t0 + s * (t1 - t0));
            if point_is_finite(&point) {
                point.distance_to(sample)
            } else {
                f64::INFINITY
            }
        })
        .fold(0.0, f64::max)
}

pub(super) fn oriented_pcurve_sample_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    samples: &[Point3],
) -> f64 {
    if samples.len() < 2 {
        return f64::INFINITY;
    }
    let cap_start = brep.vertices[brep.half_edges[half_edge].origin].position;
    let (t0, t1) = pcurve.parameter_range();
    let p0 = pcurve_surface_point(surface, pcurve, t0);
    let p1 = pcurve_surface_point(surface, pcurve, t1);
    if !point_is_finite(&p0) || !point_is_finite(&p1) {
        return f64::INFINITY;
    }
    let forward = p0.distance_to(&cap_start) <= p1.distance_to(&cap_start);
    samples
        .iter()
        .enumerate()
        .map(|(index, sample)| {
            let s = index as f64 / (samples.len() - 1) as f64;
            let local = if forward { s } else { 1.0 - s };
            let point = pcurve_surface_point(surface, pcurve, t0 + local * (t1 - t0));
            if point_is_finite(&point) {
                point.distance_to(sample)
            } else {
                f64::INFINITY
            }
        })
        .fold(0.0, f64::max)
}

pub(super) fn support_pcurve_candidate_improves_samples(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    support_surface: &SurfaceGeometry,
    candidate: &Pcurve,
    samples: &[Point3],
) -> bool {
    let candidate_gap = pcurve_surface_sample_gap(support_surface, candidate, samples);
    if !candidate_gap.is_finite() {
        return false;
    }
    let Some(existing) = brep
        .half_edges
        .get(half_edge)
        .and_then(|he| he.pcurve.as_ref())
    else {
        return true;
    };
    let existing_gap =
        oriented_pcurve_sample_gap(brep, half_edge, support_surface, existing, samples);
    candidate_gap + 1.0e-12 < existing_gap
}

fn pcurve_surface_point(surface: &SurfaceGeometry, pcurve: &Pcurve, t: f64) -> Point3 {
    let uv = pcurve.point_at(t);
    surface.point_at(uv.x, uv.y)
}

fn point_is_finite(point: &Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}
