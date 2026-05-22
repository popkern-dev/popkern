//! Measured pcurve repair for late export-readiness gates.

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::edge_triple::{EdgeAnchors, EdgeTripleOptions, compute_edge_geometry};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, pcurve::Pcurve};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::pcurve_refit::{fit_pcurve_for_half_edge, polyline_uv_curve};

pub(crate) fn refit_measured_pcurve_blockers(brep: &mut BRep, tolerance: f64) -> usize {
    let half_edges = brep.half_edges.keys().collect::<Vec<_>>();
    let mut repaired = 0;
    for half_edge in half_edges {
        let Some(current) = brep.measure_half_edge_pcurve_deviation(half_edge) else {
            continue;
        };
        if current.max_deviation <= tolerance {
            continue;
        }
        let Some(edge) = brep
            .half_edges
            .get(half_edge)
            .map(|half_edge| half_edge.edge)
        else {
            continue;
        };
        let Some(curve_on_edge) = brep.edges.get(edge).and_then(|edge| edge.curve.clone()) else {
            continue;
        };
        let Some(surface) = surface_for_half_edge(brep, half_edge) else {
            continue;
        };
        let origin = brep.vertices[brep.half_edges[half_edge].origin].position;
        let curve_start = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
        let curve_end = curve_on_edge.geometry.point_at(curve_on_edge.t_end);
        let forward = origin.distance_to(&curve_start) <= origin.distance_to(&curve_end);
        let Some(new_pcurve) = reparameterize_existing_pcurve_trace(
            brep,
            half_edge,
            &curve_on_edge,
            &surface,
            forward,
            257,
        )
        .or_else(|| {
            fit_pcurve_for_half_edge(
                &curve_on_edge.geometry,
                curve_on_edge.t_start,
                curve_on_edge.t_end,
                &surface,
                forward,
                257,
            )
        }) else {
            continue;
        };
        let previous = brep.half_edges[half_edge].pcurve.clone();
        brep.half_edges[half_edge].pcurve = Some(new_pcurve.into());
        let accepted = brep
            .measure_half_edge_pcurve_deviation(half_edge)
            .is_some_and(|report| {
                report.max_deviation <= tolerance && report.max_deviation < current.max_deviation
            });
        if accepted {
            repaired += 1;
        } else {
            brep.half_edges[half_edge].pcurve = previous;
        }
    }
    repaired
}

pub(crate) fn repair_measured_plane_cylinder_pcurves(brep: &mut BRep, tolerance: f64) -> usize {
    let edges = brep.edges.keys().collect::<Vec<_>>();
    let mut repaired = 0;
    for edge in edges {
        let Some(report) = brep.measure_edge_pcurve_deviation(edge) else {
            continue;
        };
        if report.max_deviation <= tolerance {
            continue;
        }
        if repair_plane_cylinder_edge_from_triple(brep, edge, tolerance) {
            repaired += 1;
        }
    }
    repaired
}

fn reparameterize_existing_pcurve_trace(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve_on_edge: &CurveOnEdge,
    surface: &SurfaceGeometry,
    forward: bool,
    samples: usize,
) -> Option<NurbsCurve> {
    let pcurve = brep.half_edges.get(half_edge)?.pcurve.as_ref()?;
    let (pcurve_start, pcurve_end) = pcurve.parameter_range();
    if !pcurve_start.is_finite()
        || !pcurve_end.is_finite()
        || (pcurve_end - pcurve_start).abs() <= f64::EPSILON
    {
        return None;
    }
    let samples = samples.max(2);
    let mut points = Vec::with_capacity(samples);
    let mut params = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let curve_fraction = if forward { fraction } else { 1.0 - fraction };
        let curve_t =
            curve_on_edge.t_start + (curve_on_edge.t_end - curve_on_edge.t_start) * curve_fraction;
        let target = curve_on_edge.geometry.point_at(curve_t);
        let pcurve_t = closest_pcurve_trace_parameter(surface, pcurve, &target, 512)?;
        let uv = pcurve.point_at(pcurve_t);
        points.push(Point3::new(uv.x, uv.y, 0.0));
        params.push(fraction);
    }
    polyline_uv_curve(&points, &params)
}

fn closest_pcurve_trace_parameter(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    target: &Point3,
    samples: usize,
) -> Option<f64> {
    let (start, end) = pcurve.parameter_range();
    let samples = samples.max(2);
    let mut best_index = 0usize;
    let mut best_distance = f64::INFINITY;
    for index in 0..=samples {
        let fraction = index as f64 / samples as f64;
        let parameter = start + (end - start) * fraction;
        let uv = pcurve.point_at(parameter);
        let distance = surface.point_at(uv.x, uv.y).distance_to(target);
        if distance < best_distance {
            best_distance = distance;
            best_index = index;
        }
    }
    let lo_index = best_index.saturating_sub(1);
    let hi_index = (best_index + 1).min(samples);
    let mut lo = start + (end - start) * (lo_index as f64 / samples as f64);
    let mut hi = start + (end - start) * (hi_index as f64 / samples as f64);
    for _ in 0..24 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        let d1 = pcurve_trace_distance_at(surface, pcurve, m1, target);
        let d2 = pcurve_trace_distance_at(surface, pcurve, m2, target);
        if d1 <= d2 {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    Some((lo + hi) * 0.5)
}

fn pcurve_trace_distance_at(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    parameter: f64,
    target: &Point3,
) -> f64 {
    let uv = pcurve.point_at(parameter);
    surface.point_at(uv.x, uv.y).distance_to(target)
}

fn repair_plane_cylinder_edge_from_triple(
    brep: &mut BRep,
    edge: EdgeHandle,
    tolerance: f64,
) -> bool {
    let half_edges = edge_half_edges(brep, edge);
    if half_edges.len() != 2 {
        return false;
    }
    let Some((_face_a, surface_a)) = face_surface_for_half_edge(brep, half_edges[0]) else {
        return false;
    };
    let Some((_face_b, surface_b)) = face_surface_for_half_edge(brep, half_edges[1]) else {
        return false;
    };
    if !is_plane_cylinder_pair(&surface_a, &surface_b) {
        return false;
    }
    let start = brep.vertices[brep.half_edges[half_edges[0]].origin].position;
    let end = brep.vertices[brep.half_edges[brep.half_edges[half_edges[0]].next].origin].position;
    let Ok(triple) = compute_edge_geometry(
        &surface_a,
        &surface_b,
        EdgeAnchors::Open { start, end },
        None,
        &popkern_core::tolerance::Tolerance::default(),
        &EdgeTripleOptions::default(),
    ) else {
        return false;
    };
    let previous_pcurves = half_edges
        .iter()
        .map(|half_edge| (*half_edge, brep.half_edges[*half_edge].pcurve.clone()))
        .collect::<Vec<_>>();
    brep.half_edges[half_edges[0]].pcurve = Some(triple.pcurve_a.clone());
    brep.half_edges[half_edges[1]].pcurve = Some(triple.pcurve_b.clone());
    let accepted = brep
        .measure_edge_pcurve_deviation(edge)
        .is_some_and(|updated| updated.max_deviation <= tolerance)
        && brep.validate().is_ok();
    if accepted {
        true
    } else {
        for (half_edge, pcurve) in previous_pcurves {
            brep.half_edges[half_edge].pcurve = pcurve;
        }
        false
    }
}

fn is_plane_cylinder_pair(surface_a: &SurfaceGeometry, surface_b: &SurfaceGeometry) -> bool {
    matches!(surface_a, SurfaceGeometry::Plane { .. })
        && matches!(surface_b, SurfaceGeometry::Cylinder { .. })
        || matches!(surface_a, SurfaceGeometry::Cylinder { .. })
            && matches!(surface_b, SurfaceGeometry::Plane { .. })
}

fn face_surface_for_half_edge(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<(FaceHandle, SurfaceGeometry)> {
    let face = brep.loops.get(brep.half_edges.get(half_edge)?.loop_)?.face;
    let surface = brep.faces.get(face)?.surface.as_ref()?.geometry.clone();
    Some((face, surface))
}

fn surface_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<SurfaceGeometry> {
    face_surface_for_half_edge(brep, half_edge).map(|(_, surface)| surface)
}

fn edge_half_edges(brep: &BRep, edge: EdgeHandle) -> Vec<HalfEdgeHandle> {
    let Some(primary) = brep.edges.get(edge).map(|edge| edge.half_edge) else {
        return Vec::new();
    };
    let mut out = vec![primary];
    if let Some(twin) = brep.half_edges.get(primary).map(|half_edge| half_edge.twin)
        && twin != primary
    {
        out.push(twin);
    }
    out
}
