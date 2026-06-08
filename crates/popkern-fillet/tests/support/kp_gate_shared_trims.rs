use std::collections::{HashMap, HashSet};

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle};

pub struct SharedTrimMetrics {
    pub shared_trim_count: usize,
    pub pcurve_out_of_unit_count: usize,
    pub degenerate_pcurve_count: usize,
    pub degenerate_source_pair: Option<(usize, usize)>,
    pub degenerate_uv_edges: Option<((f64, f64), (f64, f64), (f64, f64), (f64, f64))>,
    pub worst_source_pair: Option<(usize, usize)>,
    pub worst_uv_edges: Option<((f64, f64), (f64, f64), (f64, f64), (f64, f64))>,
    pub worst_pair_best_boundary_gap: f64,
    pub max_boundary_row_shared_trim_count: usize,
    pub worst_abs_uv: f64,
    pub worst_g0_gap: f64,
    pub worst_normal_disagreement: f64,
}

pub fn measure_shared_gregory_trims(
    brep: &BRep,
    faces: &[FaceHandle],
    n_samples: usize,
) -> SharedTrimMetrics {
    measure_shared_gregory_trims_with_sources(brep, faces, &[], n_samples)
}

pub fn measure_shared_gregory_trims_with_sources(
    brep: &BRep,
    faces: &[FaceHandle],
    face_sources: &[(FaceHandle, usize)],
    n_samples: usize,
) -> SharedTrimMetrics {
    let face_set = faces.iter().copied().collect::<HashSet<_>>();
    let mut visited_edges = HashSet::new();
    let mut metrics = SharedTrimMetrics {
        shared_trim_count: 0,
        pcurve_out_of_unit_count: 0,
        degenerate_pcurve_count: 0,
        degenerate_source_pair: None,
        degenerate_uv_edges: None,
        worst_source_pair: None,
        worst_uv_edges: None,
        worst_pair_best_boundary_gap: f64::INFINITY,
        max_boundary_row_shared_trim_count: 0,
        worst_abs_uv: 0.0,
        worst_g0_gap: 0.0,
        worst_normal_disagreement: 0.0,
    };
    let mut boundary_row_counts = HashMap::new();
    for &face in faces {
        let Some(face_data) = brep.faces.get(face) else {
            continue;
        };
        let Some(start) = brep
            .loops
            .get(face_data.outer_loop)
            .map(|loop_| loop_.half_edge)
        else {
            continue;
        };
        let mut cursor = start;
        for _ in 0..brep.half_edges.len().max(1) {
            let he = &brep.half_edges[cursor];
            if visited_edges.insert(he.edge) {
                let twin = &brep.half_edges[he.twin];
                let Some(twin_face) = brep.loops.get(twin.loop_).map(|loop_| loop_.face) else {
                    cursor = he.next;
                    continue;
                };
                if face_set.contains(&twin_face) {
                    measure_shared_trim_pair(
                        brep,
                        face,
                        cursor,
                        twin_face,
                        he.twin,
                        n_samples,
                        face_sources,
                        &mut metrics,
                        &mut boundary_row_counts,
                    );
                }
            }
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
    }
    metrics
}

fn measure_shared_trim_pair(
    brep: &BRep,
    face_a: FaceHandle,
    he_a: popkern_core::topology::HalfEdgeHandle,
    face_b: FaceHandle,
    he_b: popkern_core::topology::HalfEdgeHandle,
    n_samples: usize,
    face_sources: &[(FaceHandle, usize)],
    metrics: &mut SharedTrimMetrics,
    boundary_row_counts: &mut HashMap<(FaceHandle, u8), usize>,
) {
    let Some(surface_a) = brep.faces[face_a].surface.as_ref() else {
        return;
    };
    let Some(surface_b) = brep.faces[face_b].surface.as_ref() else {
        return;
    };
    if !matches!(surface_a.geometry, SurfaceGeometry::Gregory(_))
        || !matches!(surface_b.geometry, SurfaceGeometry::Gregory(_))
    {
        return;
    }
    let Some(pcurve_a) = brep.half_edges[he_a].pcurve.as_ref() else {
        return;
    };
    let Some(pcurve_b) = brep.half_edges[he_b].pcurve.as_ref() else {
        return;
    };
    metrics.shared_trim_count += 1;
    record_boundary_row_count(metrics, boundary_row_counts, face_a, pcurve_a);
    record_boundary_row_count(metrics, boundary_row_counts, face_b, pcurve_b);
    record_pcurve_span(metrics, pcurve_a, pcurve_b, face_sources, face_a, face_b);
    for index in 0..n_samples {
        let fraction = if n_samples <= 1 {
            0.5
        } else {
            index as f64 / (n_samples - 1) as f64
        };
        let uv_a = pcurve_at_fraction(pcurve_a, fraction);
        let uv_b = pcurve_at_fraction(pcurve_b, 1.0 - fraction);
        record_uv(metrics, uv_a);
        record_uv(metrics, uv_b);
        let point_a = surface_a.geometry.point_at(uv_a.x, uv_a.y);
        let point_b = surface_b.geometry.point_at(uv_b.x, uv_b.y);
        let gap = point_a.distance_to(&point_b);
        if gap > metrics.worst_g0_gap {
            metrics.worst_g0_gap = gap;
            metrics.worst_source_pair = source_pair(face_sources, face_a, face_b);
            metrics.worst_uv_edges = Some(edge_uvs(pcurve_a, pcurve_b));
            metrics.worst_pair_best_boundary_gap =
                best_surface_boundary_gap(&surface_a.geometry, &surface_b.geometry);
        }
        let normal_a = oriented_normal(surface_a, uv_a.x, uv_a.y);
        let normal_b = oriented_normal(surface_b, uv_b.x, uv_b.y);
        if let (Some(a), Some(b)) = (normal_a.normalized(), normal_b.normalized()) {
            metrics.worst_normal_disagreement =
                metrics.worst_normal_disagreement.max(1.0 - a.dot(&b).abs());
        }
    }
}

fn record_boundary_row_count(
    metrics: &mut SharedTrimMetrics,
    boundary_row_counts: &mut HashMap<(FaceHandle, u8), usize>,
    face: FaceHandle,
    pcurve: &popkern_core::geometry::pcurve::Pcurve,
) {
    let Some(row) = boundary_row(pcurve) else {
        return;
    };
    let count = boundary_row_counts
        .entry((face, row))
        .and_modify(|count| *count += 1)
        .or_insert(1);
    metrics.max_boundary_row_shared_trim_count =
        metrics.max_boundary_row_shared_trim_count.max(*count);
}

fn boundary_row(pcurve: &popkern_core::geometry::pcurve::Pcurve) -> Option<u8> {
    let (t0, t1) = pcurve.parameter_range();
    let start = pcurve.point_at(t0);
    let end = pcurve.point_at(t1);
    if close(start.x, 0.0) && close(end.x, 0.0) {
        Some(0)
    } else if close(start.x, 1.0) && close(end.x, 1.0) {
        Some(1)
    } else if close(start.y, 0.0) && close(end.y, 0.0) {
        Some(2)
    } else if close(start.y, 1.0) && close(end.y, 1.0) {
        Some(3)
    } else {
        None
    }
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-9
}

fn best_surface_boundary_gap(surface_a: &SurfaceGeometry, surface_b: &SurfaceGeometry) -> f64 {
    let mut best = f64::INFINITY;
    for side_a in 0..4 {
        for side_b in 0..4 {
            best = best.min(surface_boundary_gap(
                surface_a, side_a, surface_b, side_b, false,
            ));
            best = best.min(surface_boundary_gap(
                surface_a, side_a, surface_b, side_b, true,
            ));
        }
    }
    best
}

fn surface_boundary_gap(
    surface_a: &SurfaceGeometry,
    side_a: usize,
    surface_b: &SurfaceGeometry,
    side_b: usize,
    reversed: bool,
) -> f64 {
    let mut worst = 0.0_f64;
    for index in 0..9 {
        let t = index as f64 / 8.0;
        let other_t = if reversed { 1.0 - t } else { t };
        let point_a = surface_a.point_at(side_u(side_a, t), side_v(side_a, t));
        let point_b = surface_b.point_at(side_u(side_b, other_t), side_v(side_b, other_t));
        worst = worst.max(point_a.distance_to(&point_b));
    }
    worst
}

fn side_u(side: usize, t: f64) -> f64 {
    match side {
        0 => 0.0,
        1 => 1.0,
        _ => t,
    }
}

fn side_v(side: usize, t: f64) -> f64 {
    match side {
        2 => 0.0,
        3 => 1.0,
        _ => t,
    }
}

fn source_pair(
    face_sources: &[(FaceHandle, usize)],
    face_a: FaceHandle,
    face_b: FaceHandle,
) -> Option<(usize, usize)> {
    let source_a = face_sources
        .iter()
        .find_map(|&(face, source)| (face == face_a).then_some(source))?;
    let source_b = face_sources
        .iter()
        .find_map(|&(face, source)| (face == face_b).then_some(source))?;
    Some((source_a, source_b))
}

fn record_pcurve_span(
    metrics: &mut SharedTrimMetrics,
    pcurve_a: &popkern_core::geometry::pcurve::Pcurve,
    pcurve_b: &popkern_core::geometry::pcurve::Pcurve,
    face_sources: &[(FaceHandle, usize)],
    face_a: FaceHandle,
    face_b: FaceHandle,
) {
    let span_a = pcurve_span(pcurve_a);
    let span_b = pcurve_span(pcurve_b);
    if span_a <= 1.0e-12 || span_b <= 1.0e-12 {
        metrics.degenerate_pcurve_count += usize::from(span_a <= 1.0e-12);
        metrics.degenerate_pcurve_count += usize::from(span_b <= 1.0e-12);
        if metrics.degenerate_source_pair.is_none() {
            metrics.degenerate_source_pair = source_pair(face_sources, face_a, face_b);
            metrics.degenerate_uv_edges = Some(edge_uvs(pcurve_a, pcurve_b));
        }
    }
}

fn pcurve_span(pcurve: &popkern_core::geometry::pcurve::Pcurve) -> f64 {
    let (t0, t1) = pcurve.parameter_range();
    let start = pcurve.point_at(t0);
    let end = pcurve.point_at(t1);
    ((end.x - start.x).powi(2) + (end.y - start.y).powi(2)).sqrt()
}

fn edge_uvs(
    pcurve_a: &popkern_core::geometry::pcurve::Pcurve,
    pcurve_b: &popkern_core::geometry::pcurve::Pcurve,
) -> ((f64, f64), (f64, f64), (f64, f64), (f64, f64)) {
    let (a0, a1) = pcurve_a.parameter_range();
    let (b0, b1) = pcurve_b.parameter_range();
    let a_start = pcurve_a.point_at(a0);
    let a_end = pcurve_a.point_at(a1);
    let b_start = pcurve_b.point_at(b0);
    let b_end = pcurve_b.point_at(b1);
    (
        (a_start.x, a_start.y),
        (a_end.x, a_end.y),
        (b_start.x, b_start.y),
        (b_end.x, b_end.y),
    )
}

fn record_uv(metrics: &mut SharedTrimMetrics, uv: Vec3) {
    metrics.worst_abs_uv = metrics.worst_abs_uv.max(uv.x.abs()).max(uv.y.abs());
    if uv.x < -1e-9 || uv.x > 1.0 + 1e-9 || uv.y < -1e-9 || uv.y > 1.0 + 1e-9 {
        metrics.pcurve_out_of_unit_count += 1;
    }
}

fn pcurve_at_fraction(pcurve: &popkern_core::geometry::pcurve::Pcurve, fraction: f64) -> Vec3 {
    let t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * fraction;
    let point = pcurve.point_at(t);
    Vec3::new(point.x, point.y, point.z)
}

fn oriented_normal(
    surface: &popkern_core::geometry::surface::SurfaceOnFace,
    u: f64,
    v: f64,
) -> Vec3 {
    let normal = surface.geometry.normal_at(u, v);
    if surface.orientation {
        normal
    } else {
        normal * -1.0
    }
}
