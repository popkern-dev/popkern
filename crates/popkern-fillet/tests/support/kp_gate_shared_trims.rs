use std::collections::HashSet;

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle};

pub struct SharedTrimMetrics {
    pub shared_trim_count: usize,
    pub pcurve_out_of_unit_count: usize,
    pub worst_source_pair: Option<(usize, usize)>,
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
        worst_source_pair: None,
        worst_abs_uv: 0.0,
        worst_g0_gap: 0.0,
        worst_normal_disagreement: 0.0,
    };
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
        }
        let normal_a = oriented_normal(surface_a, uv_a.x, uv_a.y);
        let normal_b = oriented_normal(surface_b, uv_b.x, uv_b.y);
        if let (Some(a), Some(b)) = (normal_a.normalized(), normal_b.normalized()) {
            metrics.worst_normal_disagreement =
                metrics.worst_normal_disagreement.max(1.0 - a.dot(&b).abs());
        }
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
