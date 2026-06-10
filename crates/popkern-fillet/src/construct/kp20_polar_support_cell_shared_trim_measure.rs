//! Shared-trim measurement for KP20 support-cell patchwork.

use std::collections::{HashMap, HashSet};

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{GregoryPatch, normal_curvature_along};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle};

use super::{Kp20SupportCellSharedTrimQuality, Kp20SupportCellWorstTrim};

/// Parametric degeneracy guard for normal-curvature samples —
/// rejects only samples where the first fundamental form or the
/// projected trim direction is numerically meaningless. Matches
/// the ring-gate helper value; it is a dimensionless parametric
/// guard, not a model tolerance.
const G2_PARAMETRIC_DEGENERACY_TOL: f64 = 1.0e-14;

/// Pcurve fraction step used to finite-difference the trim
/// direction. Support-cell pcurves are affine, so the difference
/// quotient is direction-exact for them at any step size.
const TRIM_DIRECTION_FRACTION_STEP: f64 = 1.0e-6;

pub(in crate::construct) fn measure_support_cell_shared_trim_quality(
    brep: &BRep,
    faces: &[FaceHandle],
    n_samples: usize,
) -> Kp20SupportCellSharedTrimQuality {
    let face_set = faces.iter().copied().collect::<HashSet<_>>();
    let mut visited_edges = HashSet::new();
    let mut boundary_row_counts = HashMap::new();
    let mut quality = Kp20SupportCellSharedTrimQuality {
        shared_trim_count: 0,
        partial_shared_trim_count: 0,
        pcurve_out_of_unit_count: 0,
        degenerate_pcurve_count: 0,
        max_boundary_row_shared_trim_count: 0,
        worst_abs_uv: 0.0,
        worst_g0_gap: 0.0,
        worst_normal_disagreement: 0.0,
        worst_g2_residual: 0.0,
        g2_unmeasured_sample_count: 0,
        worst_trim: None,
    };
    for &face in faces {
        let Some(start) = brep
            .faces
            .get(face)
            .and_then(|face_data| brep.loops.get(face_data.outer_loop))
            .map(|loop_| loop_.half_edge)
        else {
            continue;
        };
        let mut cursor = start;
        for _ in 0..brep.half_edges.len().max(1) {
            let half_edge = &brep.half_edges[cursor];
            if visited_edges.insert(half_edge.edge) {
                let twin = &brep.half_edges[half_edge.twin];
                if let Some(twin_face) = brep.loops.get(twin.loop_).map(|loop_| loop_.face)
                    && face_set.contains(&twin_face)
                {
                    measure_shared_trim_pair(
                        brep,
                        face,
                        cursor,
                        twin_face,
                        half_edge.twin,
                        n_samples,
                        &mut quality,
                        &mut boundary_row_counts,
                    );
                }
            }
            cursor = half_edge.next;
            if cursor == start {
                break;
            }
        }
    }
    quality
}

fn measure_shared_trim_pair(
    brep: &BRep,
    face_a: FaceHandle,
    half_edge_a: popkern_core::topology::HalfEdgeHandle,
    face_b: FaceHandle,
    half_edge_b: popkern_core::topology::HalfEdgeHandle,
    n_samples: usize,
    quality: &mut Kp20SupportCellSharedTrimQuality,
    boundary_row_counts: &mut HashMap<(FaceHandle, u8), usize>,
) {
    let Some(surface_a) = brep.faces[face_a].surface.as_ref() else {
        return;
    };
    let Some(surface_b) = brep.faces[face_b].surface.as_ref() else {
        return;
    };
    let (SurfaceGeometry::Gregory(patch_a), SurfaceGeometry::Gregory(patch_b)) =
        (&surface_a.geometry, &surface_b.geometry)
    else {
        return;
    };
    let Some(pcurve_a) = brep.half_edges[half_edge_a].pcurve.as_ref() else {
        return;
    };
    let Some(pcurve_b) = brep.half_edges[half_edge_b].pcurve.as_ref() else {
        return;
    };
    let row_a = boundary_row(pcurve_a);
    let row_b = boundary_row(pcurve_b);
    let full_span = full_boundary_span(pcurve_a) && full_boundary_span(pcurve_b);
    quality.shared_trim_count += 1;
    quality.partial_shared_trim_count += usize::from(!full_span);
    record_boundary_row_count(quality, boundary_row_counts, face_a, pcurve_a);
    record_boundary_row_count(quality, boundary_row_counts, face_b, pcurve_b);
    record_pcurve_span(quality, pcurve_a, pcurve_b);
    for index in 0..n_samples {
        let fraction = if n_samples <= 1 {
            0.5
        } else {
            index as f64 / (n_samples - 1) as f64
        };
        let uv_a = pcurve_at_fraction(pcurve_a, fraction);
        let uv_b = pcurve_at_fraction(pcurve_b, 1.0 - fraction);
        record_uv(quality, uv_a);
        record_uv(quality, uv_b);
        let point_a = surface_a.geometry.point_at(uv_a.x, uv_a.y);
        let point_b = surface_b.geometry.point_at(uv_b.x, uv_b.y);
        let g0_gap = point_a.distance_to(&point_b);
        if g0_gap > quality.worst_g0_gap {
            quality.worst_g0_gap = g0_gap;
            quality.worst_trim = Some(Kp20SupportCellWorstTrim {
                face_a,
                row_a,
                face_b,
                row_b,
                sample_fraction: fraction,
                full_span,
            });
        }
        let normal_a = oriented_normal(surface_a, uv_a.x, uv_a.y);
        let normal_b = oriented_normal(surface_b, uv_b.x, uv_b.y);
        if let (Some(a), Some(b)) = (normal_a.normalized(), normal_b.normalized()) {
            quality.worst_normal_disagreement =
                quality.worst_normal_disagreement.max(1.0 - a.dot(&b).abs());
        }
        // Curvature is sampled on the open interior of the trim
        // only: the exact endpoints can land on Gregory corner
        // singularities (e.g. the polar extraordinary point),
        // where second derivatives are undefined by
        // construction. Same convention as the ring gates'
        // t ∈ [0.05, 0.95] sampling; G0/G1 still cover the
        // endpoints above.
        if fraction > 0.0 && fraction < 1.0 {
            record_trim_g2_residual(
                quality,
                patch_a,
                surface_a.orientation,
                pcurve_a,
                uv_a,
                patch_b,
                surface_b.orientation,
                uv_b,
                fraction,
            );
        }
    }
}

/// Normal-curvature (G2) residual at one trim sample. Both faces
/// are measured along the same ambient trim direction; each
/// curvature is expressed against its face's *oriented* normal so
/// the two values are comparable on a consistently-oriented
/// manifold. Samples that cannot be measured (unsupported patch
/// variant, parametric degeneracy) are counted — callers fail
/// closed on a non-zero unmeasured count.
#[expect(
    clippy::too_many_arguments,
    reason = "single trim-sample record taking both faces' sample state"
)]
fn record_trim_g2_residual(
    quality: &mut Kp20SupportCellSharedTrimQuality,
    patch_a: &GregoryPatch,
    orientation_a: bool,
    pcurve_a: &popkern_core::geometry::pcurve::Pcurve,
    uv_a: Vec3,
    patch_b: &GregoryPatch,
    orientation_b: bool,
    uv_b: Vec3,
    fraction: f64,
) {
    let Some(direction) = trim_direction(patch_a, pcurve_a, uv_a, fraction) else {
        quality.g2_unmeasured_sample_count += 1;
        return;
    };
    let curvature_a = normal_curvature_along(
        patch_a,
        uv_a.x,
        uv_a.y,
        direction,
        G2_PARAMETRIC_DEGENERACY_TOL,
    );
    let curvature_b = normal_curvature_along(
        patch_b,
        uv_b.x,
        uv_b.y,
        direction,
        G2_PARAMETRIC_DEGENERACY_TOL,
    );
    let (Ok(curvature_a), Ok(curvature_b)) = (curvature_a, curvature_b) else {
        quality.g2_unmeasured_sample_count += 1;
        return;
    };
    let oriented_a = oriented_curvature(curvature_a, orientation_a);
    let oriented_b = oriented_curvature(curvature_b, orientation_b);
    let scale = oriented_a.abs().max(oriented_b.abs());
    let residual = (oriented_a - oriented_b).abs() / (1.0 + scale);
    quality.worst_g2_residual = quality.worst_g2_residual.max(residual);
}

fn oriented_curvature(curvature: f64, orientation: bool) -> f64 {
    if orientation { curvature } else { -curvature }
}

/// Ambient trim direction at `fraction`, via the chain rule
/// `dS/dt = S_u · du/dt + S_v · dv/dt` with the pcurve
/// difference quotient (direction-exact for affine pcurves).
/// Uses a one-sided difference at the trim endpoints.
fn trim_direction(
    patch: &GregoryPatch,
    pcurve: &popkern_core::geometry::pcurve::Pcurve,
    uv: Vec3,
    fraction: f64,
) -> Option<Vec3> {
    let low = (fraction - TRIM_DIRECTION_FRACTION_STEP).max(0.0);
    let high = (fraction + TRIM_DIRECTION_FRACTION_STEP).min(1.0);
    let uv_low = pcurve_at_fraction(pcurve, low);
    let uv_high = pcurve_at_fraction(pcurve, high);
    let du = uv_high.x - uv_low.x;
    let dv = uv_high.y - uv_low.y;
    let tangent_u = patch.partial_u(uv.x, uv.y);
    let tangent_v = patch.partial_v(uv.x, uv.y);
    (tangent_u * du + tangent_v * dv).normalized()
}

fn record_boundary_row_count(
    quality: &mut Kp20SupportCellSharedTrimQuality,
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
    quality.max_boundary_row_shared_trim_count =
        quality.max_boundary_row_shared_trim_count.max(*count);
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

fn record_pcurve_span(
    quality: &mut Kp20SupportCellSharedTrimQuality,
    pcurve_a: &popkern_core::geometry::pcurve::Pcurve,
    pcurve_b: &popkern_core::geometry::pcurve::Pcurve,
) {
    let span_a = pcurve_span(pcurve_a);
    let span_b = pcurve_span(pcurve_b);
    quality.degenerate_pcurve_count += usize::from(span_a <= 1.0e-12);
    quality.degenerate_pcurve_count += usize::from(span_b <= 1.0e-12);
}

fn pcurve_span(pcurve: &popkern_core::geometry::pcurve::Pcurve) -> f64 {
    let (t0, t1) = pcurve.parameter_range();
    let start = pcurve.point_at(t0);
    let end = pcurve.point_at(t1);
    ((end.x - start.x).powi(2) + (end.y - start.y).powi(2)).sqrt()
}

fn full_boundary_span(pcurve: &popkern_core::geometry::pcurve::Pcurve) -> bool {
    let (t0, t1) = pcurve.parameter_range();
    let start = pcurve.point_at(t0);
    let end = pcurve.point_at(t1);
    (close(start.x, end.x) && close((start.y - end.y).abs(), 1.0))
        || (close(start.y, end.y) && close((start.x - end.x).abs(), 1.0))
}

fn record_uv(quality: &mut Kp20SupportCellSharedTrimQuality, uv: Vec3) {
    quality.worst_abs_uv = quality.worst_abs_uv.max(uv.x.abs()).max(uv.y.abs());
    if uv.x < -1.0e-9 || uv.x > 1.0 + 1.0e-9 || uv.y < -1.0e-9 || uv.y > 1.0 + 1.0e-9 {
        quality.pcurve_out_of_unit_count += 1;
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

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-9
}
