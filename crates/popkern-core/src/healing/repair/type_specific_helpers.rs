//! Private helpers for geometry-type-specific healing.

use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, FaceHandle};
use popkern_math::PortableMath;

use super::{SurfaceType, TypeDetectionResult};

/// Get faces adjacent to a given face (sharing an edge).
pub(super) fn adjacent_faces(brep: &BRep, face: FaceHandle) -> Vec<FaceHandle> {
    let Some(face_data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return Vec::new();
    };

    let start_he = loop_data.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    let mut neighbors = Vec::new();
    let mut current = start_he;
    let max_walk = 1000;
    for _ in 0..max_walk {
        if let Some(he) = brep.half_edges.get(current) {
            if let Some(twin) = brep.half_edges.get(he.twin)
                && let Some(twin_loop) = brep.loops.get(twin.loop_)
            {
                let neighbor_face = twin_loop.face;
                if neighbor_face != face && !neighbors.contains(&neighbor_face) {
                    neighbors.push(neighbor_face);
                }
            }
            current = he.next;
            if current == start_he {
                break;
            }
        } else {
            break;
        }
    }

    neighbors
}

/// Get edge handles from the outer loop of a face.
pub(super) fn edges_shared_with_neighbors(brep: &BRep, face: FaceHandle) -> Vec<EdgeHandle> {
    let Some(face_data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return Vec::new();
    };

    let start_he = loop_data.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    let mut edges = Vec::new();
    let mut current = start_he;
    let max_walk = 1000;
    for _ in 0..max_walk {
        if let Some(he) = brep.half_edges.get(current) {
            if !edges.contains(&he.edge) {
                edges.push(he.edge);
            }
            current = he.next;
            if current == start_he {
                break;
            }
        } else {
            break;
        }
    }

    edges
}

/// Try to detect a fillet surface (G1-continuous with 2+ adjacent faces).
pub(super) fn try_detect_fillet(
    brep: &BRep,
    face: FaceHandle,
    neighbors: &[FaceHandle],
    tolerance: &Tolerance,
) -> Option<TypeDetectionResult> {
    if neighbors.len() < 2 {
        return None;
    }

    let face_data = brep.faces.get(face)?;
    let sof = face_data.surface.as_ref()?;

    let mut g1_count = 0usize;
    let mut total_angle = 0.0_f64;
    let mut checked = 0usize;

    for neighbor in neighbors {
        let neighbor_data = brep.faces.get(*neighbor)?;
        let neighbor_sof = neighbor_data.surface.as_ref()?;

        let angle = sample_boundary_normal_angle(&sof.geometry, &neighbor_sof.geometry, 8);
        if angle < tolerance.angular() * 10.0 {
            g1_count += 1;
            total_angle += angle;
        }
        checked += 1;
    }

    if g1_count >= 2 {
        let avg_angle = if g1_count > 0 {
            total_angle / g1_count as f64
        } else {
            0.0
        };
        let confidence = (1.0 - avg_angle / (tolerance.angular() * 10.0)).clamp(0.0, 1.0);
        Some(TypeDetectionResult {
            face,
            surface_type: SurfaceType::Fillet,
            confidence,
            description: format!(
                "G1-continuous with {g1_count}/{checked} neighbors (avg angle {avg_angle:.6e})"
            ),
        })
    } else {
        None
    }
}

/// Try to detect an offset surface (constant distance from a neighbor).
pub(super) fn try_detect_offset(
    brep: &BRep,
    face: FaceHandle,
    surface: &SurfaceGeometry,
    neighbors: &[FaceHandle],
    tolerance: &Tolerance,
) -> Option<TypeDetectionResult> {
    for neighbor in neighbors {
        let n_data = brep.faces.get(*neighbor)?;
        let n_sof = n_data.surface.as_ref()?;
        let (median, max_dev, _mean_dev) = measure_offset_consistency(surface, &n_sof.geometry, 8);

        if median > tolerance.linear() * 5.0 && max_dev < tolerance.linear() * 10.0 {
            let confidence = (1.0 - max_dev / (tolerance.linear() * 10.0)).clamp(0.0, 1.0);
            return Some(TypeDetectionResult {
                face,
                surface_type: SurfaceType::Offset,
                confidence,
                description: format!(
                    "constant offset distance {median:.6e} (max dev {max_dev:.6e})"
                ),
            });
        }
    }
    None
}

/// Try to detect a swept surface from cross-section uniformity.
pub(super) fn try_detect_swept(
    face: FaceHandle,
    nurbs: &NurbsSurface,
    _tolerance: &Tolerance,
) -> Option<TypeDetectionResult> {
    if nurbs.num_u < 3 || nurbs.num_v < 3 {
        return None;
    }

    let (var_ratio, max_shape_dev) = measure_cross_section_uniformity(nurbs, 5);

    if var_ratio < 0.05 && max_shape_dev < 0.1 {
        let confidence = (1.0 - var_ratio / 0.05).clamp(0.0, 1.0);
        Some(TypeDetectionResult {
            face,
            surface_type: SurfaceType::Swept,
            confidence,
            description: format!(
                "uniform cross-sections (length var {var_ratio:.4}, shape dev {max_shape_dev:.4})"
            ),
        })
    } else {
        None
    }
}

/// Measure distance consistency between two surfaces.
/// Returns (median_distance, max_deviation_from_median, mean_deviation).
pub(super) fn measure_offset_consistency(
    surface: &dyn SurfaceEval,
    base: &dyn SurfaceEval,
    num_samples: usize,
) -> (f64, f64, f64) {
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let u_lo = u0.max(-100.0);
    let u_hi = u1.min(100.0);
    let v_lo = v0.max(-100.0);
    let v_hi = v1.min(100.0);

    let n = num_samples.max(2);
    let mut distances = Vec::with_capacity(n * n);

    for i in 0..n {
        let u = u_lo + (u_hi - u_lo) * (i as f64 + 0.5) / n as f64;
        for j in 0..n {
            let v = v_lo + (v_hi - v_lo) * (j as f64 + 0.5) / n as f64;
            let pt = surface.point_at(u, v);
            if let Some((bu, bv)) = base.closest_point(&pt) {
                let base_pt = base.point_at(bu, bv);
                distances.push(pt.distance_to(&base_pt));
            }
        }
    }

    if distances.is_empty() {
        return (0.0, 0.0, 0.0);
    }

    distances.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = if distances.len() % 2 == 1 {
        distances[distances.len() / 2]
    } else {
        (distances[distances.len() / 2 - 1] + distances[distances.len() / 2]) * 0.5
    };

    let max_dev = distances
        .iter()
        .map(|d| (d - median).abs())
        .fold(0.0_f64, f64::max);
    let mean_dev =
        distances.iter().map(|d| (d - median).abs()).sum::<f64>() / distances.len() as f64;

    (median, max_dev, mean_dev)
}

/// Measure cross-section uniformity of a NURBS surface.
///
/// Returns (length_variance_ratio, max_shape_deviation).
pub(super) fn measure_cross_section_uniformity(
    surface: &NurbsSurface,
    num_sections: usize,
) -> (f64, f64) {
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let n = num_sections.max(2);
    let sample_v = 20;

    let mut lengths = Vec::with_capacity(n);
    let mut sections: Vec<Vec<Point3>> = Vec::with_capacity(n);

    for i in 0..n {
        let u = u0 + (u1 - u0) * (i as f64 + 0.5) / n as f64;
        let mut pts = Vec::with_capacity(sample_v);
        for j in 0..sample_v {
            let v = v0 + (v1 - v0) * j as f64 / (sample_v - 1).max(1) as f64;
            pts.push(surface.point_at(u, v));
        }
        let length: f64 = pts.windows(2).map(|w| w[0].distance_to(&w[1])).sum();
        lengths.push(length);
        sections.push(pts);
    }

    let mean_length = lengths.iter().sum::<f64>() / lengths.len() as f64;
    let var_ratio = if mean_length > 1e-15 {
        let variance = lengths
            .iter()
            .map(|l| (l - mean_length) * (l - mean_length))
            .sum::<f64>()
            / lengths.len() as f64;
        variance.sqrt() / mean_length
    } else {
        0.0
    };

    let mut max_shape_dev = 0.0_f64;
    if let Some(ref_section) = sections.first() {
        for section in sections.iter().skip(1) {
            for (rp, sp) in ref_section.iter().zip(section.iter()) {
                let dist = rp.distance_to(sp);
                let normalized = if mean_length > 1e-15 {
                    dist / mean_length
                } else {
                    dist
                };
                max_shape_dev = max_shape_dev.max(normalized);
            }
        }
    }

    (var_ratio, max_shape_dev)
}

/// Sample boundary normal angles between two surfaces.
///
/// Returns the best (minimum) normal angle found across samples.
fn sample_boundary_normal_angle(
    surface_a: &dyn SurfaceEval,
    surface_b: &dyn SurfaceEval,
    num_samples: usize,
) -> f64 {
    let ((ua0, ua1), (va0, va1)) = surface_a.parameter_range();
    let ua_lo = ua0.max(-100.0);
    let ua_hi = ua1.min(100.0);
    let va_lo = va0.max(-100.0);
    let va_hi = va1.min(100.0);

    let n = num_samples.max(2);
    let mut best_angle = std::f64::consts::PI;

    for i in 0..n {
        let ua = ua_lo + (ua_hi - ua_lo) * (i as f64 + 0.5) / n as f64;
        let va = (va_lo + va_hi) * 0.5;
        let pt_a = surface_a.point_at(ua, va);

        if let Some((ub, vb)) = surface_b.closest_point(&pt_a) {
            let na = surface_a.normal_at(ua, va);
            let nb = surface_b.normal_at(ub, vb);
            let la = na.length();
            let lb = nb.length();
            if la > 1e-15 && lb > 1e-15 {
                let cos = (na.dot(&nb) / (la * lb)).abs().clamp(0.0, 1.0);
                let angle = cos.pacos();
                best_angle = best_angle.min(angle);
            }
        }
    }

    best_angle
}

/// Find the best offset base face among neighbors.
pub(super) fn find_best_offset_base(
    brep: &BRep,
    surface: &SurfaceGeometry,
    neighbors: &[FaceHandle],
) -> Option<FaceHandle> {
    let mut best: Option<(FaceHandle, f64)> = None;

    for neighbor in neighbors {
        let n_data = brep.faces.get(*neighbor)?;
        let n_sof = n_data.surface.as_ref()?;
        let (median, max_dev, _) = measure_offset_consistency(surface, &n_sof.geometry, 8);
        if median > 1e-12 {
            let score = max_dev / median;
            match &best {
                Some((_, bs)) if score >= *bs => {}
                _ => best = Some((*neighbor, score)),
            }
        }
    }

    best.map(|(fh, _)| fh)
}

/// Adjust NURBS control points to restore constant offset distance.
pub(super) fn adjust_offset_control_points(
    nurbs: &mut NurbsSurface,
    base: &dyn SurfaceEval,
    target_distance: f64,
    max_displacement: f64,
) -> f64 {
    let damping = 0.5_f64;
    let mut worst_disp = 0.0_f64;

    for idx in 0..nurbs.control_points.len() {
        let cp = nurbs.control_points[idx];
        let Some((bu, bv)) = base.closest_point(&cp) else {
            continue;
        };
        let base_pt = base.point_at(bu, bv);
        let current_dist = cp.distance_to(&base_pt);

        if current_dist < 1e-15 {
            continue;
        }

        let dist_error = current_dist - target_distance;
        if dist_error.abs() < 1e-12 {
            continue;
        }

        let dir = Vec3::new(
            (cp.x - base_pt.x) / current_dist,
            (cp.y - base_pt.y) / current_dist,
            (cp.z - base_pt.z) / current_dist,
        );

        let shift = -dist_error * damping;
        let shift_clamped = shift.clamp(-max_displacement, max_displacement);
        let new_cp = Point3::new(
            cp.x + dir.x * shift_clamped,
            cp.y + dir.y * shift_clamped,
            cp.z + dir.z * shift_clamped,
        );
        let disp = cp.distance_to(&new_cp);
        worst_disp = worst_disp.max(disp);
        nurbs.control_points[idx] = new_cp;
    }

    worst_disp
}

/// Adjust NURBS control points to restore cross-section uniformity.
pub(super) fn adjust_swept_control_points(nurbs: &mut NurbsSurface, max_displacement: f64) -> f64 {
    if nurbs.num_u < 3 {
        return 0.0;
    }

    let damping = 0.3_f64;
    let nv = nurbs.num_v;
    let nu = nurbs.num_u;

    // Reference profile: average of first and last u-rows.
    let mut ref_profile = Vec::with_capacity(nv);
    for j in 0..nv {
        let first = nurbs.control_points[j];
        let last = nurbs.control_points[(nu - 1) * nv + j];
        ref_profile.push(Point3::new(
            (first.x + last.x) * 0.5,
            (first.y + last.y) * 0.5,
            (first.z + last.z) * 0.5,
        ));
    }

    let mut worst_disp = 0.0_f64;

    for i in 1..nu - 1 {
        let mut row_cx = 0.0_f64;
        let mut row_cy = 0.0_f64;
        let mut row_cz = 0.0_f64;
        let mut ref_cx = 0.0_f64;
        let mut ref_cy = 0.0_f64;
        let mut ref_cz = 0.0_f64;

        for (j, rp) in ref_profile.iter().enumerate().take(nv) {
            let cp = nurbs.control_points[i * nv + j];
            row_cx += cp.x;
            row_cy += cp.y;
            row_cz += cp.z;
            ref_cx += rp.x;
            ref_cy += rp.y;
            ref_cz += rp.z;
        }

        let inv_nv = 1.0 / nv as f64;
        row_cx *= inv_nv;
        row_cy *= inv_nv;
        row_cz *= inv_nv;
        ref_cx *= inv_nv;
        ref_cy *= inv_nv;
        ref_cz *= inv_nv;

        let offset = Vec3::new(row_cx - ref_cx, row_cy - ref_cy, row_cz - ref_cz);

        for (j, rp) in ref_profile.iter().enumerate().take(nv) {
            let idx = i * nv + j;
            let cp = nurbs.control_points[idx];
            let target = Point3::new(rp.x + offset.x, rp.y + offset.y, rp.z + offset.z);
            let new_cp = Point3::new(
                cp.x + (target.x - cp.x) * damping,
                cp.y + (target.y - cp.y) * damping,
                cp.z + (target.z - cp.z) * damping,
            );
            let disp = cp.distance_to(&new_cp);
            if disp <= max_displacement {
                worst_disp = worst_disp.max(disp);
                nurbs.control_points[idx] = new_cp;
            }
        }
    }

    worst_disp
}
