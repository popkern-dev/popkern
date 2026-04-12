//! NURBS control-point adjustment helpers for G1/G2 continuity enforcement.

use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::topology::{BRep, FaceHandle};

/// Which boundary of the parameter domain an edge lies on.
#[derive(Debug, Clone, Copy)]
pub(super) enum Boundary {
    UMin,
    UMax,
    VMin,
    VMax,
}

/// Identify which parameter boundary a UV sample is closest to.
pub(super) fn identify_boundary(
    (u, v): (f64, f64),
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
) -> Boundary {
    let du_min = (u - u_min).abs();
    let du_max = (u - u_max).abs();
    let dv_min = (v - v_min).abs();
    let dv_max = (v - v_max).abs();
    let min_d = du_min.min(du_max).min(dv_min).min(dv_max);
    if (min_d - du_min).abs() < 1e-15 {
        Boundary::UMin
    } else if (min_d - du_max).abs() < 1e-15 {
        Boundary::UMax
    } else if (min_d - dv_min).abs() < 1e-15 {
        Boundary::VMin
    } else {
        Boundary::VMax
    }
}

/// Adjust NURBS control points on adjacent faces for G1 continuity.
///
/// For each NURBS face, identifies the boundary row of control points
/// nearest the shared edge, then shifts the next-to-boundary row so
/// the cross-derivative aligns with the target normal direction.
pub(super) fn adjust_for_g1(
    brep: &mut BRep,
    face1: FaceHandle,
    face2: FaceHandle,
    samples: &[(f64, f64, f64, f64)],
    max_displacement: f64,
) -> f64 {
    let damping = 0.5_f64;
    let surf_a = brep.faces.get(face1).and_then(|f| f.surface.clone());
    let surf_b = brep.faces.get(face2).and_then(|f| f.surface.clone());
    let (Some(sa), Some(sb)) = (surf_a, surf_b) else {
        return 0.0;
    };

    let targets: Vec<_> = samples
        .iter()
        .map(|&(ua, va, ub, vb)| {
            let na = sa.geometry.normal_at(ua, va);
            let nb = sb.geometry.normal_at(ub, vb);
            Vec3::new(
                (na.x + nb.x) * 0.5,
                (na.y + nb.y) * 0.5,
                (na.z + nb.z) * 0.5,
            )
            .normalized()
        })
        .collect();

    let mut worst = 0.0_f64;
    for &fh in &[face1, face2] {
        let nurbs = match extract_nurbs(brep, fh) {
            Some(n) => n,
            None => continue,
        };
        let sample_uv = if fh == face1 {
            (samples[0].0, samples[0].1)
        } else {
            (samples[0].2, samples[0].3)
        };
        let ((u_min, u_max), (v_min, v_max)) = nurbs.parameter_range();
        let boundary = identify_boundary(sample_uv, u_min, u_max, v_min, v_max);
        let mut nc = nurbs;
        let d = shift_boundary_row(&mut nc, boundary, &targets, damping, max_displacement);
        worst = worst.max(d);
        write_nurbs(brep, fh, nc);
    }
    worst
}

/// Adjust NURBS control points on adjacent faces for G2 continuity.
///
/// Targets the second row of control points to smooth curvature
/// across the boundary.
pub(super) fn adjust_for_g2(
    brep: &mut BRep,
    face1: FaceHandle,
    face2: FaceHandle,
    samples: &[(f64, f64, f64, f64)],
    max_displacement: f64,
) -> f64 {
    let damping = 0.3_f64;
    let mut worst = 0.0_f64;

    for &fh in &[face1, face2] {
        let nurbs = match extract_nurbs_g2(brep, fh) {
            Some(n) => n,
            None => continue,
        };
        let sample_uv = if fh == face1 {
            (samples[0].0, samples[0].1)
        } else {
            (samples[0].2, samples[0].3)
        };
        let ((u_min, u_max), (v_min, v_max)) = nurbs.parameter_range();
        let boundary = identify_boundary(sample_uv, u_min, u_max, v_min, v_max);
        let mut nc = nurbs;
        let d = shift_second_row(&mut nc, boundary, damping, max_displacement);
        worst = worst.max(d);
        write_nurbs(brep, fh, nc);
    }
    worst
}

/// Extract a cloned NURBS surface from a face (any degree).
fn extract_nurbs(brep: &BRep, fh: FaceHandle) -> Option<NurbsSurface> {
    let face = brep.faces.get(fh)?;
    let surf = face.surface.as_ref()?;
    match &surf.geometry {
        SurfaceGeometry::Nurbs(n) => Some(n.clone()),
        _ => None,
    }
}

/// Extract a cloned NURBS surface from a face, requiring degree >= 3.
fn extract_nurbs_g2(brep: &BRep, fh: FaceHandle) -> Option<NurbsSurface> {
    let n = extract_nurbs(brep, fh)?;
    if n.degree_u >= 3 && n.degree_v >= 3 {
        Some(n)
    } else {
        None
    }
}

/// Write a modified NURBS surface back to a face.
fn write_nurbs(brep: &mut BRep, fh: FaceHandle, nurbs: NurbsSurface) {
    if let Some(face) = brep.faces.get_mut(fh)
        && let Some(ref mut s) = face.surface
    {
        s.geometry = SurfaceGeometry::Nurbs(nurbs);
    }
}

/// Shift the next-to-boundary row toward target normals.
fn shift_boundary_row(
    nurbs: &mut NurbsSurface,
    boundary: Boundary,
    targets: &[Option<Vec3>],
    damping: f64,
    max_disp: f64,
) -> f64 {
    let mut worst = 0.0_f64;
    let nt = targets.len().max(1);
    let pairs = boundary_index_pairs(nurbs, boundary);
    for (k, (idx_b, idx_n)) in pairs.iter().enumerate() {
        if let Some((d, pt)) = compute_shift(
            &nurbs.control_points[*idx_b],
            &nurbs.control_points[*idx_n],
            &targets[k % nt],
            damping,
            max_disp,
        ) {
            worst = worst.max(d);
            nurbs.control_points[*idx_n] = pt;
        }
    }
    worst
}

/// Returns (boundary_idx, next_idx) pairs for the first row of a boundary.
fn boundary_index_pairs(nurbs: &NurbsSurface, boundary: Boundary) -> Vec<(usize, usize)> {
    match boundary {
        Boundary::UMin if nurbs.num_u >= 2 => {
            (0..nurbs.num_v).map(|j| (j, nurbs.num_v + j)).collect()
        }
        Boundary::UMax if nurbs.num_u >= 2 => {
            let lr = nurbs.num_u - 1;
            let pr = nurbs.num_u - 2;
            (0..nurbs.num_v)
                .map(|j| (lr * nurbs.num_v + j, pr * nurbs.num_v + j))
                .collect()
        }
        Boundary::VMin if nurbs.num_v >= 2 => (0..nurbs.num_u)
            .map(|i| (i * nurbs.num_v, i * nurbs.num_v + 1))
            .collect(),
        Boundary::VMax if nurbs.num_v >= 2 => {
            let lc = nurbs.num_v - 1;
            let pc = nurbs.num_v - 2;
            (0..nurbs.num_u)
                .map(|i| (i * nurbs.num_v + lc, i * nurbs.num_v + pc))
                .collect()
        }
        _ => Vec::new(),
    }
}

/// Compute a clamped shift for a next-to-boundary control point.
fn compute_shift(
    boundary_pt: &Point3,
    next_pt: &Point3,
    target_normal: &Option<Vec3>,
    damping: f64,
    max_disp: f64,
) -> Option<(f64, Point3)> {
    let target = (*target_normal)?;
    let diff = *next_pt - *boundary_pt;
    if diff.length() < 1e-15 {
        return None;
    }
    let nc = diff.dot(&target);
    let desired = Vec3::new(target.x * nc, target.y * nc, target.z * nc);
    let new_diff = Vec3::new(
        diff.x + (desired.x - diff.x) * damping,
        diff.y + (desired.y - diff.y) * damping,
        diff.z + (desired.z - diff.z) * damping,
    );
    let new_pt = Point3::new(
        boundary_pt.x + new_diff.x,
        boundary_pt.y + new_diff.y,
        boundary_pt.z + new_diff.z,
    );
    let disp = next_pt.distance_to(&new_pt);
    if disp > max_disp {
        None
    } else {
        Some((disp, new_pt))
    }
}

/// Shift the second row for G2 curvature smoothing.
fn shift_second_row(
    nurbs: &mut NurbsSurface,
    boundary: Boundary,
    damping: f64,
    max_disp: f64,
) -> f64 {
    let mut worst = 0.0_f64;
    let triples = second_row_triples(nurbs, boundary);
    for (idx_a, idx_b, idx_c) in triples {
        let pa = nurbs.control_points[idx_a];
        let pb = nurbs.control_points[idx_b];
        let pc = nurbs.control_points[idx_c];
        // Linear extrapolation: zero second derivative target.
        let target = Point3::new(2.0 * pb.x - pa.x, 2.0 * pb.y - pa.y, 2.0 * pb.z - pa.z);
        let new_c = Point3::new(
            pc.x + (target.x - pc.x) * damping,
            pc.y + (target.y - pc.y) * damping,
            pc.z + (target.z - pc.z) * damping,
        );
        let d = pc.distance_to(&new_c);
        if d <= max_disp {
            worst = worst.max(d);
            nurbs.control_points[idx_c] = new_c;
        }
    }
    worst
}

/// Adjust NURBS control points on adjacent faces for G3 continuity.
///
/// Targets the third row of control points to smooth the rate of
/// curvature change across the boundary. Uses quadratic extrapolation
/// from rows 0-2 to compute target positions for row 3, then shifts
/// toward the target with damping.
pub(super) fn adjust_for_g3(
    brep: &mut BRep,
    face1: FaceHandle,
    face2: FaceHandle,
    samples: &[(f64, f64, f64, f64)],
    max_displacement: f64,
) -> f64 {
    let damping = 0.3_f64;
    let mut worst = 0.0_f64;

    for &fh in &[face1, face2] {
        let nurbs = match extract_nurbs_g3(brep, fh) {
            Some(n) => n,
            None => continue,
        };
        let sample_uv = if fh == face1 {
            (samples[0].0, samples[0].1)
        } else {
            (samples[0].2, samples[0].3)
        };
        let ((u_min, u_max), (v_min, v_max)) = nurbs.parameter_range();
        let boundary = identify_boundary(sample_uv, u_min, u_max, v_min, v_max);
        let mut nc = nurbs;
        let d = shift_third_row(&mut nc, boundary, damping, max_displacement);
        worst = worst.max(d);
        write_nurbs(brep, fh, nc);
    }
    worst
}

/// Extract a cloned NURBS surface from a face, requiring degree >= 5.
fn extract_nurbs_g3(brep: &BRep, fh: FaceHandle) -> Option<NurbsSurface> {
    let n = extract_nurbs(brep, fh)?;
    if n.degree_u >= 5 && n.degree_v >= 5 {
        Some(n)
    } else {
        None
    }
}

/// Shift the third row for G3 curvature-rate smoothing.
///
/// Computes a target for each third-row control point via quadratic
/// extrapolation from rows 0, 1, and 2:
///   target = 3 * P[2] - 3 * P[1] + P[0]
/// This is the cubic extension that produces a zero third-derivative
/// discontinuity at the boundary.
fn shift_third_row(
    nurbs: &mut NurbsSurface,
    boundary: Boundary,
    damping: f64,
    max_disp: f64,
) -> f64 {
    let mut worst = 0.0_f64;
    let quads = third_row_quads(nurbs, boundary);
    for (idx_a, idx_b, idx_c, idx_d) in quads {
        let pa = nurbs.control_points[idx_a];
        let pb = nurbs.control_points[idx_b];
        let pc = nurbs.control_points[idx_c];
        let pd = nurbs.control_points[idx_d];
        // Quadratic extrapolation: target = 3*P[2] - 3*P[1] + P[0].
        let target = Point3::new(
            3.0 * pc.x - 3.0 * pb.x + pa.x,
            3.0 * pc.y - 3.0 * pb.y + pa.y,
            3.0 * pc.z - 3.0 * pb.z + pa.z,
        );
        let new_d = Point3::new(
            pd.x + (target.x - pd.x) * damping,
            pd.y + (target.y - pd.y) * damping,
            pd.z + (target.z - pd.z) * damping,
        );
        let d = pd.distance_to(&new_d);
        if d <= max_disp {
            worst = worst.max(d);
            nurbs.control_points[idx_d] = new_d;
        }
    }
    worst
}

/// Returns (boundary, next, second, third) index quads for a boundary.
///
/// Indices are: row 0 (boundary), row 1, row 2, row 3 (target).
fn third_row_quads(nurbs: &NurbsSurface, boundary: Boundary) -> Vec<(usize, usize, usize, usize)> {
    let nv = nurbs.num_v;
    match boundary {
        Boundary::UMin if nurbs.num_u >= 4 => (0..nv)
            .map(|j| (j, nv + j, 2 * nv + j, 3 * nv + j))
            .collect(),
        Boundary::UMax if nurbs.num_u >= 4 => {
            let r = nurbs.num_u - 1;
            (0..nv)
                .map(|j| {
                    (
                        r * nv + j,
                        (r - 1) * nv + j,
                        (r - 2) * nv + j,
                        (r - 3) * nv + j,
                    )
                })
                .collect()
        }
        Boundary::VMin if nurbs.num_v >= 4 => (0..nurbs.num_u)
            .map(|i| (i * nv, i * nv + 1, i * nv + 2, i * nv + 3))
            .collect(),
        Boundary::VMax if nurbs.num_v >= 4 => {
            let c = nv - 1;
            (0..nurbs.num_u)
                .map(|i| (i * nv + c, i * nv + c - 1, i * nv + c - 2, i * nv + c - 3))
                .collect()
        }
        _ => Vec::new(),
    }
}

/// Returns (boundary, next, second) index triples for a boundary.
fn second_row_triples(nurbs: &NurbsSurface, boundary: Boundary) -> Vec<(usize, usize, usize)> {
    let nv = nurbs.num_v;
    match boundary {
        Boundary::UMin if nurbs.num_u >= 3 => (0..nv).map(|j| (j, nv + j, 2 * nv + j)).collect(),
        Boundary::UMax if nurbs.num_u >= 3 => {
            let r = nurbs.num_u - 1;
            (0..nv)
                .map(|j| (r * nv + j, (r - 1) * nv + j, (r - 2) * nv + j))
                .collect()
        }
        Boundary::VMin if nurbs.num_v >= 3 => (0..nurbs.num_u)
            .map(|i| (i * nv, i * nv + 1, i * nv + 2))
            .collect(),
        Boundary::VMax if nurbs.num_v >= 3 => {
            let c = nv - 1;
            (0..nurbs.num_u)
                .map(|i| (i * nv + c, i * nv + c - 1, i * nv + c - 2))
                .collect()
        }
        _ => Vec::new(),
    }
}
