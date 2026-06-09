//! Coupled worst-neighborhood fitting for KP20 support-cell boundaries.

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, FaceHandle};

use super::{
    BoundarySide, FullBoundaryTrim, group_reconcile, rectangular_binonic_on_face,
    set_rectangular_binonic_on_face, top_boundary_trims, trim_gap,
};

pub(super) fn reconcile_worst_boundary_neighborhood(
    brep: &mut BRep,
    trims: &[FullBoundaryTrim],
) -> Result<(), CoreError> {
    let Some(worst) = top_boundary_trims(brep, trims, 1)?.first().copied() else {
        return Ok(());
    };
    let cluster = worst_neighborhood(trims, worst);
    let before_gap = worst_trim_gap(brep, trims)?;
    let mut candidate = brep.clone();
    for trim in &cluster {
        super::reconcile_full_boundary_trim(&mut candidate, *trim)?;
    }
    group_reconcile::reconcile_grouped_boundary_rows(&mut candidate, &cluster)?;
    if candidate.validate().is_err() {
        return Ok(());
    }
    if worst_trim_gap(&candidate, trims)? + 1.0e-12 < before_gap {
        *brep = candidate;
    }
    Ok(())
}

pub(super) fn reconcile_worst_endpoint_cluster(
    brep: &mut BRep,
    trims: &[FullBoundaryTrim],
) -> Result<(), CoreError> {
    let Some(worst) = worst_trim_endpoint(brep, trims)? else {
        return Ok(());
    };
    let mut cluster = vec![worst.corner_a, worst.corner_b];
    grow_endpoint_cluster(&mut cluster, trims);
    let before_gap = worst_trim_gap(brep, trims)?;
    let mut candidate = brep.clone();
    average_patch_corners(&mut candidate, &cluster)?;
    if candidate.validate().is_err() {
        return Ok(());
    }
    if worst_trim_gap(&candidate, trims)? + 1.0e-12 < before_gap {
        *brep = candidate;
    }
    Ok(())
}

fn worst_neighborhood(
    trims: &[FullBoundaryTrim],
    worst: FullBoundaryTrim,
) -> Vec<FullBoundaryTrim> {
    trims
        .iter()
        .copied()
        .filter(|trim| shares_boundary_row(*trim, worst))
        .collect()
}

fn shares_boundary_row(trim: FullBoundaryTrim, worst: FullBoundaryTrim) -> bool {
    same_row(trim.face_a, trim.side_a, worst.face_a, worst.side_a)
        || same_row(trim.face_a, trim.side_a, worst.face_b, worst.side_b)
        || same_row(trim.face_b, trim.side_b, worst.face_a, worst.side_a)
        || same_row(trim.face_b, trim.side_b, worst.face_b, worst.side_b)
}

fn same_row(
    face_a: popkern_core::topology::FaceHandle,
    side_a: super::BoundarySide,
    face_b: popkern_core::topology::FaceHandle,
    side_b: super::BoundarySide,
) -> bool {
    face_a == face_b && side_a == side_b
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PatchCornerKey {
    face: FaceHandle,
    u_max: bool,
    v_max: bool,
}

#[derive(Debug, Clone, Copy)]
struct WorstEndpoint {
    corner_a: PatchCornerKey,
    corner_b: PatchCornerKey,
}

fn worst_trim_endpoint(
    brep: &BRep,
    trims: &[FullBoundaryTrim],
) -> Result<Option<WorstEndpoint>, CoreError> {
    let mut worst = None;
    let mut worst_gap = 0.0_f64;
    for &trim in trims {
        let patch_a = rectangular_binonic_on_face(brep, trim.face_a)?;
        let patch_b = rectangular_binonic_on_face(brep, trim.face_b)?;
        for index in 0..9 {
            let fraction = index as f64 / 8.0;
            let param_a = lerp_scalar(trim.start_a, trim.end_a, fraction);
            let param_b = lerp_scalar(trim.start_b, trim.end_b, 1.0 - fraction);
            let gap = super::row_fit::eval_boundary(&patch_a, trim.side_a, param_a).distance_to(
                &super::row_fit::eval_boundary(&patch_b, trim.side_b, param_b),
            );
            if gap > worst_gap {
                worst_gap = gap;
                worst = corner_pair_at(trim, param_a, param_b);
            }
        }
    }
    Ok(worst)
}

fn corner_pair_at(trim: FullBoundaryTrim, param_a: f64, param_b: f64) -> Option<WorstEndpoint> {
    Some(WorstEndpoint {
        corner_a: patch_corner(trim.face_a, trim.side_a, param_a)?,
        corner_b: patch_corner(trim.face_b, trim.side_b, param_b)?,
    })
}

fn grow_endpoint_cluster(cluster: &mut Vec<PatchCornerKey>, trims: &[FullBoundaryTrim]) {
    let mut changed = true;
    while changed {
        changed = false;
        for &trim in trims {
            for (corner_a, corner_b) in trim_endpoint_pairs(trim) {
                if cluster.contains(&corner_a) && !cluster.contains(&corner_b) {
                    cluster.push(corner_b);
                    changed = true;
                }
                if cluster.contains(&corner_b) && !cluster.contains(&corner_a) {
                    cluster.push(corner_a);
                    changed = true;
                }
            }
        }
    }
}

fn trim_endpoint_pairs(trim: FullBoundaryTrim) -> Vec<(PatchCornerKey, PatchCornerKey)> {
    [(trim.start_a, trim.end_b), (trim.end_a, trim.start_b)]
        .into_iter()
        .filter_map(|(param_a, param_b)| {
            Some((
                patch_corner(trim.face_a, trim.side_a, param_a)?,
                patch_corner(trim.face_b, trim.side_b, param_b)?,
            ))
        })
        .collect()
}

fn patch_corner(face: FaceHandle, side: BoundarySide, parameter: f64) -> Option<PatchCornerKey> {
    let at_min = close(parameter, 0.0);
    let at_max = close(parameter, 1.0);
    if !at_min && !at_max {
        return None;
    }
    let (u_max, v_max) = match side {
        BoundarySide::UMin => (false, at_max),
        BoundarySide::UMax => (true, at_max),
        BoundarySide::VMin => (at_max, false),
        BoundarySide::VMax => (at_max, true),
    };
    Some(PatchCornerKey { face, u_max, v_max })
}

fn average_patch_corners(brep: &mut BRep, cluster: &[PatchCornerKey]) -> Result<(), CoreError> {
    if cluster.len() < 2 {
        return Ok(());
    }
    let mut sum = Point3::new(0.0, 0.0, 0.0);
    for key in cluster {
        let patch = rectangular_binonic_on_face(brep, key.face)?;
        let point = corner_point(&patch, *key);
        sum.x += point.x;
        sum.y += point.y;
        sum.z += point.z;
    }
    let scale = 1.0 / cluster.len() as f64;
    let average = Point3::new(sum.x * scale, sum.y * scale, sum.z * scale);
    for key in cluster {
        let mut patch = rectangular_binonic_on_face(brep, key.face)?;
        set_corner_point(&mut patch, *key, average);
        set_rectangular_binonic_on_face(brep, key.face, patch)?;
    }
    Ok(())
}

fn corner_point(
    patch: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    key: PatchCornerKey,
) -> Point3 {
    patch.grid[usize::from(key.u_max) * 9][usize::from(key.v_max) * 9]
}

fn set_corner_point(
    patch: &mut popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    key: PatchCornerKey,
    point: Point3,
) {
    patch.grid[usize::from(key.u_max) * 9][usize::from(key.v_max) * 9] = point;
}

fn worst_trim_gap(brep: &BRep, trims: &[FullBoundaryTrim]) -> Result<f64, CoreError> {
    let mut worst = 0.0_f64;
    for &trim in trims {
        let patch_a = rectangular_binonic_on_face(brep, trim.face_a)?;
        let patch_b = rectangular_binonic_on_face(brep, trim.face_b)?;
        worst = worst.max(trim_gap(&patch_a, &patch_b, trim));
    }
    Ok(worst)
}

fn lerp_scalar(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-9
}
