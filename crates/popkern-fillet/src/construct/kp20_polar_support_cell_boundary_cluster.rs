//! Coupled worst-neighborhood fitting for KP20 support-cell boundaries.

use std::collections::HashMap;

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
    let mut best = endpoint_cluster_candidate(brep, trims, &cluster, false, before_gap)?;
    if let Some(candidate) = endpoint_cluster_candidate(brep, trims, &cluster, true, before_gap)? {
        let improves_best = best
            .as_ref()
            .map(|(_, best_gap)| candidate.1 < *best_gap)
            .unwrap_or(true);
        if improves_best {
            best = Some(candidate);
        }
    }
    if let Some((candidate, _)) = best {
        *brep = candidate;
    }
    Ok(())
}

fn endpoint_cluster_candidate(
    brep: &BRep,
    trims: &[FullBoundaryTrim],
    cluster: &[PatchCornerKey],
    include_adjacent_controls: bool,
    before_gap: f64,
) -> Result<Option<(BRep, f64)>, CoreError> {
    let mut candidate = brep.clone();
    average_patch_corners(&mut candidate, cluster)?;
    if include_adjacent_controls {
        average_adjacent_boundary_controls(&mut candidate, cluster, trims)?;
    }
    if candidate.validate().is_err() {
        return Ok(None);
    }
    let after_gap = worst_trim_gap(&candidate, trims)?;
    if after_gap + 1.0e-12 < before_gap {
        Ok(Some((candidate, after_gap)))
    } else {
        Ok(None)
    }
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct PatchControlKey {
    face: FaceHandle,
    u: usize,
    v: usize,
}

#[derive(Debug, Clone, Copy)]
struct ControlAccumulator {
    x: f64,
    y: f64,
    z: f64,
    count: usize,
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

fn average_adjacent_boundary_controls(
    brep: &mut BRep,
    cluster: &[PatchCornerKey],
    trims: &[FullBoundaryTrim],
) -> Result<(), CoreError> {
    let mut targets: HashMap<PatchControlKey, ControlAccumulator> = HashMap::new();
    for &trim in trims {
        for (corner_a, corner_b) in trim_endpoint_pairs(trim) {
            if !cluster.contains(&corner_a) || !cluster.contains(&corner_b) {
                continue;
            }
            let Some(control_a) = adjacent_boundary_control(trim.face_a, trim.side_a, corner_a)
            else {
                continue;
            };
            let Some(control_b) = adjacent_boundary_control(trim.face_b, trim.side_b, corner_b)
            else {
                continue;
            };
            let point_a = control_point(brep, control_a)?;
            let point_b = control_point(brep, control_b)?;
            let target = average_point(point_a, point_b);
            push_target(&mut targets, control_a, target);
            push_target(&mut targets, control_b, target);
        }
    }
    apply_control_targets(brep, &targets)
}

fn adjacent_boundary_control(
    face: FaceHandle,
    side: BoundarySide,
    corner: PatchCornerKey,
) -> Option<PatchControlKey> {
    if corner.face != face {
        return None;
    }
    let u = usize::from(corner.u_max) * 9;
    let v = usize::from(corner.v_max) * 9;
    let key = match side {
        BoundarySide::UMin | BoundarySide::UMax => PatchControlKey {
            face,
            u,
            v: adjacent_index(v)?,
        },
        BoundarySide::VMin | BoundarySide::VMax => PatchControlKey {
            face,
            u: adjacent_index(u)?,
            v,
        },
    };
    Some(key)
}

fn adjacent_index(index: usize) -> Option<usize> {
    match index {
        0 => Some(1),
        9 => Some(8),
        _ => None,
    }
}

fn push_target(
    targets: &mut HashMap<PatchControlKey, ControlAccumulator>,
    key: PatchControlKey,
    point: Point3,
) {
    targets
        .entry(key)
        .and_modify(|target| {
            target.x += point.x;
            target.y += point.y;
            target.z += point.z;
            target.count += 1;
        })
        .or_insert(ControlAccumulator {
            x: point.x,
            y: point.y,
            z: point.z,
            count: 1,
        });
}

fn apply_control_targets(
    brep: &mut BRep,
    targets: &HashMap<PatchControlKey, ControlAccumulator>,
) -> Result<(), CoreError> {
    for (&key, target) in targets {
        let mut patch = rectangular_binonic_on_face(brep, key.face)?;
        let scale = 1.0 / target.count as f64;
        patch.grid[key.u][key.v] =
            Point3::new(target.x * scale, target.y * scale, target.z * scale);
        set_rectangular_binonic_on_face(brep, key.face, patch)?;
    }
    Ok(())
}

fn control_point(brep: &BRep, key: PatchControlKey) -> Result<Point3, CoreError> {
    let patch = rectangular_binonic_on_face(brep, key.face)?;
    Ok(patch.grid[key.u][key.v])
}

fn average_point(a: Point3, b: Point3) -> Point3 {
    Point3::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5)
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
