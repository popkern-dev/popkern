//! Adjacent tangent-row alignment for KP20 support-cell shared boundaries.

use std::collections::HashMap;

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, FaceHandle};

use crate::construct::kp20_polar_support_cell_shared_trim_quality as shared_trim_quality;

use super::{
    BoundarySide, FullBoundaryTrim, boundary_row_group, rectangular_binonic_on_face,
    set_rectangular_binonic_on_face, trim_ownership,
};

const TANGENT_ROW_DAMPING_FACTORS: [f64; 4] = [1.0, 0.5, 0.25, 0.125];

#[derive(Debug, Clone, Copy)]
pub(super) enum TangentRowAlignment {
    SharedDirection,
    MirroredDirection,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct ControlKey {
    face: FaceHandle,
    u: usize,
    v: usize,
}

#[derive(Debug, Clone, Copy)]
struct TargetAccumulator {
    point: Point3,
    count: usize,
}

fn align_worst_row_group_tangent_band_with_damping(
    brep: &mut BRep,
    trims: &[FullBoundaryTrim],
    alignment: TangentRowAlignment,
    damping: f64,
) -> Result<(), CoreError> {
    let row_group = boundary_row_group::worst_boundary_row_group_trims(brep, trims)?;
    let mut targets = HashMap::new();
    for trim in row_group {
        push_trim_targets(brep, trim, alignment, &mut targets)?;
    }
    apply_targets(brep, &targets, damping)
}

pub(super) fn reconcile_worst_row_group_if_normal_improves(
    brep: &mut BRep,
    faces: &[FaceHandle],
    trims: &[FullBoundaryTrim],
) -> Result<Option<TangentRowAlignment>, CoreError> {
    reconcile_tangent_band_if_normal_improves(brep, faces, |candidate, alignment, damping| {
        align_worst_row_group_tangent_band_with_damping(candidate, trims, alignment, damping)?;
        Ok(true)
    })
}

fn align_segment_owned_tangent_band_with_damping(
    brep: &mut BRep,
    trims: &[FullBoundaryTrim],
    alignment: TangentRowAlignment,
    damping: f64,
) -> Result<bool, CoreError> {
    let segment_trims = trim_ownership::segment_owned_boundary_trims(trims);
    if segment_trims.is_empty() {
        return Ok(false);
    }
    let mut targets = HashMap::new();
    for record in segment_trims {
        push_trim_targets(brep, record.trim, alignment, &mut targets)?;
    }
    apply_targets(brep, &targets, damping)?;
    Ok(true)
}

pub(super) fn reconcile_segment_owned_if_normal_improves(
    brep: &mut BRep,
    faces: &[FaceHandle],
    trims: &[FullBoundaryTrim],
) -> Result<Option<TangentRowAlignment>, CoreError> {
    reconcile_tangent_band_if_normal_improves(brep, faces, |candidate, alignment, damping| {
        align_segment_owned_tangent_band_with_damping(candidate, trims, alignment, damping)
    })
}

fn reconcile_tangent_band_if_normal_improves(
    brep: &mut BRep,
    faces: &[FaceHandle],
    align: impl Fn(&mut BRep, TangentRowAlignment, f64) -> Result<bool, CoreError>,
) -> Result<Option<TangentRowAlignment>, CoreError> {
    let before = shared_trim_quality::measure_support_cell_shared_trim_quality(brep, faces, 9);
    let mut best: Option<(
        TangentRowAlignment,
        shared_trim_quality::Kp20SupportCellSharedTrimQuality,
        BRep,
    )> = None;
    for alignment in [
        TangentRowAlignment::SharedDirection,
        TangentRowAlignment::MirroredDirection,
    ] {
        for damping in TANGENT_ROW_DAMPING_FACTORS {
            let mut candidate = brep.clone();
            if !align(&mut candidate, alignment, damping)? || candidate.validate().is_err() {
                continue;
            }
            let after =
                shared_trim_quality::measure_support_cell_shared_trim_quality(&candidate, faces, 9);
            if after.worst_g0_gap > before.worst_g0_gap + 1.0e-12
                || after.worst_normal_disagreement + 1.0e-12 >= before.worst_normal_disagreement
            {
                continue;
            }
            let improves_best = match best.as_ref() {
                Some((_, quality, _)) => {
                    after.worst_normal_disagreement < quality.worst_normal_disagreement
                }
                None => true,
            };
            if improves_best {
                best = Some((alignment, after, candidate));
            }
        }
    }
    let Some((alignment, _, candidate)) = best else {
        return Ok(None);
    };
    *brep = candidate;
    Ok(Some(alignment))
}

fn push_trim_targets(
    brep: &BRep,
    trim: FullBoundaryTrim,
    alignment: TangentRowAlignment,
    targets: &mut HashMap<ControlKey, TargetAccumulator>,
) -> Result<(), CoreError> {
    for index in 0..10 {
        let fraction = index as f64 / 9.0;
        let key_a = adjacent_key(trim.face_a, trim.side_a, trim.start_a, trim.end_a, fraction);
        let key_b = adjacent_key(
            trim.face_b,
            trim.side_b,
            trim.start_b,
            trim.end_b,
            1.0 - fraction,
        );
        let boundary_a = control_point(brep, boundary_key(key_a, trim.side_a))?;
        let boundary_b = control_point(brep, boundary_key(key_b, trim.side_b))?;
        let dir_a = control_point(brep, key_a)? - boundary_a;
        let dir_b = control_point(brep, key_b)? - boundary_b;
        let (target_a, target_b) = match alignment {
            TangentRowAlignment::SharedDirection => {
                let direction = (dir_a + dir_b) * 0.5;
                (boundary_a + direction, boundary_b + direction)
            }
            TangentRowAlignment::MirroredDirection => {
                let direction = (dir_a - dir_b) * 0.5;
                (boundary_a + direction, boundary_b - direction)
            }
        };
        push_target(targets, key_a, target_a);
        push_target(targets, key_b, target_b);
    }
    Ok(())
}

fn adjacent_key(
    face: FaceHandle,
    side: BoundarySide,
    start: f64,
    end: f64,
    fraction: f64,
) -> ControlKey {
    let index = control_index(start + (end - start) * fraction);
    let (u, v) = match side {
        BoundarySide::UMin => (1, index),
        BoundarySide::UMax => (8, index),
        BoundarySide::VMin => (index, 1),
        BoundarySide::VMax => (index, 8),
    };
    ControlKey { face, u, v }
}

fn boundary_key(key: ControlKey, side: BoundarySide) -> ControlKey {
    let (u, v) = match side {
        BoundarySide::UMin => (0, key.v),
        BoundarySide::UMax => (9, key.v),
        BoundarySide::VMin => (key.u, 0),
        BoundarySide::VMax => (key.u, 9),
    };
    ControlKey {
        face: key.face,
        u,
        v,
    }
}

fn control_index(parameter: f64) -> usize {
    (parameter.clamp(0.0, 1.0) * 9.0).round() as usize
}

fn control_point(brep: &BRep, key: ControlKey) -> Result<Point3, CoreError> {
    Ok(rectangular_binonic_on_face(brep, key.face)?.grid[key.u][key.v])
}

fn push_target(
    targets: &mut HashMap<ControlKey, TargetAccumulator>,
    key: ControlKey,
    point: Point3,
) {
    targets
        .entry(key)
        .and_modify(|target| {
            target.point =
                target.point + (point - target.point) * (1.0 / (target.count as f64 + 1.0));
            target.count += 1;
        })
        .or_insert(TargetAccumulator { point, count: 1 });
}

fn apply_targets(
    brep: &mut BRep,
    targets: &HashMap<ControlKey, TargetAccumulator>,
    damping: f64,
) -> Result<(), CoreError> {
    for (&key, target) in targets {
        let mut patch = rectangular_binonic_on_face(brep, key.face)?;
        let point = patch.grid[key.u][key.v];
        patch.grid[key.u][key.v] = point + (target.point - point) * damping;
        set_rectangular_binonic_on_face(brep, key.face, patch)?;
    }
    Ok(())
}
