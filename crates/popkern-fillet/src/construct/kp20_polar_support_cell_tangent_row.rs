//! Adjacent tangent-row alignment for KP20 support-cell shared boundaries.

use std::collections::HashMap;

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, FaceHandle};

use super::{
    BoundarySide, FullBoundaryTrim, boundary_row_group, rectangular_binonic_on_face,
    set_rectangular_binonic_on_face,
};

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

pub(super) fn align_worst_row_group_tangent_band(
    brep: &mut BRep,
    trims: &[FullBoundaryTrim],
    alignment: TangentRowAlignment,
) -> Result<(), CoreError> {
    let row_group = boundary_row_group::worst_boundary_row_group_trims(brep, trims)?;
    let mut targets = HashMap::new();
    for trim in row_group {
        push_trim_targets(brep, trim, alignment, &mut targets)?;
    }
    apply_targets(brep, &targets)
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
) -> Result<(), CoreError> {
    for (&key, target) in targets {
        let mut patch = rectangular_binonic_on_face(brep, key.face)?;
        patch.grid[key.u][key.v] = target.point;
        set_rectangular_binonic_on_face(brep, key.face, patch)?;
    }
    Ok(())
}
