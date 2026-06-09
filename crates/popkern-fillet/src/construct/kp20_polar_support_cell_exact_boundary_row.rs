//! Exact shared boundary-row ownership for KP20 support-cell patches.

use std::collections::HashMap;

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::topology::{BRep, FaceHandle};

use super::{
    BoundarySide, FullBoundaryTrim, full_boundary_trims, rectangular_binonic_on_face,
    set_rectangular_binonic_on_face,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct BoundaryControlKey {
    face: FaceHandle,
    side: BoundarySide,
    index: usize,
}

#[derive(Debug, Clone, Copy)]
struct BoundaryTarget {
    point: Point3,
    count: usize,
}

pub(super) fn exact_full_span_shared_boundary_rows(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<&'static str, CoreError> {
    let trims = full_boundary_trims(brep, slots)?;
    let mut targets = HashMap::new();
    let mut compatible_trim_count = 0;
    for trim in trims {
        if !full_span(trim.start_a, trim.end_a) || !full_span(trim.start_b, trim.end_b) {
            continue;
        }
        push_trim_targets(brep, trim, &mut targets)?;
        compatible_trim_count += 1;
    }
    if compatible_trim_count == 0 {
        return Ok("exact_boundary_no_full_span_trims");
    }
    apply_targets(brep, &targets)?;
    Ok("exact_boundary_full_span_rows_authored")
}

fn push_trim_targets(
    brep: &BRep,
    trim: FullBoundaryTrim,
    targets: &mut HashMap<BoundaryControlKey, BoundaryTarget>,
) -> Result<(), CoreError> {
    let patch_a = rectangular_binonic_on_face(brep, trim.face_a)?;
    let patch_b = rectangular_binonic_on_face(brep, trim.face_b)?;
    let controls_a = boundary_controls(&patch_a, trim.side_a);
    let controls_b = boundary_controls(&patch_b, trim.side_b);
    for fraction_index in 0..10 {
        let point_a =
            controls_a[natural_index_for_fraction(trim.start_a, trim.end_a, fraction_index)];
        let point_b =
            controls_b[natural_index_for_fraction(trim.start_b, trim.end_b, 9 - fraction_index)];
        let target = point_a + (point_b - point_a) * 0.5;
        push_target(
            targets,
            BoundaryControlKey {
                face: trim.face_a,
                side: trim.side_a,
                index: natural_index_for_fraction(trim.start_a, trim.end_a, fraction_index),
            },
            target,
        );
        push_target(
            targets,
            BoundaryControlKey {
                face: trim.face_b,
                side: trim.side_b,
                index: natural_index_for_fraction(trim.start_b, trim.end_b, 9 - fraction_index),
            },
            target,
        );
    }
    Ok(())
}

fn full_span(start: f64, end: f64) -> bool {
    (close(start, 0.0) && close(end, 1.0)) || (close(start, 1.0) && close(end, 0.0))
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-9
}

fn natural_index_for_fraction(start: f64, end: f64, fraction_index: usize) -> usize {
    if start <= end {
        fraction_index
    } else {
        9 - fraction_index
    }
}

fn boundary_controls(patch: &RectangularGregoryBiNonic, side: BoundarySide) -> [Point3; 10] {
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in out.iter_mut().enumerate() {
        *point = match side {
            BoundarySide::UMin => patch.grid[0][index],
            BoundarySide::UMax => patch.grid[9][index],
            BoundarySide::VMin => patch.grid[index][0],
            BoundarySide::VMax => patch.grid[index][9],
        };
    }
    out
}

fn push_target(
    targets: &mut HashMap<BoundaryControlKey, BoundaryTarget>,
    key: BoundaryControlKey,
    point: Point3,
) {
    targets
        .entry(key)
        .and_modify(|target| {
            target.point =
                target.point + (point - target.point) * (1.0 / (target.count as f64 + 1.0));
            target.count += 1;
        })
        .or_insert(BoundaryTarget { point, count: 1 });
}

fn apply_targets(
    brep: &mut BRep,
    targets: &HashMap<BoundaryControlKey, BoundaryTarget>,
) -> Result<(), CoreError> {
    let mut patches = HashMap::new();
    for key in targets.keys() {
        patches
            .entry(key.face)
            .or_insert(rectangular_binonic_on_face(brep, key.face)?);
    }
    for (key, target) in targets {
        let Some(patch) = patches.get_mut(&key.face) else {
            continue;
        };
        set_boundary_control(patch, key.side, key.index, target.point);
    }
    for (face, patch) in patches {
        set_rectangular_binonic_on_face(brep, face, patch)?;
    }
    Ok(())
}

fn set_boundary_control(
    patch: &mut RectangularGregoryBiNonic,
    side: BoundarySide,
    index: usize,
    point: Point3,
) {
    match side {
        BoundarySide::UMin => patch.grid[0][index] = point,
        BoundarySide::UMax => patch.grid[9][index] = point,
        BoundarySide::VMin => patch.grid[index][0] = point,
        BoundarySide::VMax => patch.grid[index][9] = point,
    }
}
