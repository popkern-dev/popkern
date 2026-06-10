//! Exact shared boundary-row ownership for KP20 support-cell patches.

use std::collections::HashMap;

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::topology::{BRep, FaceHandle};

use super::{
    BoundarySide, FullBoundaryTrim, full_boundary_trims, rectangular_binonic_on_face, row_fit,
    set_rectangular_binonic_on_face, trim_ownership,
};

const PARTIAL_SEGMENT_SAMPLE_COUNT: usize = 16;
const PARTIAL_SEGMENT_ANCHOR_SAMPLE_COUNT: usize = 25;
const PARTIAL_SEGMENT_ANCHOR_REPEAT_COUNT: usize = 4;

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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct BoundaryRowKey {
    face: FaceHandle,
    side: BoundarySide,
}

pub(super) fn exact_full_span_shared_boundary_rows(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<&'static str, CoreError> {
    let trims = trim_ownership::full_row_owned_trims(&full_boundary_trims(brep, slots)?);
    let mut targets = HashMap::new();
    for trim in trims {
        push_trim_targets(brep, trim, &mut targets)?;
    }
    if targets.is_empty() {
        return Ok("exact_boundary_no_full_span_trims");
    }
    apply_targets(brep, &targets)?;
    Ok("exact_boundary_full_span_rows_authored")
}

pub(super) fn fit_partial_shared_boundary_segments(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<&'static str, CoreError> {
    let segment_records =
        trim_ownership::segment_owned_boundary_trims(&full_boundary_trims(brep, slots)?);
    if segment_records
        .iter()
        .any(|record| !record.has_topology_address())
    {
        return Ok("partial_boundary_segment_topology_address_invalid");
    }
    let mut samples = HashMap::new();
    let mut edited_intervals = HashMap::new();
    for record in segment_records {
        let trim = record.trim;
        push_partial_trim_samples(brep, trim, &mut samples, &mut edited_intervals)?;
    }
    if samples.is_empty() {
        return Ok("partial_boundary_no_partial_trims");
    }
    push_original_boundary_anchor_samples(brep, &mut samples, &edited_intervals)?;
    fit_boundary_row_samples(brep, &samples)?;
    Ok("bounded_partial_boundary_segments_fit")
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

fn push_partial_trim_samples(
    brep: &BRep,
    trim: FullBoundaryTrim,
    samples: &mut HashMap<BoundaryRowKey, Vec<row_fit::BoundaryRowFitSample>>,
    edited_intervals: &mut HashMap<BoundaryRowKey, Vec<(f64, f64)>>,
) -> Result<(), CoreError> {
    let patch_a = rectangular_binonic_on_face(brep, trim.face_a)?;
    let patch_b = rectangular_binonic_on_face(brep, trim.face_b)?;
    push_interval(
        edited_intervals,
        trim.face_a,
        trim.side_a,
        trim.start_a,
        trim.end_a,
    );
    push_interval(
        edited_intervals,
        trim.face_b,
        trim.side_b,
        trim.start_b,
        trim.end_b,
    );
    for index in 0..PARTIAL_SEGMENT_SAMPLE_COUNT {
        let fraction = index as f64 / (PARTIAL_SEGMENT_SAMPLE_COUNT - 1) as f64;
        let param_a = lerp_scalar(trim.start_a, trim.end_a, fraction);
        let param_b = lerp_scalar(trim.start_b, trim.end_b, 1.0 - fraction);
        let point_a = row_fit::eval_boundary(&patch_a, trim.side_a, param_a);
        let point_b = row_fit::eval_boundary(&patch_b, trim.side_b, param_b);
        let target = point_a + (point_b - point_a) * 0.5;
        push_boundary_sample(samples, trim.face_a, trim.side_a, param_a, target);
        push_boundary_sample(samples, trim.face_b, trim.side_b, param_b, target);
    }
    Ok(())
}

fn push_original_boundary_anchor_samples(
    brep: &BRep,
    samples: &mut HashMap<BoundaryRowKey, Vec<row_fit::BoundaryRowFitSample>>,
    edited_intervals: &HashMap<BoundaryRowKey, Vec<(f64, f64)>>,
) -> Result<(), CoreError> {
    for (key, intervals) in edited_intervals {
        let patch = rectangular_binonic_on_face(brep, key.face)?;
        for index in 0..PARTIAL_SEGMENT_ANCHOR_SAMPLE_COUNT {
            let parameter = index as f64 / (PARTIAL_SEGMENT_ANCHOR_SAMPLE_COUNT - 1) as f64;
            if intervals
                .iter()
                .any(|&(start, end)| parameter >= start - 1.0e-9 && parameter <= end + 1.0e-9)
            {
                continue;
            }
            let target = row_fit::eval_boundary(&patch, key.side, parameter);
            for _ in 0..PARTIAL_SEGMENT_ANCHOR_REPEAT_COUNT {
                push_boundary_sample(samples, key.face, key.side, parameter, target);
            }
        }
    }
    Ok(())
}

fn push_interval(
    edited_intervals: &mut HashMap<BoundaryRowKey, Vec<(f64, f64)>>,
    face: FaceHandle,
    side: BoundarySide,
    start: f64,
    end: f64,
) {
    edited_intervals
        .entry(BoundaryRowKey { face, side })
        .or_default()
        .push(if start <= end {
            (start, end)
        } else {
            (end, start)
        });
}

fn natural_index_for_fraction(start: f64, end: f64, fraction_index: usize) -> usize {
    if start <= end {
        fraction_index
    } else {
        9 - fraction_index
    }
}

fn lerp_scalar(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
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

fn push_boundary_sample(
    samples: &mut HashMap<BoundaryRowKey, Vec<row_fit::BoundaryRowFitSample>>,
    face: FaceHandle,
    side: BoundarySide,
    parameter: f64,
    target: Point3,
) {
    samples
        .entry(BoundaryRowKey { face, side })
        .or_default()
        .push(row_fit::BoundaryRowFitSample { parameter, target });
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

fn fit_boundary_row_samples(
    brep: &mut BRep,
    samples: &HashMap<BoundaryRowKey, Vec<row_fit::BoundaryRowFitSample>>,
) -> Result<(), CoreError> {
    let mut patches = HashMap::new();
    for key in samples.keys() {
        patches
            .entry(key.face)
            .or_insert(rectangular_binonic_on_face(brep, key.face)?);
    }
    for (key, row_samples) in samples {
        let Some(patch) = patches.get_mut(&key.face) else {
            continue;
        };
        row_fit::fit_boundary_samples(patch, key.side, row_samples);
    }
    for (face, patch) in patches {
        set_rectangular_binonic_on_face(brep, face, patch)?;
    }
    Ok(())
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
