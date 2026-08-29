//! Exact shared boundary-row ownership for KP20 support-cell patches.

use std::collections::HashMap;

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::topology::{BRep, FaceHandle};

use super::row_segment::bezier_row_segment;
use super::{
    BoundarySide, FullBoundaryTrim, full_boundary_trims, rectangular_binonic_on_face, row_fit,
    segment_topology, set_rectangular_binonic_on_face, trim_ownership,
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
    let mut classes = SharedControlClasses::default();
    for trim in &trims {
        push_trim_control_classes(*trim, &mut classes);
    }
    let targets = classes.class_average_targets(brep)?;
    if targets.is_empty() {
        return Ok("exact_boundary_no_full_span_trims");
    }
    apply_targets(brep, &targets)?;
    Ok("exact_boundary_full_span_rows_authored")
}

/// Equivalence classes of boundary controls that must share one
/// position for every full-span trim to be a CP identity: the 10
/// paired controls along each trim union with the corner aliases
/// within one face (keys normalized to one canonical side). Each
/// class takes its members' average (sorted accumulation).
#[derive(Default)]
struct SharedControlClasses {
    parents: HashMap<BoundaryControlKey, BoundaryControlKey>,
}

impl SharedControlClasses {
    fn find(&mut self, key: BoundaryControlKey) -> BoundaryControlKey {
        let parent = *self.parents.entry(key).or_insert(key);
        if parent == key {
            return key;
        }
        let root = self.find(parent);
        self.parents.insert(key, root);
        root
    }

    fn union(&mut self, a: BoundaryControlKey, b: BoundaryControlKey) {
        let root_a = self.find(a);
        let root_b = self.find(b);
        if root_a != root_b {
            self.parents.insert(root_b, root_a);
        }
    }

    fn class_average_targets(
        &mut self,
        brep: &BRep,
    ) -> Result<HashMap<BoundaryControlKey, Point3>, CoreError> {
        // Sorted: the per-class float sums accumulate key by
        // key, so map order would reach the averaged controls.
        let mut keys = self.parents.keys().copied().collect::<Vec<_>>();
        keys.sort_by_key(|key| format!("{key:?}"));
        let mut patches = HashMap::new();
        for key in &keys {
            patches
                .entry(key.face)
                .or_insert(rectangular_binonic_on_face(brep, key.face)?);
        }
        let mut sums: HashMap<BoundaryControlKey, (Point3, usize)> = HashMap::new();
        for key in &keys {
            let root = self.find(*key);
            let Some(patch) = patches.get(&key.face) else {
                continue;
            };
            let position = boundary_controls(patch, key.side)[key.index];
            sums.entry(root)
                .and_modify(|(sum, count)| {
                    *sum = *sum + (position - Point3::new(0.0, 0.0, 0.0));
                    *count += 1;
                })
                .or_insert((position, 1));
        }
        let mut targets = HashMap::new();
        for key in keys {
            let root = self.find(key);
            let Some(&(sum, count)) = sums.get(&root) else {
                continue;
            };
            let inverse = 1.0 / count as f64;
            targets.insert(
                key,
                Point3::new(sum.x * inverse, sum.y * inverse, sum.z * inverse),
            );
        }
        Ok(targets)
    }
}

fn push_trim_control_classes(trim: FullBoundaryTrim, classes: &mut SharedControlClasses) {
    for fraction_index in 0..10 {
        let key_a = normalize_corner_alias(BoundaryControlKey {
            face: trim.face_a,
            side: trim.side_a,
            index: natural_index_for_fraction(trim.start_a, trim.end_a, fraction_index),
        });
        let key_b = normalize_corner_alias(BoundaryControlKey {
            face: trim.face_b,
            side: trim.side_b,
            index: natural_index_for_fraction(trim.start_b, trim.end_b, 9 - fraction_index),
        });
        classes.union(key_a, key_b);
    }
}

/// One grid corner is addressed by two (side, index) pairs; both
/// rows share the control, so corner keys normalize to the U-row
/// representative before union (`boundary_controls` layout:
/// UMin = grid[0][j], UMax = grid[9][j], VMin = grid[i][0],
/// VMax = grid[i][9]).
fn normalize_corner_alias(key: BoundaryControlKey) -> BoundaryControlKey {
    let (side, index) = match (key.side, key.index) {
        (BoundarySide::VMin, 0) => (BoundarySide::UMin, 0),
        (BoundarySide::VMin, 9) => (BoundarySide::UMax, 0),
        (BoundarySide::VMax, 0) => (BoundarySide::UMin, 9),
        (BoundarySide::VMax, 9) => (BoundarySide::UMax, 9),
        other => other,
    };
    BoundaryControlKey {
        face: key.face,
        side,
        index,
    }
}

pub(super) fn fit_partial_shared_boundary_segments(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<&'static str, CoreError> {
    let segment_records =
        trim_ownership::segment_owned_boundary_trims(&full_boundary_trims(brep, slots)?);
    match segment_topology::commit_segment_boundary_topology(brep, &segment_records)? {
        segment_topology::SegmentBoundaryCommitStatus::AlreadySegmented => {
            fit_valid_segment_boundary_records(brep, segment_records)
        }
        segment_topology::SegmentBoundaryCommitStatus::SplitCommitted => {
            let refreshed_records =
                trim_ownership::segment_owned_boundary_trims(&full_boundary_trims(brep, slots)?);
            fit_valid_segment_boundary_records(brep, refreshed_records)
        }
        segment_topology::SegmentBoundaryCommitStatus::SplitDeclinedNoParameters => {
            Ok("partial_boundary_segment_topology_address_invalid")
        }
        segment_topology::SegmentBoundaryCommitStatus::SplitRolledBack => {
            Ok("partial_boundary_segment_edge_split_rolled_back")
        }
    }
}

fn fit_valid_segment_boundary_records(
    brep: &mut BRep,
    segment_records: Vec<trim_ownership::SegmentOwnedBoundaryTrim>,
) -> Result<&'static str, CoreError> {
    let mut samples = HashMap::new();
    let mut edited_intervals = HashMap::new();
    let mut exact_rows: HashMap<BoundaryRowKey, [Point3; 10]> = HashMap::new();
    let mut fallback_records = Vec::new();
    for record in &segment_records {
        if !push_exact_segment_row(brep, record.trim, &mut exact_rows)? {
            fallback_records.push(record.trim);
        }
    }
    apply_exact_rows(brep, &exact_rows)?;
    for trim in fallback_records {
        // Skip rows already authored exactly — the polynomial
        // restriction is authoritative; a sample fit on top of
        // it would only reintroduce approximation error.
        if exact_rows.contains_key(&BoundaryRowKey {
            face: trim.face_a,
            side: trim.side_a,
        }) || exact_rows.contains_key(&BoundaryRowKey {
            face: trim.face_b,
            side: trim.side_b,
        }) {
            continue;
        }
        push_partial_trim_samples(brep, trim, &mut samples, &mut edited_intervals)?;
    }
    if samples.is_empty() {
        if exact_rows.is_empty() {
            return Ok("partial_boundary_no_partial_trims");
        }
        return Ok("segment_topology_already_segmented_then_partial_boundary_segments_exact");
    }
    push_original_boundary_anchor_samples(brep, &mut samples, &edited_intervals)?;
    fit_boundary_row_samples(brep, &samples)?;
    if exact_rows.is_empty() {
        Ok("segment_topology_already_segmented_then_bounded_partial_boundary_segments_fit")
    } else {
        Ok("segment_topology_already_segmented_then_partial_boundary_segments_exact_and_fit")
    }
}

/// Exact T-junction ownership: when one trim side covers its
/// full unit row and the other a sub-interval, the sub-interval
/// side owns the geometry and the full-row side receives the
/// exact de Casteljau restriction of the owner's row over the
/// shared interval. Polynomial restriction is exact, so the trim
/// becomes a CP identity instead of a bounded least-squares fit.
/// Returns false (caller falls back to the sample fit) when both
/// sides are sub-intervals or the interval is degenerate.
fn push_exact_segment_row(
    brep: &BRep,
    trim: FullBoundaryTrim,
    exact_rows: &mut HashMap<BoundaryRowKey, [Point3; 10]>,
) -> Result<bool, CoreError> {
    let full_a = covers_full_row(trim.start_a, trim.end_a);
    let full_b = covers_full_row(trim.start_b, trim.end_b);
    let (owner_face, owner_side, owner_start, owner_end, receiver_face, receiver_side, reversed) =
        match (full_a, full_b) {
            (true, false) => (
                trim.face_b,
                trim.side_b,
                trim.start_b,
                trim.end_b,
                trim.face_a,
                trim.side_a,
                trim.start_a > trim.end_a,
            ),
            (false, true) => (
                trim.face_a,
                trim.side_a,
                trim.start_a,
                trim.end_a,
                trim.face_b,
                trim.side_b,
                trim.start_b > trim.end_b,
            ),
            _ => return Ok(false),
        };
    if (owner_end - owner_start).abs() < 1.0e-9 {
        return Ok(false);
    }
    let owner_patch = rectangular_binonic_on_face(brep, owner_face)?;
    let owner_controls = boundary_controls(&owner_patch, owner_side);
    // The trim walks the owner from owner_start to owner_end and
    // the receiver from receiver end to start (`9 - index`
    // correspondence, as in the full-span path); a reversed
    // receiver interval flips it back.
    let mut segment = bezier_row_segment(&owner_controls, owner_start, owner_end);
    if !reversed {
        segment.reverse();
    }
    exact_rows.insert(
        BoundaryRowKey {
            face: receiver_face,
            side: receiver_side,
        },
        segment,
    );
    Ok(true)
}

fn covers_full_row(start: f64, end: f64) -> bool {
    let lo = start.min(end);
    let hi = start.max(end);
    lo.abs() <= 1.0e-9 && (hi - 1.0).abs() <= 1.0e-9
}

fn apply_exact_rows(
    brep: &mut BRep,
    exact_rows: &HashMap<BoundaryRowKey, [Point3; 10]>,
) -> Result<(), CoreError> {
    let mut patches = HashMap::new();
    for key in exact_rows.keys() {
        patches
            .entry(key.face)
            .or_insert(rectangular_binonic_on_face(brep, key.face)?);
    }
    let mut exact_keys: Vec<_> = exact_rows.keys().collect();
    exact_keys.sort_by_key(|key| format!("{key:?}"));
    for key in exact_keys {
        let row = &exact_rows[key];
        let Some(patch) = patches.get_mut(&key.face) else {
            continue;
        };
        for (index, point) in row.iter().enumerate() {
            set_boundary_control(patch, key.side, index, *point);
        }
    }
    for (face, patch) in patches {
        set_rectangular_binonic_on_face(brep, face, patch)?;
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
    let mut sample_keys: Vec<_> = samples.keys().collect();
    sample_keys.sort_by_key(|key| format!("{key:?}"));
    for key in sample_keys {
        let row_samples = &samples[key];
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
    targets: &HashMap<BoundaryControlKey, Point3>,
) -> Result<(), CoreError> {
    let mut patches = HashMap::new();
    for key in targets.keys() {
        patches
            .entry(key.face)
            .or_insert(rectangular_binonic_on_face(brep, key.face)?);
    }
    let mut target_keys: Vec<_> = targets.keys().collect();
    target_keys.sort_by_key(|key| format!("{key:?}"));
    for key in target_keys {
        let target = &targets[key];
        let Some(patch) = patches.get_mut(&key.face) else {
            continue;
        };
        // Corner keys arrive normalized to the U-row alias; the
        // V-row alias addresses the same grid slot, so one write
        // covers both rows.
        set_boundary_control(patch, key.side, key.index, *target);
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
