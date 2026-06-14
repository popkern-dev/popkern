//! Target selection helpers for the measured worst-G1 normal solve.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::topology::FaceHandle;

use crate::construct::kp20_polar_support_cell_shared_trim_quality as shared_trim_quality;

use super::{BoundarySide, FullBoundaryTrim};

const LOCAL_NORMAL_SAMPLE_OFFSETS: [(f64, f64); 5] = [
    (0.0, 1.0),
    (-0.03125, 0.75),
    (0.03125, 0.75),
    (-0.0625, 0.375),
    (0.0625, 0.375),
];
const MAX_LOCAL_NORMAL_CONTROL_FRACTION: f64 = 0.04;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) struct ControlKey {
    pub(super) face: FaceHandle,
    pub(super) u: usize,
    pub(super) v: usize,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct ControlVariable {
    pub(super) key: ControlKey,
    pub(super) axis: usize,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct MatchedWorstTrim {
    pub(super) face_a: FaceHandle,
    pub(super) side_a: BoundarySide,
    pub(super) start_a: f64,
    pub(super) end_a: f64,
    pub(super) face_b: FaceHandle,
    pub(super) side_b: BoundarySide,
    pub(super) start_b: f64,
    pub(super) end_b: f64,
}

pub(super) fn control_variables(
    matched: MatchedWorstTrim,
    worst_fraction: f64,
) -> Vec<ControlVariable> {
    let mut controls = Vec::new();
    push_control_keys(
        &mut controls,
        matched.face_a,
        matched.side_a,
        matched.start_a,
        matched.end_a,
        worst_fraction,
    );
    push_control_keys(
        &mut controls,
        matched.face_b,
        matched.side_b,
        matched.start_b,
        matched.end_b,
        1.0 - worst_fraction,
    );
    controls
        .into_iter()
        .flat_map(|key| (0..3).map(move |axis| ControlVariable { key, axis }))
        .collect()
}

pub(super) fn match_worst_trim(
    trim: FullBoundaryTrim,
    worst: shared_trim_quality::Kp20SupportCellWorstTrim,
) -> Option<MatchedWorstTrim> {
    if matches_side(trim.face_a, trim.side_a, trim.face_b, trim.side_b, worst) {
        return Some(MatchedWorstTrim {
            face_a: trim.face_a,
            side_a: trim.side_a,
            start_a: trim.start_a,
            end_a: trim.end_a,
            face_b: trim.face_b,
            side_b: trim.side_b,
            start_b: trim.start_b,
            end_b: trim.end_b,
        });
    }
    if matches_side(trim.face_b, trim.side_b, trim.face_a, trim.side_a, worst) {
        return Some(MatchedWorstTrim {
            face_a: trim.face_b,
            side_a: trim.side_b,
            start_a: trim.start_b,
            end_a: trim.end_b,
            face_b: trim.face_a,
            side_b: trim.side_a,
            start_b: trim.start_a,
            end_b: trim.end_a,
        });
    }
    None
}

pub(super) fn local_sample_fractions(center: f64) -> Vec<(f64, f64)> {
    let mut samples = Vec::new();
    for (offset, weight) in LOCAL_NORMAL_SAMPLE_OFFSETS {
        let fraction = (center + offset).clamp(0.0, 1.0);
        if samples
            .iter()
            .all(|(existing, _): &(f64, f64)| (*existing - fraction).abs() > 1.0e-12)
        {
            samples.push((fraction, weight));
        }
    }
    samples
}

pub(super) fn local_control_limit(
    patch_a: &RectangularGregoryBiNonic,
    patch_b: &RectangularGregoryBiNonic,
    matched: MatchedWorstTrim,
) -> f64 {
    let length_a = side_point(patch_a, matched.side_a, matched.start_a).distance_to(&side_point(
        patch_a,
        matched.side_a,
        matched.end_a,
    ));
    let length_b = side_point(patch_b, matched.side_b, matched.start_b).distance_to(&side_point(
        patch_b,
        matched.side_b,
        matched.end_b,
    ));
    MAX_LOCAL_NORMAL_CONTROL_FRACTION * length_a.max(length_b).max(1.0e-9)
}

pub(super) fn uv_on_side(side: BoundarySide, start: f64, end: f64, fraction: f64) -> (f64, f64) {
    uv_from_parameter(side, start + (end - start) * fraction)
}

fn push_control_keys(
    controls: &mut Vec<ControlKey>,
    face: FaceHandle,
    side: BoundarySide,
    start: f64,
    end: f64,
    fraction: f64,
) {
    let parameter = start + (end - start) * fraction;
    let center = (parameter.clamp(0.0, 1.0) * 9.0).round() as isize;
    let mut inserted = false;
    for index in center - 2..=center + 2 {
        if (2..=7).contains(&index) {
            controls.push(adjacent_key(face, side, index as usize));
            inserted = true;
        }
    }
    if !inserted {
        controls.push(adjacent_key(face, side, center.clamp(2, 7) as usize));
    }
}

fn adjacent_key(face: FaceHandle, side: BoundarySide, index: usize) -> ControlKey {
    let (u, v) = match side {
        BoundarySide::UMin => (1, index),
        BoundarySide::UMax => (8, index),
        BoundarySide::VMin => (index, 1),
        BoundarySide::VMax => (index, 8),
    };
    ControlKey { face, u, v }
}

fn matches_side(
    face_a: FaceHandle,
    side_a: BoundarySide,
    face_b: FaceHandle,
    side_b: BoundarySide,
    worst: shared_trim_quality::Kp20SupportCellWorstTrim,
) -> bool {
    face_a == worst.face_a
        && row_for_side(side_a) == worst.row_a
        && face_b == worst.face_b
        && row_for_side(side_b) == worst.row_b
}

fn row_for_side(side: BoundarySide) -> Option<u8> {
    Some(match side {
        BoundarySide::UMin => 0,
        BoundarySide::UMax => 1,
        BoundarySide::VMin => 2,
        BoundarySide::VMax => 3,
    })
}

fn side_point(patch: &RectangularGregoryBiNonic, side: BoundarySide, parameter: f64) -> Point3 {
    let (u, v) = uv_from_parameter(side, parameter);
    patch.point_at(u, v)
}

fn uv_from_parameter(side: BoundarySide, parameter: f64) -> (f64, f64) {
    match side {
        BoundarySide::UMin => (0.0, parameter),
        BoundarySide::UMax => (1.0, parameter),
        BoundarySide::VMin => (parameter, 0.0),
        BoundarySide::VMax => (parameter, 1.0),
    }
}
