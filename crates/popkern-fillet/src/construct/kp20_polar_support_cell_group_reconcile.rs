//! Transactional grouped row fitting for KP20 support-cell boundaries.

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::topology::{BRep, FaceHandle};
use std::collections::HashMap;

use super::{
    BoundarySide, FullBoundaryTrim, average_point, lerp_scalar, rectangular_binonic_on_face,
    row_fit, set_rectangular_binonic_on_face, trim_gap,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct BoundaryRowKey {
    face: FaceHandle,
    side: BoundarySide,
}

pub(super) fn reconcile_grouped_boundary_rows(
    brep: &mut BRep,
    trims: &[FullBoundaryTrim],
) -> Result<(), CoreError> {
    if trims.len() < 2 {
        return Ok(());
    }
    let before_gap = worst_trim_gap_on_brep(brep, trims)?;
    let originals = trim_patches(brep, trims)?;
    let mut candidates = originals.clone();
    let samples = grouped_boundary_samples(&originals, trims)?;
    for (key, row_samples) in &samples {
        let Some((_, patch)) = candidates.iter_mut().find(|(face, _)| *face == key.face) else {
            continue;
        };
        row_fit::fit_boundary_samples(patch, key.side, row_samples);
    }
    if !samples.iter().all(|(key, _)| {
        patch_for(&originals, key.face)
            .zip(patch_for(&candidates, key.face))
            .is_some_and(|(before, after)| {
                row_fit::max_boundary_displacement(before, after, key.side) <= 0.1
            })
    }) {
        return Ok(());
    }
    let after_gap = worst_trim_gap_on_patches(&candidates, trims)?;
    if after_gap + 1.0e-12 >= before_gap {
        return Ok(());
    }
    set_trim_patches(brep, &candidates)?;
    if brep.validate().is_err() {
        set_trim_patches(brep, &originals)?;
    }
    Ok(())
}

fn grouped_boundary_samples(
    patches: &[(FaceHandle, RectangularGregoryBiNonic)],
    trims: &[FullBoundaryTrim],
) -> Result<HashMap<BoundaryRowKey, Vec<row_fit::BoundaryRowFitSample>>, CoreError> {
    let mut out: HashMap<BoundaryRowKey, Vec<row_fit::BoundaryRowFitSample>> = HashMap::new();
    for trim in trims {
        let patch_a = patch_for(patches, trim.face_a).ok_or_else(|| {
            CoreError::PreconditionFailed("kp20 support-cell grouped row missing face A".into())
        })?;
        let patch_b = patch_for(patches, trim.face_b).ok_or_else(|| {
            CoreError::PreconditionFailed("kp20 support-cell grouped row missing face B".into())
        })?;
        for index in 0..10 {
            let fraction = index as f64 / 9.0;
            let param_a = lerp_scalar(trim.start_a, trim.end_a, fraction);
            let param_b = lerp_scalar(trim.start_b, trim.end_b, 1.0 - fraction);
            let target = average_point(
                row_fit::eval_boundary(patch_a, trim.side_a, param_a),
                row_fit::eval_boundary(patch_b, trim.side_b, param_b),
            );
            push_boundary_sample(&mut out, trim.face_a, trim.side_a, param_a, target);
            push_boundary_sample(&mut out, trim.face_b, trim.side_b, param_b, target);
        }
    }
    Ok(out)
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

fn trim_patches(
    brep: &BRep,
    trims: &[FullBoundaryTrim],
) -> Result<Vec<(FaceHandle, RectangularGregoryBiNonic)>, CoreError> {
    let mut out = Vec::new();
    for face in trims.iter().flat_map(|trim| [trim.face_a, trim.face_b]) {
        if out.iter().any(|(existing, _)| *existing == face) {
            continue;
        }
        out.push((face, rectangular_binonic_on_face(brep, face)?));
    }
    Ok(out)
}

fn set_trim_patches(
    brep: &mut BRep,
    patches: &[(FaceHandle, RectangularGregoryBiNonic)],
) -> Result<(), CoreError> {
    for (face, patch) in patches {
        set_rectangular_binonic_on_face(brep, *face, patch.clone())?;
    }
    Ok(())
}

fn patch_for(
    patches: &[(FaceHandle, RectangularGregoryBiNonic)],
    face: FaceHandle,
) -> Option<&RectangularGregoryBiNonic> {
    patches
        .iter()
        .find(|(candidate, _)| *candidate == face)
        .map(|(_, patch)| patch)
}

fn worst_trim_gap_on_brep(brep: &BRep, trims: &[FullBoundaryTrim]) -> Result<f64, CoreError> {
    let patches = trim_patches(brep, trims)?;
    worst_trim_gap_on_patches(&patches, trims)
}

fn worst_trim_gap_on_patches(
    patches: &[(FaceHandle, RectangularGregoryBiNonic)],
    trims: &[FullBoundaryTrim],
) -> Result<f64, CoreError> {
    let mut worst = 0.0_f64;
    for &trim in trims {
        let patch_a = patch_for(patches, trim.face_a).ok_or_else(|| {
            CoreError::PreconditionFailed("kp20 support-cell trim gap missing face A".into())
        })?;
        let patch_b = patch_for(patches, trim.face_b).ok_or_else(|| {
            CoreError::PreconditionFailed("kp20 support-cell trim gap missing face B".into())
        })?;
        worst = worst.max(trim_gap(patch_a, patch_b, trim));
    }
    Ok(worst)
}
