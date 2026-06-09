//! Coupled worst-neighborhood fitting for KP20 support-cell boundaries.

use popkern_core::error::CoreError;
use popkern_core::topology::BRep;

use super::{
    FullBoundaryTrim, group_reconcile, rectangular_binonic_on_face, top_boundary_trims, trim_gap,
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

fn worst_trim_gap(brep: &BRep, trims: &[FullBoundaryTrim]) -> Result<f64, CoreError> {
    let mut worst = 0.0_f64;
    for &trim in trims {
        let patch_a = rectangular_binonic_on_face(brep, trim.face_a)?;
        let patch_b = rectangular_binonic_on_face(brep, trim.face_b)?;
        worst = worst.max(trim_gap(&patch_a, &patch_b, trim));
    }
    Ok(worst)
}
