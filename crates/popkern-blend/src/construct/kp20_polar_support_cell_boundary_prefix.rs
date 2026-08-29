//! Bounded ranked-prefix fitting for KP20 support-cell boundaries.

use popkern_core::error::CoreError;
use popkern_core::topology::BRep;

use super::{FullBoundaryTrim, group_reconcile, rectangular_binonic_on_face, trim_gap};

pub(super) fn reconcile_ranked_boundary_prefix(
    brep: &mut BRep,
    trims: &[FullBoundaryTrim],
) -> Result<(), CoreError> {
    let ranked = ranked_boundary_trims(brep, trims)?;
    let before_gap = worst_trim_gap(brep, trims)?;
    let mut best_brep = None;
    let mut best_gap = before_gap;
    for count in prefix_counts(ranked.len()) {
        let prefix = ranked.iter().copied().take(count).collect::<Vec<_>>();
        let mut candidate = brep.clone();
        for trim in &prefix {
            super::reconcile_full_boundary_trim(&mut candidate, *trim)?;
        }
        group_reconcile::reconcile_grouped_boundary_rows(&mut candidate, &prefix)?;
        if candidate.validate().is_err() {
            continue;
        }
        let after_gap = worst_trim_gap(&candidate, trims)?;
        if after_gap + 1.0e-12 < best_gap {
            best_gap = after_gap;
            best_brep = Some(candidate);
        }
    }
    if let Some(candidate) = best_brep {
        *brep = candidate;
    }
    Ok(())
}

fn ranked_boundary_trims(
    brep: &BRep,
    trims: &[FullBoundaryTrim],
) -> Result<Vec<FullBoundaryTrim>, CoreError> {
    let mut ranked = Vec::with_capacity(trims.len());
    for &trim in trims {
        let patch_a = rectangular_binonic_on_face(brep, trim.face_a)?;
        let patch_b = rectangular_binonic_on_face(brep, trim.face_b)?;
        ranked.push((trim, trim_gap(&patch_a, &patch_b, trim)));
    }
    ranked.sort_by(|(_, a), (_, b)| b.total_cmp(a));
    Ok(ranked.into_iter().map(|(trim, _)| trim).collect())
}

fn prefix_counts(trim_count: usize) -> Vec<usize> {
    let mut counts = [3, 6, 12, trim_count / 2, trim_count]
        .into_iter()
        .filter(|count| *count > 0 && *count <= trim_count)
        .collect::<Vec<_>>();
    counts.sort_unstable();
    counts.dedup();
    counts
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
