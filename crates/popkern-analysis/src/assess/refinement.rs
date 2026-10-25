//! Singularities by refinement growth: a stress singularity (a sharp
//! internal corner, a point-like support) shows as a peak that keeps
//! growing when the cells shrink, while a resolved concentration (a
//! fillet) settles. The caller solves the same problem at h and h/2;
//! this module only compares the two finished results.

use serde::Serialize;

use super::{AssessError, AssessOptions, Singularity, SolidPeaks, surface_peaks};
use crate::model::{AnalysisModel, FaceId};
use crate::structural::StructuralResult;

/// How one solid's von Mises peak changed from h to h/2, measured in
/// place: the largest fine value within `edge_cells` coarse cells of the
/// coarse peak (a global maximum can jump between symmetric features).
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct PeakGrowth {
    /// Local fine peak over the coarse peak.
    pub ratio: f64,
    /// Distance from the coarse peak to the local fine peak, in coarse
    /// cells.
    pub moved_cells: f64,
    /// A supported face within `edge_cells` coarse cells of the peak:
    /// the growth comes from the idealised support, not the geometry.
    pub support: Option<FaceId>,
}

/// A fine surface sample: position, von Mises, face (parent ids).
pub type FineSample = ([f64; 3], f64, FaceId);

/// The peak growth of each solid between `coarse` (cell h) and `fine`
/// (cell h/2) results of the same problem, whose supports are on
/// `supported` faces.
pub fn peak_growth(
    model: &AnalysisModel,
    coarse: &StructuralResult,
    fine: &StructuralResult,
    supported: &[FaceId],
    options: &AssessOptions,
) -> Result<Vec<PeakGrowth>, AssessError> {
    let samples: Vec<Vec<FineSample>> = fine
        .stresses
        .iter()
        .map(|s| {
            s.surface
                .iter()
                .map(|x| (x.position, x.stress.von_mises, x.face))
                .collect()
        })
        .collect();
    peak_growth_from(model, coarse, &samples, supported, options)
}

/// [`peak_growth`] from fine surface samples per solid (all of a
/// whole-part solve's, or a submodel's around the peak).
pub fn peak_growth_from(
    model: &AnalysisModel,
    coarse: &StructuralResult,
    fine: &[Vec<FineSample>],
    supported: &[FaceId],
    options: &AssessOptions,
) -> Result<Vec<PeakGrowth>, AssessError> {
    let h = coarse.report.cell_size.ok_or(AssessError::NoCellSize)?;
    if fine.len() != model.solids.len() {
        return Err(AssessError::Mismatch("solid count"));
    }
    let peaks = surface_peaks(model, coarse, options)?;
    let mut out = Vec::with_capacity(peaks.len());
    for (peak, samples) in peaks.iter().zip(fine) {
        let p = peak.von_mises.position;
        let mut best: Option<(f64, f64)> = None;
        let mut support: Option<(f64, FaceId)> = None;
        for &(q, vm, face) in samples {
            let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
            let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() / h;
            if dist > options.edge_cells {
                continue;
            }
            if best.is_none_or(|(v, _)| vm > v) {
                best = Some((vm, dist));
            }
            if supported.contains(&face) && support.is_none_or(|(e, _)| dist < e) {
                support = Some((dist, face));
            }
        }
        let (value, moved_cells) = best.unwrap_or((0.0, f64::INFINITY));
        let ratio = if peak.von_mises.value > 0.0 {
            value / peak.von_mises.value
        } else {
            1.0
        };
        out.push(PeakGrowth {
            ratio,
            moved_cells,
            support: support.map(|(_, f)| f),
        });
    }
    Ok(out)
}

/// Labels the peaks (from [`super::assess`] on either result) whose
/// growth exceeds `options.refinement_growth`: at a support's edge as
/// [`Singularity::SupportEdge`], elsewhere as
/// [`Singularity::RefinementGrowth`]. A peak already labelled at a
/// concave edge keeps that label.
pub fn label_by_growth(peaks: &mut [SolidPeaks], growth: &[PeakGrowth], options: &AssessOptions) {
    for (p, g) in peaks.iter_mut().zip(growth) {
        if p.singularity.is_some() || g.ratio <= options.refinement_growth {
            continue;
        }
        p.singularity = Some(match g.support {
            Some(face) => Singularity::SupportEdge {
                ratio: g.ratio,
                face,
            },
            None => Singularity::RefinementGrowth { ratio: g.ratio },
        });
    }
}

/// A local refinement check reaches this many coarse cells around each
/// peak: the comparison window (`edge_cells`, 2), the submodel's cap
/// clearance (1.5 at h/2) and room for the coarse cap data to fade.
pub const LOCAL_CHECK_CELLS: f64 = 6.0;

/// [`peak_growth`] with the h/2 solve local (PH-047): a submodel at h/2
/// around each solid's coarse peak (reaching [`LOCAL_CHECK_CELLS`]
/// coarse cells), driven by the coarse field, in place of the whole
/// part. Also returns the submodels' dofs in all.
pub fn local_peak_growth(
    model: &AnalysisModel,
    problem: &crate::structural::ElasticProblem,
    coarse: &StructuralResult,
    supported: &[FaceId],
    solve: &crate::options::AnalysisOptions,
    options: &AssessOptions,
) -> Result<(Vec<PeakGrowth>, usize), crate::error::AnalysisError> {
    let assess_err = |e: AssessError| crate::error::AnalysisError::Internal {
        stage: "assess",
        detail: e.to_string(),
    };
    let h = coarse
        .report
        .cell_size
        .ok_or(AssessError::NoCellSize)
        .map_err(assess_err)?;
    let peaks = surface_peaks(model, coarse, options).map_err(assess_err)?;
    let fine_options = crate::options::AnalysisOptions {
        target_cell_size: Some(0.5 * h),
        max_active_cells: crate::structural::submodel::SUBMODEL_CELL_BUDGET,
        grid_shift: [0.0; 3],
        ..*solve
    };
    let mut fine: Vec<Vec<FineSample>> = vec![Vec::new(); model.solids.len()];
    let mut dofs = 0;
    for (solid, peak) in peaks.iter().enumerate() {
        let p = peak.von_mises.position;
        let half = LOCAL_CHECK_CELLS * h;
        let sub = crate::structural::submodel::solve_submodel(
            model,
            problem,
            coarse,
            crate::model::SolidId(solid as u32),
            p.map(|c| c - half),
            p.map(|c| c + half),
            &fine_options,
            crate::structural::submodel::Caps::Mixed,
        )?;
        dofs += sub.result.report.dof;
        if let Some(slot) = fine.get_mut(solid) {
            *slot = sub.parent_samples(crate::structural::submodel::CAP_CLEARANCE_CELLS * 0.5 * h);
        }
    }
    let growth = peak_growth_from(model, coarse, &fine, supported, options).map_err(assess_err)?;
    Ok((growth, dofs))
}
