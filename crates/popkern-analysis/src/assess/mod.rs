//! What a structural result means for the part (map §1.6 steps 5–6):
//! surface peaks with their location, the singularity label at sharp
//! internal corners, the fillet-resolution check, and the safety factor
//! with its band. Everything here reads a finished
//! [`StructuralResult`]; nothing is re-solved.

mod build_axes;
pub mod fillet;
pub mod neuber;
pub mod refinement;
pub mod safety;

use serde::Serialize;

use crate::model::{AnalysisModel, FaceId, SolidId};
use crate::recover::VertexStress;
use crate::recover::sample::SurfaceSample;
use crate::structural::StructuralResult;

pub use fillet::{FilletCheck, bore_checks, concave_raisers, fillet_checks};
pub use refinement::{PeakGrowth, label_by_growth, local_peak_growth, peak_growth};
pub use safety::{Check, CheckResult, SafetyAssessment, safety};

/// Thresholds and policies of the assessment; [`AssessOptions::default`]
/// holds the measured Stage 1 values.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct AssessOptions {
    /// Fillets with r/h below this many cells warn that their peak is
    /// under-resolved. Measured (PH-027): the C8/C9 fillet and C10
    /// bending peaks read within 10% at 4 cells over 4 grid shifts
    /// (worst 0.933), and up to 14% low at 3.
    pub fillet_resolved_cells: f64,
    /// Bores, holes and slot ends (concave cylinders without a role)
    /// with r/h below this many cells warn that their peak is
    /// under-resolved. Measured: C10 tension reads within 10% only from
    /// 6 cells (0.92–0.95; PH-027), and bored_block's bore in tension
    /// 0.875 of CalculiX at 4 cells (PH-050).
    pub bore_resolved_cells: f64,
    /// Fillets (and bores) with r/h below this are treated as sharp
    /// corners
    /// (cross-cutting Thread 5).
    pub fillet_sharp_cells: f64,
    /// A face is in a stressed region when its peak von Mises is at
    /// least this fraction of the solid's.
    pub stressed_fraction: f64,
    /// Quantile of the area-weighted surface von Mises reported next to
    /// the peak (information only; see [`Singularity`]).
    pub quantile: f64,
    /// A peak within this many cells of a concave sharp edge is at that
    /// edge. Measured (PH-027): on the L-bracket at 2 mm cells the
    /// recovered peak sits 1.5 cells from the corner line.
    pub edge_cells: f64,
    /// The band's low allowable is the mean minus this many datasheet
    /// standard deviations.
    pub scatter_sigmas: f64,
    /// Build direction for the layer checks; `None` checks the
    /// worst-case orientation.
    pub build_direction: Option<[f64; 3]>,
    /// A peak growing by more than this factor from h to h/2 is a
    /// singularity ([`refinement`]). Measured (PH-038): the L-bracket's
    /// 90° re-entrant corner grows 1.196 (2 → 1 mm) and 1.361 (1 → 0.5
    /// mm), against about 1.37 from its r^(−0.456) field; a resolved
    /// shoulder fillet 1.08, 1.02 and 1.07 at 3, 4 and 5 cells.
    pub refinement_growth: f64,
}

impl Default for AssessOptions {
    fn default() -> Self {
        AssessOptions {
            fillet_resolved_cells: 4.0,
            bore_resolved_cells: 6.0,
            fillet_sharp_cells: 0.3,
            stressed_fraction: 0.5,
            quantile: 0.995,
            edge_cells: 2.0,
            scatter_sigmas: 2.0,
            refinement_growth: 1.15,
            build_direction: None,
        }
    }
}

/// One surface peak.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Peak {
    /// The solid.
    pub solid: SolidId,
    /// Position, mm.
    pub position: [f64; 3],
    /// The face it lies on.
    pub face: FaceId,
    /// The peak value, MPa.
    pub value: f64,
    /// Local relative ZZ error there; `None` where no estimate is
    /// honest.
    pub error: Option<f64>,
}

/// Why a peak is labelled a singularity: "stress concentration at a
/// sharp internal corner; add a fillet" (notes 04).
///
/// The statistical test of map §1.6 (max / P99.5 > 1 + 2·local error)
/// was measured on one solve and does not separate the cases (PH-027):
/// a resolved C8 fillet scores 1.47–1.58 at 4–5 cells, the sharp
/// L-bracket corner 1.20–1.50 at 1 and 0.5 mm cells. Two tests are
/// used instead: the BRep's concave sharp edges (one solve), and growth
/// under refinement (two solves, [`refinement`]).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub enum Singularity {
    /// The peak lies on (within `edge_cells` of) a concave sharp edge.
    ConcaveEdge {
        /// The faces meeting at the edge.
        faces: (FaceId, FaceId),
    },
    /// The peak grew by `ratio` from h to h/2 in place: a sharp corner
    /// the BRep does not report ("add a fillet").
    RefinementGrowth {
        /// Fine peak over coarse peak.
        ratio: f64,
    },
    /// The peak grew by `ratio` at the edge of a support on `face`: the
    /// idealised support concentrates stress there, so read stresses
    /// away from it or spread the support; a fillet will not help.
    SupportEdge {
        /// Fine peak over coarse peak.
        ratio: f64,
        /// The supported face.
        face: FaceId,
    },
}

/// The surface peaks of one solid.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SolidPeaks {
    /// Largest von Mises.
    pub von_mises: Peak,
    /// Largest principal stress.
    pub max_principal: Peak,
    /// Area-weighted `quantile` of von Mises over the surface, MPa.
    pub von_mises_quantile: f64,
    /// Peak von Mises over that quantile (information only).
    pub peak_ratio: Option<f64>,
    /// Set when the von Mises peak is a singularity.
    pub singularity: Option<Singularity>,
}

/// The assessment of a structural result.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Assessment {
    /// Options used.
    pub options: AssessOptions,
    /// Peaks per solid.
    pub peaks: Vec<SolidPeaks>,
    /// Concave fillets in stressed regions, with their resolution.
    pub fillets: Vec<FilletCheck>,
    /// Bores, holes and slot ends in stressed regions, with their
    /// resolution.
    pub bores: Vec<FilletCheck>,
    /// Safety factor; always an estimate.
    pub safety: SafetyAssessment,
    /// A ductile peak above yield: Neuber's local estimate, or gross
    /// yielding (`neuber`).
    pub local_yielding: Option<neuber::LocalYielding>,
}

/// Why an assessment could not be made.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum AssessError {
    /// The result does not belong to the model (solid or vertex counts
    /// differ).
    #[error("result does not match the model: {0}")]
    Mismatch(&'static str),
    /// The result has no cell size.
    #[error("result has no cell size")]
    NoCellSize,
    /// A solid has no vertices.
    #[error("solid {0:?} has no vertices")]
    EmptySolid(SolidId),
    /// The material has no in-plane strength to check against.
    #[error("material `{0}` has no in-plane strength")]
    NoStrength(&'static str),
}

/// The area-weighted `q` quantile of `values`.
fn weighted_quantile(values: &[f64], weights: &[f64], q: f64) -> f64 {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&a, &b| values[a].total_cmp(&values[b]).then(a.cmp(&b)));
    let total: f64 = weights.iter().sum();
    let mut acc = 0.0;
    for &i in &order {
        acc += weights[i];
        if acc >= q * total {
            return values[i];
        }
    }
    order.last().map_or(0.0, |&i| values[i])
}

/// The sample maximising `key` (lowest index on ties).
fn argmax(samples: &[SurfaceSample], key: impl Fn(&VertexStress) -> f64) -> usize {
    let mut best = 0;
    for (i, s) in samples.iter().enumerate() {
        if key(&s.stress) > key(&samples[best].stress) {
            best = i;
        }
    }
    best
}

fn distance_to_polyline(p: [f64; 3], line: &[[f64; 3]]) -> f64 {
    let mut best = f64::INFINITY;
    for w in line.windows(2) {
        let (a, b) = (w[0], w[1]);
        let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let len2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        let t = if len2 > 0.0 {
            (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1] + (p[2] - a[2]) * d[2]) / len2)
                .clamp(0.0, 1.0)
        } else {
            0.0
        };
        let q = [a[0] + t * d[0], a[1] + t * d[1], a[2] + t * d[2]];
        let e = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
        best = best.min((e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt());
    }
    best
}

/// Surface peaks of every solid, with the singularity label.
pub fn surface_peaks(
    model: &AnalysisModel,
    result: &StructuralResult,
    options: &AssessOptions,
) -> Result<Vec<SolidPeaks>, AssessError> {
    let h = result.report.cell_size.ok_or(AssessError::NoCellSize)?;
    if result.stresses.len() != model.solids.len() {
        return Err(AssessError::Mismatch("solid count"));
    }
    let mut out = Vec::new();
    for (mesh, stresses) in model.solids.iter().zip(&result.stresses) {
        let samples = &stresses.surface;
        if samples.is_empty() {
            return Err(AssessError::EmptySolid(mesh.id));
        }
        let peak = |i: usize, value: f64| Peak {
            solid: mesh.id,
            position: samples[i].position,
            face: samples[i].face,
            value,
            error: samples[i].stress.error,
        };
        let iv = argmax(samples, |s| s.von_mises);
        let ip = argmax(samples, |s| s.max_principal);
        let von_mises = peak(iv, samples[iv].stress.von_mises);
        let vm: Vec<f64> = samples.iter().map(|s| s.stress.von_mises).collect();
        let area: Vec<f64> = samples.iter().map(|s| s.area).collect();
        let quantile = weighted_quantile(&vm, &area, options.quantile);
        let own = |f: FaceId| model.faces.get(f.0 as usize).map(|i| i.solid) == Some(mesh.id);
        let at_edge = model.concave_edges.iter().find(|e| {
            own(e.faces.0)
                && distance_to_polyline(von_mises.position, &e.polyline) <= options.edge_cells * h
        });
        let singularity = at_edge.map(|e| Singularity::ConcaveEdge { faces: e.faces });
        let peak_ratio = (quantile > 0.0).then(|| von_mises.value / quantile);
        out.push(SolidPeaks {
            max_principal: peak(ip, samples[ip].stress.max_principal),
            von_mises,
            von_mises_quantile: quantile,
            peak_ratio,
            singularity,
        });
    }
    Ok(out)
}

/// The full assessment of `result` for `material`.
pub fn assess(
    model: &AnalysisModel,
    result: &StructuralResult,
    material: &crate::material::Material,
    options: &AssessOptions,
) -> Result<Assessment, AssessError> {
    let peaks = surface_peaks(model, result, options)?;
    let fillets = fillet_checks(model, result, options)?;
    let bores = bore_checks(model, result, options)?;
    let safety = safety(model, result, material, options)?;
    let local_yielding = material
        .youngs_modulus
        .value()
        .and_then(|e| neuber::local_yielding(result, material, e));
    Ok(Assessment {
        options: *options,
        peaks,
        fillets,
        bores,
        safety,
        local_yielding,
    })
}

#[cfg(test)]
#[path = "assess_tests.rs"]
mod assess_tests;

#[cfg(test)]
#[path = "strength_tests.rs"]
mod strength_tests;
