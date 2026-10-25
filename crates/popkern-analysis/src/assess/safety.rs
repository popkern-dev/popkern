//! Safety factor per the Stage 1 rule (notes 11), with the band of
//! cross-cutting Thread 6:
//! SF ∈ [σ_allow,low / (σ*(1+η)), σ_allow,mean / (σ*(1−η))],
//! where η is the local ZZ error at the governing peak and σ_allow,low
//! is the datasheet mean minus `scatter_sigmas` standard deviations.
//!
//! - In-plane: when the material yields, σ_inplane/σ_vM, or Raghava's
//!   paraboloidal criterion σ_vM² + (σ_c − σ_t) I₁ ≤ σ_c σ_t where a
//!   compressive strength is cited (PH-109); when brittle, Coulomb–Mohr
//!   σ₁⁺/σ_t − σ₃⁻/σ_c ≤ 1 where σ_c is cited (PH-084), else
//!   σ_inplane/σ₁⁺.
//! - Printed polymers (an inter-layer strength exists) add the
//!   inter-layer tension check σ_layer/σ₁⁺ and the inter-layer shear
//!   check S_layer/τ_max (PH-080) at the worst-case build orientation.
//! - With a build direction z_b (notes 11, orientation-aware rule), σ is
//!   read in build axes (`build_axes`): tension across the layers
//!   (z_b·σ·z_b)⁺, shear on them |σz_b − (z_b·σ·z_b)z_b|, and the
//!   in-plane check on the in-layer block.
//! - An infill core (PH-079) is checked against Bean 2025's gyroid
//!   yield strength (`structural::infill::core_yield`).
//! - Metals: σ_y/σ_vM.
//!
//! A mode without a sourced strength is listed as unchecked. Always an
//! estimate. The band excludes print defects (travel paths, cold joins;
//! Özen 2021) and slicer settings.

use serde::Serialize;

use crate::material::{Behaviour, Material};
use crate::model::AnalysisModel;
use crate::recover::core::CorePeak;
use crate::structural::StructuralResult;
use crate::structural::infill::{GYROID_YIELD_SCATTER, core_yield, gyroid_yield};

use super::build_axes::{Measures, measures, unit};
use super::{AssessError, AssessOptions};

/// The plain note shipped with every safety factor.
pub const ESTIMATE_NOTE: &str = "Estimated. The range covers the mesh error estimate and datasheet \
     scatter only; slicer settings and print quality (travel paths, cold joins between print \
     regions) can lower strength further.";

/// One strength check.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum Check {
    /// In-plane yield: σ_inplane / σ_vM.
    InPlaneYield,
    /// In-plane yield with a compressive yield strength: Raghava,
    /// Caddell & Yeh's paraboloid σ_vM² + (σ_c − σ_t) I₁ = σ_c σ_t,
    /// reported as σ_t over its equivalent stress
    /// [(σ_c − σ_t) I₁ + √((σ_c − σ_t)² I₁² + 4 σ_c σ_t σ_vM²)] / (2σ_c),
    /// which reads σ in uniaxial tension and (σ_t/σ_c)|σ| in compression.
    InPlaneParaboloidalYield,
    /// In-plane brittle: σ_inplane / σ₁⁺.
    InPlaneBrittle,
    /// In-plane brittle with a compressive strength: Coulomb–Mohr,
    /// σ_t / (σ₁⁺ − (σ_t/σ_c) σ₃⁻).
    InPlaneCoulombMohr,
    /// Between layers at the worst-case orientation: σ_layer / σ₁⁺
    /// ("limited by layer adhesion if printed with layers across this
    /// stress").
    InterLayer,
    /// Shear on the layers: S_layer / τ, τ the largest shear on any
    /// plane at the worst-case orientation.
    InterLayerShear,
    /// The infill core: the gyroid's yield strength over the core's
    /// equivalent stress (`recover::core::core_equivalent`).
    InfillCore,
}

/// What the band leaves out.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum BandGap {
    /// The datasheet gives no scatter for the allowable.
    ScatterUnsourced,
    /// The peak vertex has no honest error estimate (in-cell fallback
    /// or a thin wall).
    NoErrorEstimate,
    /// The local error is 100% or more: no upper bound.
    ErrorTooLarge,
}

/// One check's result.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CheckResult {
    /// The check.
    pub check: Check,
    /// Mean allowable, MPa.
    pub allowable: f64,
    /// Low allowable (mean − k·sd), MPa.
    pub allowable_low: f64,
    /// The governing stress over all solids, MPa: the check's
    /// equivalent stress (σ_vM, σ₁⁺, the Coulomb–Mohr σ₁⁺ − (σ_t/σ_c)σ₃⁻,
    /// a shear, or the core's equivalent).
    pub stress: f64,
    /// allowable / stress; `None` when the stress is not positive.
    pub factor: Option<f64>,
    /// Band [low, high]; high is `None` without an upper bound.
    pub band: Option<(f64, Option<f64>)>,
    /// Local relative error used for the band.
    pub error: Option<f64>,
    /// Where the stress peaks: (solid index, position).
    pub at: (usize, [f64; 3]),
    /// What the band leaves out.
    pub gaps: Vec<BandGap>,
}

/// A failure mode the assessment could not check.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum Unchecked {
    /// Shear between layers: the material has no sourced S_layer.
    InterLayerShear,
    /// An infill core no recovered sample reads: every core cell is
    /// mostly wall, so its stresses are the walls'.
    InfillCore,
}

/// The safety assessment.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SafetyAssessment {
    /// The unit build direction the layer checks used; `None` for the
    /// worst-case orientation (or a material without layers).
    pub build_direction: Option<[f64; 3]>,
    /// Modes left unchecked, with why.
    pub unchecked: Vec<Unchecked>,
    /// Every check that applies to the material.
    pub checks: Vec<CheckResult>,
    /// Index into `checks` of the lowest factor; `None` when no check
    /// sees a positive stress.
    pub governing: Option<usize>,
    /// Always true: a safety factor here is an estimate.
    pub estimated: bool,
    /// The plain note.
    pub note: &'static str,
}

/// A peak: (value, solid, position, local error, measures there).
type Peak = (f64, usize, [f64; 3], Option<f64>, Option<Measures>);

/// Largest value of `key` over every surface sample and, in an
/// infilled part, every interior wall sample (PH-113).
fn worst(
    result: &StructuralResult,
    direction: Option<[f64; 3]>,
    key: impl Fn(&Measures) -> f64,
) -> Peak {
    let mut best: Peak = (f64::NEG_INFINITY, 0, [0.0; 3], None, None);
    for (s, solid) in result.stresses.iter().enumerate() {
        let surface = solid.surface.iter().map(|x| (x.position, &x.stress));
        let interior = solid.interior.iter().map(|x| (x.position, &x.stress));
        for (position, stress) in surface.chain(interior) {
            let m = measures(&stress.stress, direction);
            let value = key(&m);
            if value > best.0 {
                best = (value, s, position, stress.error, Some(m));
            }
        }
    }
    best
}

/// A strength's mean, standard deviation (when sourced) and the low
/// value mean − k·sd.
#[derive(Debug, Clone, Copy)]
struct Allowable {
    mean: f64,
    sd: Option<f64>,
}

impl Allowable {
    fn low(&self, k: f64) -> f64 {
        (self.mean - k * self.sd.unwrap_or(0.0)).max(0.0)
    }
}

/// The check of `peak` against `allowable`; `low_stress` is the stress
/// the low allowable meets (differs from the peak's only when the
/// equivalent stress itself depends on the strengths, Coulomb–Mohr).
fn check(
    check: Check,
    allowable: Allowable,
    peak: &Peak,
    low_stress: f64,
    k: f64,
    sourced: bool,
) -> CheckResult {
    let (stress, s, v, error, _) = *peak;
    let mut gaps = Vec::new();
    if !sourced {
        gaps.push(BandGap::ScatterUnsourced);
    }
    let allowable_low = allowable.low(k);
    if error.is_none() {
        gaps.push(BandGap::NoErrorEstimate);
    }
    let eta = error.unwrap_or(0.0);
    let mean = allowable.mean;
    let (factor, band) = if stress > 0.0 {
        let high = if eta < 1.0 {
            Some(mean / (stress * (1.0 - eta)))
        } else {
            gaps.push(BandGap::ErrorTooLarge);
            None
        };
        (
            Some(mean / stress),
            Some((allowable_low / (low_stress * (1.0 + eta)), high)),
        )
    } else {
        (None, None)
    };
    CheckResult {
        check,
        allowable: mean,
        allowable_low,
        stress,
        factor,
        band,
        error,
        at: (s, v),
        gaps,
    }
}

/// A simple check: the equivalent stress does not depend on the
/// strength.
fn plain(check_kind: Check, allowable: Allowable, peak: &Peak, k: f64) -> CheckResult {
    check(
        check_kind,
        allowable,
        peak,
        peak.0,
        k,
        allowable.sd.is_some(),
    )
}

/// Coulomb–Mohr's equivalent stress σ₁⁺ − (σ_t/σ_c) σ₃⁻.
fn coulomb_mohr(m: &Measures, t: f64, c: f64) -> f64 {
    m.tension - (t / c) * m.compression
}

/// Raghava's paraboloidal equivalent stress (see
/// [`Check::InPlaneParaboloidalYield`]).
fn paraboloid(m: &Measures, t: f64, c: f64) -> f64 {
    let b = (c - t) * m.trace;
    (b + (b * b + 4.0 * c * t * m.von_mises * m.von_mises).sqrt()) / (2.0 * c)
}

/// A tension–compression check at the worst sample of `equivalent`
/// for the mean strengths; the band's low end re-reads that sample with
/// the low strengths.
fn two_strength(
    kind: Check,
    equivalent: fn(&Measures, f64, f64) -> f64,
    result: &StructuralResult,
    direction: Option<[f64; 3]>,
    (t, c): (Allowable, Allowable),
    k: f64,
) -> CheckResult {
    let peak = worst(result, direction, |m| equivalent(m, t.mean, c.mean));
    let (t_low, c_low) = (t.low(k), c.low(k));
    let low_stress = match peak.4 {
        Some(m) if c_low > 0.0 => equivalent(&m, t_low, c_low),
        _ => peak.0,
    };
    check(
        kind,
        t,
        &peak,
        low_stress,
        k,
        t.sd.is_some() && c.sd.is_some(),
    )
}

/// The worst core sample over every solid.
fn core_peak(result: &StructuralResult) -> Option<(usize, CorePeak)> {
    let mut best: Option<(usize, CorePeak)> = None;
    for (s, solid) in result.stresses.iter().enumerate() {
        if let Some(p) = solid.core_peak
            && best.is_none_or(|(_, b)| p.equivalent > b.equivalent)
        {
            best = Some((s, p));
        }
    }
    best
}

/// The safety assessment of `result` for `material`.
pub fn safety(
    model: &AnalysisModel,
    result: &StructuralResult,
    material: &Material,
    options: &AssessOptions,
) -> Result<SafetyAssessment, AssessError> {
    if result.stresses.len() != model.solids.len() {
        return Err(AssessError::Mismatch("solid count"));
    }
    let st = &material.strength;
    let inplane = st
        .inplane
        .value()
        .ok_or(AssessError::NoStrength(material.name))?;
    let allowable = |mean: f64, sd: &crate::material::Property| Allowable {
        mean,
        sd: sd.value(),
    };
    let inplane = allowable(inplane, &st.inplane_scatter);
    let k = options.scatter_sigmas;
    // A build direction only matters for a layered (printed) material.
    let direction = options
        .build_direction
        .filter(|_| st.interlayer.value().is_some())
        .and_then(unit);
    let compressive = st
        .ultimate_compressive
        .value()
        .map(|c| allowable(c, &st.ultimate_compressive_scatter));
    let mut checks = vec![match (st.inplane_behaviour, compressive) {
        (Behaviour::Yielding, Some(c)) => two_strength(
            Check::InPlaneParaboloidalYield,
            paraboloid,
            result,
            direction,
            (inplane, c),
            k,
        ),
        (Behaviour::Yielding, None) => plain(
            Check::InPlaneYield,
            inplane,
            &worst(result, direction, |m| m.von_mises),
            k,
        ),
        (Behaviour::Brittle, Some(c)) => two_strength(
            Check::InPlaneCoulombMohr,
            coulomb_mohr,
            result,
            direction,
            (inplane, c),
            k,
        ),
        (Behaviour::Brittle, None) => plain(
            Check::InPlaneBrittle,
            inplane,
            &worst(result, direction, |m| m.tension),
            k,
        ),
    }];
    let mut unchecked = Vec::new();
    if let Some(layer) = st.interlayer.value() {
        checks.push(plain(
            Check::InterLayer,
            allowable(layer, &st.interlayer_scatter),
            &worst(result, direction, |m| m.across),
            k,
        ));
        match st.interlayer_shear.value() {
            Some(shear) => checks.push(plain(
                Check::InterLayerShear,
                allowable(shear, &st.interlayer_shear_scatter),
                &worst(result, direction, |m| m.layer_shear),
                k,
            )),
            None => unchecked.push(Unchecked::InterLayerShear),
        }
    }
    if let Some(d) = result.report.infill_density {
        match core_peak(result) {
            Some((s, p)) => {
                // The gyroid's strength carries the parent's scatter
                // scaled by its ratio and the gyroid tests' own.
                let (compressive, shear) = gyroid_yield(d);
                let mean = core_yield(d, inplane.mean);
                let sd = inplane.sd.map(|sd| {
                    let parent = compressive.min(shear) * sd;
                    let gyroid = GYROID_YIELD_SCATTER * mean;
                    (parent * parent + gyroid * gyroid).sqrt()
                });
                let peak = (p.equivalent, s, p.position, p.error, None);
                checks.push(plain(Check::InfillCore, Allowable { mean, sd }, &peak, k));
            }
            None => unchecked.push(Unchecked::InfillCore),
        }
    }
    let governing = checks
        .iter()
        .enumerate()
        .filter_map(|(i, c)| c.factor.map(|f| (i, f)))
        .min_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)))
        .map(|(i, _)| i);
    Ok(SafetyAssessment {
        build_direction: direction,
        unchecked,
        checks,
        governing,
        estimated: true,
        note: ESTIMATE_NOTE,
    })
}
