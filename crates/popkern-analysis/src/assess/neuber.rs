//! Local yielding (map Stage 2, item 7): where a ductile part's elastic
//! peak passes yield over a small region, Neuber's rule σε = σ_e²/E
//! (σ_e the linear-elastic peak) estimates the local stress and strain
//! the surrounding elastic material allows. No preset has a sourced
//! stress–strain curve, so the material is taken elastic–perfectly
//! plastic: σ = σ_y and ε = σ_e²/(E σ_y). Where yield spreads over more
//! than [`LOCAL_FRACTION`] of the surface the yielding is not local, the
//! linear solution does not describe the part, and the result says so
//! instead of a Neuber value. Brittle materials never get this: a peak
//! above their strength is fracture, which the safety factor reports.

use serde::Serialize;

use crate::material::{Behaviour, Material};
use crate::structural::StructuralResult;

/// Yielded surface area above this fraction of the part's surface is
/// gross yielding, not local.
pub const LOCAL_FRACTION: f64 = 0.01;

/// The note shipped with every local-yielding estimate.
pub const NEUBER_NOTE: &str = "Estimated by Neuber's rule with an elastic–perfectly-plastic \
     material (no sourced stress–strain curve); valid while the yielded zone stays small.";

/// A ductile peak above yield.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LocalYielding {
    /// The linear-elastic peak von Mises, MPa.
    pub elastic_peak_mpa: f64,
    /// Yield strength, MPa.
    pub yield_mpa: f64,
    /// Where: (solid index, position).
    pub at: (usize, [f64; 3]),
    /// Surface area above yield over the part's surface.
    pub yielded_fraction: f64,
    /// Whether the yielding is local (`yielded_fraction` ≤
    /// [`LOCAL_FRACTION`]); when not, the estimates below are absent.
    pub local: bool,
    /// Neuber's local stress, MPa (σ_y).
    pub local_stress_mpa: Option<f64>,
    /// Neuber's local strain σ_e²/(E σ_y).
    pub local_strain: Option<f64>,
    /// Its plastic part, local strain − σ_y/E.
    pub plastic_strain: Option<f64>,
    /// Always true.
    pub estimated: bool,
    /// The plain note.
    pub note: &'static str,
}

/// The local yielding of `result` for `material` (Young's modulus
/// `youngs`, MPa), when its in-plane check yields and the peak passes
/// yield.
pub fn local_yielding(
    result: &StructuralResult,
    material: &Material,
    youngs: f64,
) -> Option<LocalYielding> {
    if material.strength.inplane_behaviour != Behaviour::Yielding {
        return None;
    }
    let sigma_y = material.strength.inplane.value()?;
    let (mut peak, mut at) = (f64::NEG_INFINITY, (0, [0.0; 3]));
    let (mut yielded, mut total) = (0.0, 0.0);
    for (s, solid) in result.stresses.iter().enumerate() {
        for x in &solid.surface {
            total += x.area;
            if x.stress.von_mises > sigma_y {
                yielded += x.area;
            }
            if x.stress.von_mises > peak {
                peak = x.stress.von_mises;
                at = (s, x.position);
            }
        }
    }
    if peak.is_nan() || peak <= sigma_y || youngs.is_nan() || youngs <= 0.0 {
        return None;
    }
    let yielded_fraction = if total > 0.0 { yielded / total } else { 0.0 };
    let local = yielded_fraction <= LOCAL_FRACTION;
    let strain = peak * peak / (youngs * sigma_y);
    Some(LocalYielding {
        elastic_peak_mpa: peak,
        yield_mpa: sigma_y,
        at,
        yielded_fraction,
        local,
        local_stress_mpa: local.then_some(sigma_y),
        local_strain: local.then_some(strain),
        plastic_strain: local.then_some(strain - sigma_y / youngs),
        estimated: true,
        note: NEUBER_NOTE,
    })
}

#[cfg(test)]
#[path = "neuber_tests.rs"]
mod tests;
