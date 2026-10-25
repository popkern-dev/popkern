//! Steady thermal setup: fixed temperatures, heat flux, convection and
//! volumetric heat (map Stage 1b, notes 13).
//!
//! Convection coefficients are effective values (convection plus
//! linearised radiation). Presets take the low end of each range: a
//! lower h predicts a hotter part, the conservative side for the
//! service-temperature warning. The coefficient actually used is
//! always reported.

use crate::face_ref::Tolerance;
use crate::material::{MaterialError, Property, preset};
use crate::options::AnalysisOptions;
use crate::setup::{FaceDef, SetupError};

/// Most steps a transient may take.
pub const MAX_TRANSIENT_STEPS: usize = 100_000;

/// Stefan–Boltzmann constant, W/(m²·K⁴).
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;

/// Absolute zero in °C.
pub const ABSOLUTE_ZERO_C: f64 = -273.15;

/// A complete steady thermal setup.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalSetup {
    /// Material preset name.
    pub material: String,
    /// Default selector tolerances.
    pub tolerance: Tolerance,
    /// Named face references, sorted by name.
    pub faces: Vec<FaceDef>,
    /// Faces held at a temperature.
    pub fixed: Vec<FixedTemperature>,
    /// Heat flowing in through faces.
    pub flux: Vec<HeatFlux>,
    /// Convection to an ambient.
    pub convection: Vec<Convection>,
    /// Heat generated in the part.
    pub sources: Vec<VolumeSource>,
    /// A transient from a uniform start (`None`: steady).
    pub transient: Option<TransientSettings>,
    /// Discretisation and solver options.
    pub options: AnalysisOptions,
}

/// A transient run: from a uniform initial temperature, in fixed steps.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransientSettings {
    /// Initial temperature, °C.
    pub initial_c: f64,
    /// Time step, s.
    pub step_s: f64,
    /// Duration, s (rounded to whole steps).
    pub duration_s: f64,
}

impl TransientSettings {
    /// Number of steps: the duration in whole steps, at least one.
    pub fn steps(&self) -> usize {
        (self.duration_s / self.step_s).round().max(1.0) as usize
    }
}

/// Faces held at a fixed temperature.
#[derive(Debug, Clone, PartialEq)]
pub struct FixedTemperature {
    /// Names of face definitions.
    pub faces: Vec<String>,
    /// Temperature, °C.
    pub temperature_c: f64,
}

/// A total heat flow into the part, spread uniformly over faces.
#[derive(Debug, Clone, PartialEq)]
pub struct HeatFlux {
    /// Names of face definitions.
    pub faces: Vec<String>,
    /// Total power, W (negative draws heat out).
    pub power_w: f64,
}

/// Which faces a convection condition covers.
#[derive(Debug, Clone, PartialEq)]
pub enum ConvectionFaces {
    /// These named faces.
    Named(Vec<String>),
    /// Every face no other condition names.
    AllOther,
}

/// A heat transfer coefficient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HeatTransfer {
    /// A preset.
    Preset(ConvectionPreset),
    /// An effective coefficient, W/(m²·K).
    Value(f64),
}

/// Effective-coefficient presets (map Stage 1b).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConvectionPreset {
    /// Still air, matte or printed surface: 10 W/(m²·K).
    StillAirMatte,
    /// Still air, bright metal: 5 W/(m²·K).
    StillAirBrightMetal,
    /// Fan cooling, range 25–50: 25 W/(m²·K).
    Fan,
    /// Strong airflow, range 50–100: 50 W/(m²·K).
    StrongAirflow,
}

impl HeatTransfer {
    /// The coefficient in W/(m²·K).
    pub fn value(&self) -> f64 {
        match self {
            HeatTransfer::Value(h) => *h,
            HeatTransfer::Preset(ConvectionPreset::StillAirMatte) => 10.0,
            HeatTransfer::Preset(ConvectionPreset::StillAirBrightMetal) => 5.0,
            HeatTransfer::Preset(ConvectionPreset::Fan) => 25.0,
            HeatTransfer::Preset(ConvectionPreset::StrongAirflow) => 50.0,
        }
    }
}

/// Convection from faces to an ambient temperature.
#[derive(Debug, Clone, PartialEq)]
pub struct Convection {
    /// The faces.
    pub faces: ConvectionFaces,
    /// The coefficient: an effective one (convection plus linearised
    /// radiation) without `emittance`, the convective part alone with
    /// it.
    pub h: HeatTransfer,
    /// Ambient temperature, °C (also the radiating surroundings).
    pub ambient_c: f64,
    /// Grey radiation to the surroundings, solved exactly by fixed
    /// point (`None`: folded into an effective `h`).
    pub emittance: Option<Emittance>,
}

/// A surface emittance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Emittance {
    /// A value in (0, 1].
    Value(f64),
    /// The material's default.
    Material,
}

/// Heat generated uniformly in the whole part, or in one named body.
#[derive(Debug, Clone, PartialEq)]
pub struct VolumeSource {
    /// Total power, W.
    pub power_w: f64,
    /// The body (solid name) generating it; `None` for the whole part.
    pub body: Option<String>,
}

/// Linearised radiation coefficient h_rad = 4σT_m³ε, W/(m²·K), for mean
/// absolute temperature `t_mean_k` (Lienhard eqn 2.31).
pub fn radiation_coefficient(t_mean_k: f64, emittance: f64) -> f64 {
    4.0 * STEFAN_BOLTZMANN * t_mean_k * t_mean_k * t_mean_k * emittance
}

/// Why a thermal setup is invalid.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ThermalSetupError {
    /// A shared setup check failed (names, tolerances, options).
    #[error(transparent)]
    Setup(#[from] SetupError),
    /// No fixed temperature and no convection: the temperature level is
    /// undetermined (a constant null mode).
    #[error("no fixed temperature or convection gives the temperature a reference")]
    NoTemperatureReference,
    /// A heat source names a body the model does not have.
    #[error("heat source body `{name}` is not in the model (bodies: {known:?})")]
    UnknownBody {
        /// The name given.
        name: String,
        /// The model's body names.
        known: Vec<String>,
    },
    /// A face name carries more than one condition.
    #[error("face `{0}` carries more than one thermal condition")]
    ConflictingBoundary(String),
    /// A temperature is non-finite or below absolute zero.
    #[error("{what} temperature {value} °C is not valid")]
    BadTemperature {
        /// Which condition.
        what: String,
        /// The value.
        value: f64,
    },
    /// A power is non-finite or zero.
    #[error("{what} power {value} W is not valid")]
    BadPower {
        /// Which condition.
        what: String,
        /// The value.
        value: f64,
    },
    /// A transient's step or duration is not positive and finite, or it
    /// would take more than [`MAX_TRANSIENT_STEPS`] steps.
    #[error("transient: step {step_s} s over {duration_s} s is not valid")]
    BadTransient {
        /// Step, s.
        step_s: f64,
        /// Duration, s.
        duration_s: f64,
    },
    /// An emittance is outside (0, 1].
    #[error("convection {index}: emittance {value} is not in (0, 1]")]
    BadEmittance {
        /// Convection entry.
        index: usize,
        /// The value.
        value: f64,
    },
    /// A heat transfer coefficient is not positive and finite.
    #[error("convection {index}: h = {value} W/m²K is not positive")]
    BadCoefficient {
        /// Convection entry.
        index: usize,
        /// The value.
        value: f64,
    },
    /// The material cannot be used for a thermal solve.
    #[error(transparent)]
    Material(#[from] MaterialError),
}

fn temperature(what: String, value: f64) -> Result<(), ThermalSetupError> {
    if value.is_finite() && value > ABSOLUTE_ZERO_C {
        Ok(())
    } else {
        Err(ThermalSetupError::BadTemperature { what, value })
    }
}

impl ThermalSetup {
    /// Checks the setup before any model is seen; returns the material's
    /// conductivity in W/(m·K).
    pub fn validate(&self) -> Result<f64, ThermalSetupError> {
        let material = preset(&self.material)?;
        let conductivity = match material.conductivity {
            Property::Cited { value, .. } | Property::Typical { value, .. } => value,
            _ => {
                return Err(MaterialError::MissingProperty {
                    material: material.name,
                    property: "conductivity",
                }
                .into());
            }
        };
        let names = |list: &[String], used_by: String| -> Result<(), ThermalSetupError> {
            if list.is_empty() {
                return Err(SetupError::NoFaces(used_by).into());
            }
            match list
                .iter()
                .find(|n| !self.faces.iter().any(|d| &d.name == *n))
            {
                Some(name) => Err(SetupError::UndefinedFace {
                    used_by,
                    name: name.clone(),
                }
                .into()),
                None => Ok(()),
            }
        };
        let mut used: Vec<String> = Vec::new();
        let mut claim = |list: &[String]| -> Result<(), ThermalSetupError> {
            for n in list {
                if used.contains(n) {
                    return Err(ThermalSetupError::ConflictingBoundary(n.clone()));
                }
                used.push(n.clone());
            }
            Ok(())
        };
        for (i, f) in self.fixed.iter().enumerate() {
            names(&f.faces, format!("fixed temperature {i}"))?;
            temperature(format!("fixed temperature {i}"), f.temperature_c)?;
            claim(&f.faces)?;
        }
        for (i, f) in self.flux.iter().enumerate() {
            names(&f.faces, format!("heat flux {i}"))?;
            if !(f.power_w.is_finite() && f.power_w != 0.0) {
                return Err(ThermalSetupError::BadPower {
                    what: format!("heat flux {i}"),
                    value: f.power_w,
                });
            }
            claim(&f.faces)?;
        }
        let mut all_other = 0;
        for (i, c) in self.convection.iter().enumerate() {
            match &c.faces {
                ConvectionFaces::Named(list) => {
                    names(list, format!("convection {i}"))?;
                    claim(list)?;
                }
                ConvectionFaces::AllOther => all_other += 1,
            }
            if let Some(Emittance::Value(e)) = c.emittance
                && !(e.is_finite() && e > 0.0 && e <= 1.0)
            {
                return Err(ThermalSetupError::BadEmittance { index: i, value: e });
            }
            let h = c.h.value();
            if !(h.is_finite() && h > 0.0) {
                return Err(ThermalSetupError::BadCoefficient { index: i, value: h });
            }
            temperature(format!("convection {i} ambient"), c.ambient_c)?;
        }
        if all_other > 1 {
            return Err(ThermalSetupError::ConflictingBoundary(
                "all other faces".to_string(),
            ));
        }
        for (i, s) in self.sources.iter().enumerate() {
            if !(s.power_w.is_finite() && s.power_w != 0.0) {
                return Err(ThermalSetupError::BadPower {
                    what: format!("heat source {i}"),
                    value: s.power_w,
                });
            }
        }
        if let Some(t) = self.transient {
            temperature("transient initial".to_string(), t.initial_c)?;
            let positive = |x: f64| x.is_finite() && x > 0.0;
            if !(positive(t.step_s) && positive(t.duration_s))
                || t.duration_s / t.step_s > MAX_TRANSIENT_STEPS as f64
            {
                return Err(ThermalSetupError::BadTransient {
                    step_s: t.step_s,
                    duration_s: t.duration_s,
                });
            }
        } else if self.fixed.is_empty() && self.convection.is_empty() {
            // A transient's capacity fixes the level; a steady solve
            // needs a reference.
            return Err(ThermalSetupError::NoTemperatureReference);
        }
        // Shared checks: tolerances, face definitions, options.
        let shared = crate::setup::StructuralSetup {
            material: self.material.clone(),
            build_direction: None,
            tolerance: self.tolerance,
            faces: self.faces.clone(),
            supports: Vec::new(),
            loads: Vec::new(),
            thermal: None,
            modal: None,
            infill: None,
            generative: None,
            options: self.options,
        };
        shared.validate()?;
        Ok(conductivity)
    }
}

#[cfg(test)]
#[path = "thermal_setup_tests.rs"]
mod thermal_setup_tests;
