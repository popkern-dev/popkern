//! The thermal setup file (TOML): `analysis = "thermal"` plus faces,
//! fixed temperatures, heat flux, convection and heat sources. Every
//! table rejects unknown keys; a parsed file is a validated setup.

use std::collections::BTreeMap;

use popkern_analysis::thermal_setup::{
    Convection, ConvectionFaces, ConvectionPreset, Emittance, FixedTemperature, HeatFlux,
    HeatTransfer, ThermalSetup, ThermalSetupError, TransientSettings, VolumeSource,
};
use serde::{Deserialize, Serialize};

use crate::setup_file::shared::{
    check_units, face_defs, face_entries, options_file, required, tolerance, tolerance_file,
};
use crate::setup_file::{FaceEntry, OptionsFile, SetupFileError, ToleranceFile};

/// The whole thermal setup file.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ThermalFile {
    /// Must be `"thermal"`.
    pub analysis: Option<String>,
    /// Length unit; only `"mm"`.
    pub units: Option<String>,
    /// Material preset name.
    pub material: Option<String>,
    /// Default selector tolerances.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tolerance: Option<ToleranceFile>,
    /// Named face references.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub faces: BTreeMap<String, FaceEntry>,
    /// `[[fixed_temperature]]` entries.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub fixed_temperature: Vec<FixedFile>,
    /// `[[heat_flux]]` entries.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub heat_flux: Vec<FluxFile>,
    /// `[[convection]]` entries.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub convection: Vec<ConvectionFile>,
    /// `[[heat_source]]` entries.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub heat_source: Vec<SourceFile>,
    /// `[transient]`: a run in time from a uniform start.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub transient: Option<TransientFile>,
    /// Option overrides.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub options: Option<OptionsFile>,
}

/// `[transient]`.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct TransientFile {
    /// Initial temperature, °C.
    pub initial_c: f64,
    /// Time step, s.
    pub step_s: f64,
    /// Duration, s.
    pub duration_s: f64,
}

/// `[[fixed_temperature]]`.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct FixedFile {
    /// Face names.
    pub faces: Vec<String>,
    /// Temperature, °C.
    pub temperature_c: f64,
}

/// `[[heat_flux]]`.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct FluxFile {
    /// Face names.
    pub faces: Vec<String>,
    /// Total power into the part, W.
    pub power_w: f64,
}

/// Convection faces: a list of names or `"all_other"`.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(untagged)]
pub enum ConvectionFacesFile {
    /// Named faces.
    Named(Vec<String>),
    /// `"all_other"`.
    Keyword(String),
}

/// `[[convection]]`: exactly one of `h` and `preset`.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ConvectionFile {
    /// Face names or `"all_other"`.
    pub faces: ConvectionFacesFile,
    /// Effective coefficient, W/(m²·K).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub h: Option<f64>,
    /// Preset name.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub preset: Option<String>,
    /// Ambient temperature, °C.
    pub ambient_c: f64,
    /// Grey radiation: an emittance in (0, 1] or `"material"`; `h` is
    /// then the convective part alone.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub emittance: Option<EmittanceFile>,
}

/// An emittance: a value, or `"material"` for the material's default.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(untagged)]
pub enum EmittanceFile {
    /// A value.
    Value(f64),
    /// A keyword (`"material"`).
    Keyword(String),
}

/// `[[heat_source]]`.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SourceFile {
    /// Total power, W.
    pub power_w: f64,
    /// The body generating it (a solid name); the whole part when
    /// absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub body: Option<String>,
}

/// Why a thermal file was refused.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ThermalFileError {
    /// A structural-format problem (keys, units, faces).
    #[error(transparent)]
    File(#[from] SetupFileError),
    /// `analysis` is missing or not `"thermal"`.
    #[error("a thermal setup needs analysis = \"thermal\"")]
    NotThermal,
    /// A convection entry gives both or neither of `h` and `preset`.
    #[error("convection {0}: give exactly one of `h` and `preset`")]
    CoefficientCount(usize),
    /// An unknown preset or faces keyword.
    #[error("convection {index}: unknown {what} `{value}`")]
    UnknownKeyword {
        /// Convection entry.
        index: usize,
        /// `preset`, `faces` or `emittance`.
        what: &'static str,
        /// The value given.
        value: String,
    },
    /// The setup parsed but is not valid.
    #[error(transparent)]
    Invalid(#[from] ThermalSetupError),
}

const PRESETS: [(&str, ConvectionPreset); 4] = [
    ("still_air_matte", ConvectionPreset::StillAirMatte),
    (
        "still_air_bright_metal",
        ConvectionPreset::StillAirBrightMetal,
    ),
    ("fan", ConvectionPreset::Fan),
    ("strong_airflow", ConvectionPreset::StrongAirflow),
];

/// The setup-file name of a convection preset.
pub fn preset_name(p: ConvectionPreset) -> &'static str {
    PRESETS
        .iter()
        .find(|(_, q)| *q == p)
        .map_or("preset", |(name, _)| name)
}

fn convection(index: usize, c: &ConvectionFile) -> Result<Convection, ThermalFileError> {
    let faces = match &c.faces {
        ConvectionFacesFile::Named(list) => ConvectionFaces::Named(list.clone()),
        ConvectionFacesFile::Keyword(k) if k == "all_other" => ConvectionFaces::AllOther,
        ConvectionFacesFile::Keyword(k) => {
            return Err(ThermalFileError::UnknownKeyword {
                index,
                what: "faces",
                value: k.clone(),
            });
        }
    };
    let h = match (c.h, &c.preset) {
        (Some(h), None) => HeatTransfer::Value(h),
        (None, Some(p)) => {
            let found = PRESETS.iter().find(|(n, _)| n == p).ok_or_else(|| {
                ThermalFileError::UnknownKeyword {
                    index,
                    what: "preset",
                    value: p.clone(),
                }
            })?;
            HeatTransfer::Preset(found.1)
        }
        _ => return Err(ThermalFileError::CoefficientCount(index)),
    };
    let emittance = match &c.emittance {
        None => None,
        Some(EmittanceFile::Value(v)) => Some(Emittance::Value(*v)),
        Some(EmittanceFile::Keyword(k)) if k == "material" => Some(Emittance::Material),
        Some(EmittanceFile::Keyword(k)) => {
            return Err(ThermalFileError::UnknownKeyword {
                index,
                what: "emittance",
                value: k.clone(),
            });
        }
    };
    Ok(Convection {
        faces,
        h,
        ambient_c: c.ambient_c,
        emittance,
    })
}

/// Parses and validates a thermal setup file's text.
pub fn parse_thermal_setup(text: &str) -> Result<ThermalSetup, ThermalFileError> {
    let f: ThermalFile = toml::from_str(text).map_err(|e| SetupFileError::from_toml(&e))?;
    if f.analysis.as_deref() != Some("thermal") {
        return Err(ThermalFileError::NotThermal);
    }
    check_units(f.units.as_deref())?;
    let setup = ThermalSetup {
        material: required("material", f.material)?,
        tolerance: tolerance(f.tolerance),
        faces: face_defs(&f.faces)?,
        fixed: f
            .fixed_temperature
            .iter()
            .map(|x| FixedTemperature {
                faces: x.faces.clone(),
                temperature_c: x.temperature_c,
            })
            .collect(),
        flux: f
            .heat_flux
            .iter()
            .map(|x| HeatFlux {
                faces: x.faces.clone(),
                power_w: x.power_w,
            })
            .collect(),
        convection: f
            .convection
            .iter()
            .enumerate()
            .map(|(i, c)| convection(i, c))
            .collect::<Result<_, _>>()?,
        sources: f
            .heat_source
            .iter()
            .map(|s| VolumeSource {
                power_w: s.power_w,
                body: s.body.clone(),
            })
            .collect(),
        transient: f.transient.map(|t| TransientSettings {
            initial_c: t.initial_c,
            step_s: t.step_s,
            duration_s: t.duration_s,
        }),
        options: crate::setup_file::options_of(f.options.as_ref()),
    };
    setup.validate()?;
    Ok(setup)
}

/// Writes a thermal setup as file text that [`parse_thermal_setup`]
/// reads back to an equal setup.
pub fn write_thermal_setup(setup: &ThermalSetup) -> Result<String, ThermalFileError> {
    let file = ThermalFile {
        analysis: Some("thermal".to_string()),
        units: Some("mm".to_string()),
        material: Some(setup.material.clone()),
        tolerance: Some(tolerance_file(&setup.tolerance)),
        faces: face_entries(&setup.faces),
        fixed_temperature: setup
            .fixed
            .iter()
            .map(|x| FixedFile {
                faces: x.faces.clone(),
                temperature_c: x.temperature_c,
            })
            .collect(),
        heat_flux: setup
            .flux
            .iter()
            .map(|x| FluxFile {
                faces: x.faces.clone(),
                power_w: x.power_w,
            })
            .collect(),
        convection: setup
            .convection
            .iter()
            .map(|c| ConvectionFile {
                faces: match &c.faces {
                    ConvectionFaces::Named(l) => ConvectionFacesFile::Named(l.clone()),
                    ConvectionFaces::AllOther => {
                        ConvectionFacesFile::Keyword("all_other".to_string())
                    }
                },
                h: match c.h {
                    HeatTransfer::Value(v) => Some(v),
                    HeatTransfer::Preset(_) => None,
                },
                preset: match c.h {
                    HeatTransfer::Preset(p) => PRESETS
                        .iter()
                        .find(|(_, q)| *q == p)
                        .map(|(n, _)| n.to_string()),
                    HeatTransfer::Value(_) => None,
                },
                ambient_c: c.ambient_c,
                emittance: c.emittance.map(|e| match e {
                    Emittance::Value(v) => EmittanceFile::Value(v),
                    Emittance::Material => EmittanceFile::Keyword("material".to_string()),
                }),
            })
            .collect(),
        heat_source: setup
            .sources
            .iter()
            .map(|s| SourceFile {
                power_w: s.power_w,
                body: s.body.clone(),
            })
            .collect(),
        transient: setup.transient.map(|t| TransientFile {
            initial_c: t.initial_c,
            step_s: t.step_s,
            duration_s: t.duration_s,
        }),
        options: Some(options_file(&setup.options)),
    };
    toml::to_string(&file).map_err(|e| {
        SetupFileError::Write {
            detail: e.to_string(),
        }
        .into()
    })
}

#[cfg(test)]
#[path = "thermal_file_tests.rs"]
mod thermal_file_tests;
