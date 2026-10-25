//! Stage 1 material presets (docs/notes/analysis/11 and 13).
//!
//! Printed polymers come from the Ultimaker technical data sheets
//! (2022): one printer, 0.15 mm layers, 100% infill, three build
//! orientations, so every number shares one protocol. E is the XY
//! flexural modulus (ISO 178) and the stiff bound the XY tensile
//! modulus (ASTM D3039). Inter-layer shear (Wang 2023) and compressive
//! strength (Kholil 2022, Valvez 2022) come from other printers and
//! filaments: the nearest measurements, not the TDS protocol. Metals use
//! handbook values and Lienhard, *A Heat Transfer Textbook*, Appendix A.

use super::{Behaviour, Material, MaterialError, Property, Strength};

use Property::{Cited, NotApplicable, Typical, Unsourced};

const TDS_PLA: &str = "Ultimaker PLA TDS (2022)";
const TDS_ABS: &str = "Ultimaker ABS TDS (2022)";
const TDS_PETG: &str = "Ultimaker PETG TDS (2022)";
const TDS_PC: &str = "Ultimaker PC TDS (2022)";
const HANDBOOK: &str = "standard handbook value (notes 11)";

const POLYMER_EMITTANCE: Property = Typical {
    value: 0.9,
    note: "non-metals 0.8-0.95 (Lienhard Table 10.1); printed polymers not listed",
};
const NO_POLYMER_CTE: Property = Unsourced {
    note: "CTE: not in Matsika-Klossa 2024 (PLA, ABS, PA, PEBA only), the Ultimaker TDS or Lienhard; vendor web values only (notes 11)",
};
const CONDUCTIVITY_FFF: &str =
    "Tychanicz-Kwiecien 2025 Table 5, printed at 100% infill, 22 C (ASTM E1530)";
const CTE_FBG: &str = "Matsika-Klossa 2024 Table 2, printed, 0 deg raster, 25-50 C (FBG)";
const NO_COMPRESSIVE: Property = Unsourced {
    note: "compressive strength: no source in the reading set",
};
/// Wang et al. 2023 (eSUN filament, Raise3D Pro3), necking-shaped shear
/// specimens failing in pure shear on the layer plane.
const SHEAR_WANG: &str = "Wang 2023 Table 2, interlayer shear, necking-shaped specimen";
/// Kholil et al. 2022, 0.15 mm layers (the TDS layer height), 100%
/// infill, ASTM D695 loaded across the layers.
const COMPRESSION_KHOLIL: &str = "Kholil 2022 Table 2, yield, 0.15 mm layers, ASTM D695 along Z";

/// Printed PLA. Brittle in-plane: the datasheet says it breaks without
/// yielding upright, and Z elongation is 2%.
// Kholil's σ_c scatter is 3.14 MPa, a measurement, not π.
#[allow(clippy::approx_constant)]
const PLA: Material = Material {
    name: "PLA (FDM)",
    youngs_modulus: Cited {
        value: 3019.0,
        source: "Ultimaker PLA TDS, flexural modulus XY, ISO 178",
    },
    youngs_modulus_stiff_bound: Cited {
        value: 3250.0,
        source: "Ultimaker PLA TDS, tensile modulus XY, ASTM D3039",
    },
    poisson_ratio: Typical {
        value: 0.36,
        note: "usual 0.35-0.36; Krupnin 2023 DIC 0.32",
    },
    strength: Strength {
        inplane: Cited {
            value: 52.5,
            source: "Ultimaker PLA TDS, tensile yield XY",
        },
        inplane_behaviour: Behaviour::Brittle,
        interlayer: Cited {
            value: 33.1,
            source: "Ultimaker PLA TDS, tensile break Z",
        },
        interlayer_shear: Cited {
            value: 25.74,
            source: SHEAR_WANG,
        },
        ultimate_compressive: Cited {
            value: 66.78,
            source: COMPRESSION_KHOLIL,
        },
        inplane_scatter: Cited {
            value: 0.9,
            source: TDS_PLA,
        },
        interlayer_scatter: Cited {
            value: 2.8,
            source: TDS_PLA,
        },
        interlayer_shear_scatter: Cited {
            value: 1.03,
            source: SHEAR_WANG,
        },
        ultimate_compressive_scatter: Cited {
            value: 3.14,
            source: COMPRESSION_KHOLIL,
        },
    },
    density: Cited {
        value: 1240.0,
        source: "Ultimaker PLA TDS, specific gravity 1.24 (ASTM D1505)",
    },
    conductivity: Cited {
        value: 0.182,
        source: CONDUCTIVITY_FFF,
    },
    specific_heat: Unsourced {
        note: "PLA c_p: Blanco 2022 gives DSC fits above Tg only; no room-temperature value sourced",
    },
    emittance_default: POLYMER_EMITTANCE,
    heat_deflection_temp: Cited {
        value: 58.8,
        source: TDS_PLA,
    },
    thermal_expansion: Cited {
        value: 96.0e-6,
        source: CTE_FBG,
    },
};

/// Printed ABS. Yields in-plane; the layer check is brittle.
const ABS: Material = Material {
    name: "ABS (FDM)",
    youngs_modulus: Cited {
        value: 1430.0,
        source: "Ultimaker ABS TDS, flexural modulus XY, ISO 178",
    },
    youngs_modulus_stiff_bound: Cited {
        value: 1962.0,
        source: "Ultimaker ABS TDS, tensile modulus XY, ASTM D3039",
    },
    poisson_ratio: Typical {
        value: 0.35,
        note: "typical for ABS (notes 11)",
    },
    strength: Strength {
        inplane: Cited {
            value: 38.1,
            source: "Ultimaker ABS TDS, tensile yield XY",
        },
        inplane_behaviour: Behaviour::Yielding,
        interlayer: Cited {
            value: 19.0,
            source: "Ultimaker ABS TDS, tensile break Z",
        },
        interlayer_shear: Cited {
            value: 23.01,
            source: SHEAR_WANG,
        },
        ultimate_compressive: Cited {
            value: 36.38,
            source: COMPRESSION_KHOLIL,
        },
        inplane_scatter: Cited {
            value: 0.3,
            source: TDS_ABS,
        },
        interlayer_scatter: Cited {
            value: 0.6,
            source: TDS_ABS,
        },
        interlayer_shear_scatter: Cited {
            value: 1.63,
            source: SHEAR_WANG,
        },
        ultimate_compressive_scatter: Cited {
            value: 2.25,
            source: COMPRESSION_KHOLIL,
        },
    },
    density: Cited {
        value: 1100.0,
        source: "Ultimaker ABS TDS, specific gravity 1.1 (ISO 1183)",
    },
    conductivity: Cited {
        value: 0.14,
        source: "Lienhard Table A.2, ABS 0.14-0.31; lower bound (hotter, conservative for HDT)",
    },
    specific_heat: Unsourced {
        note: "ABS c_p: not in the notes' transcription of Lienhard A.2",
    },
    emittance_default: POLYMER_EMITTANCE,
    heat_deflection_temp: Cited {
        value: 86.6,
        source: TDS_ABS,
    },
    thermal_expansion: Cited {
        value: 86.0e-6,
        source: CTE_FBG,
    },
};

/// Printed PETG. Treated as brittle (Özen 2021: brittle at 2.8%
/// unidirectional; Z elongation 1.8%).
const PETG: Material = Material {
    name: "PETG (FDM)",
    youngs_modulus: Cited {
        value: 1882.0,
        source: "Ultimaker PETG TDS, flexural modulus XY, ISO 178",
    },
    youngs_modulus_stiff_bound: Cited {
        value: 1939.0,
        source: "Ultimaker PETG TDS, tensile modulus XY, ASTM D3039",
    },
    poisson_ratio: Typical {
        value: 0.38,
        note: "typical for PETG (notes 11)",
    },
    strength: Strength {
        inplane: Cited {
            value: 46.2,
            source: "Ultimaker PETG TDS, tensile yield XY",
        },
        inplane_behaviour: Behaviour::Brittle,
        interlayer: Cited {
            value: 19.0,
            source: "Ultimaker PETG TDS, tensile break Z (CoV 34%)",
        },
        interlayer_shear: Unsourced {
            note: "PETG S_layer: Wang 2023 covers PLA, ABS and PC only (notes 11)",
        },
        ultimate_compressive: Cited {
            value: 35.02,
            source: "Valvez 2022 Table 2, yield, 0.40 mm layers, ASTM D695 along Z",
        },
        inplane_scatter: Cited {
            value: 0.8,
            source: TDS_PETG,
        },
        interlayer_scatter: Cited {
            value: 6.4,
            source: TDS_PETG,
        },
        interlayer_shear_scatter: Unsourced { note: "no S_layer" },
        ultimate_compressive_scatter: Cited {
            value: 2.35,
            source: "Valvez 2022 Table 2",
        },
    },
    density: Cited {
        value: 1270.0,
        source: "Ultimaker PETG TDS, specific gravity 1.27 (ISO 1183)",
    },
    conductivity: Cited {
        value: 0.198,
        source: CONDUCTIVITY_FFF,
    },
    specific_heat: Unsourced {
        note: "PETG c_p: no source yet",
    },
    emittance_default: POLYMER_EMITTANCE,
    heat_deflection_temp: Cited {
        value: 76.2,
        source: TDS_PETG,
    },
    thermal_expansion: NO_POLYMER_CTE,
};

/// Printed polycarbonate. Yields in-plane; the layer check is brittle.
const PC: Material = Material {
    name: "PC (FDM)",
    youngs_modulus: Cited {
        value: 1619.0,
        source: "Ultimaker PC TDS, flexural modulus XY, ISO 178",
    },
    youngs_modulus_stiff_bound: Cited {
        value: 2394.0,
        source: "Ultimaker PC TDS, tensile modulus XY, ASTM D3039",
    },
    poisson_ratio: Typical {
        value: 0.37,
        note: "typical for PC (notes 11)",
    },
    strength: Strength {
        inplane: Cited {
            value: 53.3,
            source: "Ultimaker PC TDS, tensile yield XY",
        },
        inplane_behaviour: Behaviour::Yielding,
        interlayer: Cited {
            value: 25.2,
            source: "Ultimaker PC TDS, tensile break Z",
        },
        interlayer_shear: Cited {
            value: 30.16,
            source: SHEAR_WANG,
        },
        ultimate_compressive: NO_COMPRESSIVE,
        inplane_scatter: Cited {
            value: 1.2,
            source: TDS_PC,
        },
        interlayer_scatter: Cited {
            value: 2.1,
            source: TDS_PC,
        },
        interlayer_shear_scatter: Cited {
            value: 1.64,
            source: SHEAR_WANG,
        },
        ultimate_compressive_scatter: Unsourced {
            note: "no compressive strength",
        },
    },
    density: Cited {
        value: 1200.0,
        source: "Lienhard Table A.2, PC (TDS specific gravity 1.18-1.20)",
    },
    conductivity: Cited {
        value: 0.29,
        source: "Lienhard Table A.2, PC",
    },
    specific_heat: Cited {
        value: 1250.0,
        source: "Lienhard Table A.2, PC",
    },
    emittance_default: POLYMER_EMITTANCE,
    heat_deflection_temp: Cited {
        value: 104.5,
        source: TDS_PC,
    },
    thermal_expansion: NO_POLYMER_CTE,
};

/// Aluminium 6061-T6.
const AL_6061_T6: Material = Material {
    name: "Aluminium 6061-T6",
    youngs_modulus: Cited {
        value: 69_000.0,
        source: HANDBOOK,
    },
    youngs_modulus_stiff_bound: NotApplicable,
    poisson_ratio: Cited {
        value: 0.33,
        source: HANDBOOK,
    },
    strength: Strength {
        inplane: Cited {
            value: 276.0,
            source: "standard handbook value, tensile yield (notes 11)",
        },
        inplane_behaviour: Behaviour::Yielding,
        interlayer: NotApplicable,
        interlayer_shear: NotApplicable,
        ultimate_compressive: NO_COMPRESSIVE,
        inplane_scatter: Unsourced {
            note: "handbook typical yield; no scatter given",
        },
        interlayer_scatter: NotApplicable,
        interlayer_shear_scatter: NotApplicable,
        ultimate_compressive_scatter: Unsourced {
            note: "no compressive strength",
        },
    },
    density: Cited {
        value: 2700.0,
        source: "Lienhard Table A.1, Al 6061-T6",
    },
    conductivity: Cited {
        value: 167.0,
        source: "Lienhard Table A.1, Al 6061-T6",
    },
    specific_heat: Cited {
        value: 896.0,
        source: "Lienhard Table A.1, Al 6061-T6",
    },
    emittance_default: Cited {
        value: 0.09,
        source: "Lienhard Table 10.1, commercial aluminium sheet",
    },
    heat_deflection_temp: NotApplicable,
    thermal_expansion: Cited {
        value: 23.0e-6,
        source: "EN 1999-1-1:2007 §3.2.5, aluminium alloys (design value)",
    },
};

/// Structural steel S235 (mild steel).
const STEEL_S235: Material = Material {
    name: "Steel S235",
    youngs_modulus: Cited {
        value: 200_000.0,
        source: HANDBOOK,
    },
    youngs_modulus_stiff_bound: NotApplicable,
    poisson_ratio: Cited {
        value: 0.30,
        source: HANDBOOK,
    },
    strength: Strength {
        inplane: Cited {
            value: 235.0,
            source: "EN 10025 S235 minimum yield (notes 11)",
        },
        inplane_behaviour: Behaviour::Yielding,
        interlayer: NotApplicable,
        interlayer_shear: NotApplicable,
        ultimate_compressive: NO_COMPRESSIVE,
        inplane_scatter: Cited {
            value: 0.0,
            source: "EN 10025 minimum yield: already a lower bound",
        },
        interlayer_scatter: NotApplicable,
        interlayer_shear_scatter: NotApplicable,
        ultimate_compressive_scatter: Unsourced {
            note: "no compressive strength",
        },
    },
    density: Cited {
        value: 7830.0,
        source: "Lienhard Table A.1, AISI 1010 (0.1% C, hot-rolled): nearest low-carbon row to S235",
    },
    conductivity: Cited {
        value: 64.0,
        source: "Lienhard Table A.1, AISI 1010 mild steel",
    },
    specific_heat: Cited {
        value: 434.0,
        source: "Lienhard Table A.1, AISI 1010 (0.1% C, hot-rolled): nearest low-carbon row to S235",
    },
    emittance_default: Cited {
        value: 0.66,
        source: "Lienhard Table 10.1, rolled sheet steel",
    },
    heat_deflection_temp: NotApplicable,
    thermal_expansion: Cited {
        value: 12.0e-6,
        source: "EN 1993-1-1:2005 §3.2.6, structural steel (design value)",
    },
};

/// Every preset, in display order.
pub static PRESETS: [Material; 6] = [PLA, ABS, PETG, PC, AL_6061_T6, STEEL_S235];

/// Materials known but refused, with the reason.
const EXCLUDED: [(&str, &str); 3] = [
    (
        "TPU 95A",
        "hyperelastic; the linear-elastic model is invalid beyond small strain",
    ),
    (
        "TPU (FDM)",
        "hyperelastic; the linear-elastic model is invalid beyond small strain",
    ),
    ("Nylon (FDM)", "no printed-orientation strength source yet"),
];

/// The preset called `name`.
pub fn preset(name: &str) -> Result<&'static Material, MaterialError> {
    if let Some(m) = PRESETS.iter().find(|m| m.name == name) {
        return Ok(m);
    }
    if let Some(&(_, reason)) = EXCLUDED.iter().find(|(n, _)| *n == name) {
        return Err(MaterialError::Excluded {
            name: name.to_string(),
            reason,
        });
    }
    Err(MaterialError::Unknown {
        name: name.to_string(),
        known: PRESETS.iter().map(|m| m.name).collect(),
    })
}
