//! Materials: elastic, strength and thermal properties, each carrying
//! where its number came from.
//!
//! Units are those of a millimetre model with forces in newtons:
//! moduli and strengths in MPa (N/mm²). Thermal properties are in SI
//! (W/m·K, kg/m³, J/kg·K, 1/K, °C) and are converted where they are
//! used. A property no source gives is [`Property::Unsourced`], never a
//! guessed number.

mod presets;

pub use presets::{PRESETS, preset};

/// One material property with its provenance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Property {
    /// A value from a cited source.
    Cited {
        /// The value.
        value: f64,
        /// Source and table or test standard.
        source: &'static str,
    },
    /// A typical value, not measured for this material; reports label
    /// it "typical".
    Typical {
        /// The value.
        value: f64,
        /// Why this value, and from where.
        note: &'static str,
    },
    /// No source gives this property yet.
    Unsourced {
        /// What is missing.
        note: &'static str,
    },
    /// The property does not apply to this material.
    NotApplicable,
}

impl Property {
    /// The value, when there is one.
    pub fn value(&self) -> Option<f64> {
        match self {
            Property::Cited { value, .. } | Property::Typical { value, .. } => Some(*value),
            Property::Unsourced { .. } | Property::NotApplicable => None,
        }
    }

    /// The citation or note, when there is one.
    pub fn source(&self) -> Option<&'static str> {
        match self {
            Property::Cited { source, .. } => Some(source),
            Property::Typical { note, .. } | Property::Unsourced { note } => Some(note),
            Property::NotApplicable => None,
        }
    }
}

/// How a material fails in-plane.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Behaviour {
    /// Marked yield: check von Mises against `inplane`.
    Yielding,
    /// Little elongation: check the largest principal stress.
    Brittle,
}

/// Strength properties.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Strength {
    /// In-plane (XY) strength, MPa: yield where a yield exists, else
    /// break.
    pub inplane: Property,
    /// Which in-plane check applies.
    pub inplane_behaviour: Behaviour,
    /// Inter-layer (upright, Z) break strength, MPa, for printed
    /// polymers; always a brittle check.
    pub interlayer: Property,
    /// Inter-layer shear strength S_layer, MPa, for printed polymers:
    /// shear on the layer planes.
    pub interlayer_shear: Property,
    /// Compressive strength, MPa: the yield (or failure) stress of a
    /// printed compression test, which brittle checks use through
    /// Coulomb–Mohr.
    pub ultimate_compressive: Property,
    /// Scatter of `inplane`, MPa: the datasheet's ± (taken as one
    /// standard deviation); 0 where `inplane` is already a guaranteed
    /// minimum.
    pub inplane_scatter: Property,
    /// Scatter of `interlayer`, MPa, as for `inplane_scatter`.
    pub interlayer_scatter: Property,
    /// Scatter of `interlayer_shear`, MPa (one standard deviation).
    pub interlayer_shear_scatter: Property,
    /// Scatter of `ultimate_compressive`, MPa (one standard deviation).
    pub ultimate_compressive_scatter: Property,
}

/// A material.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Material {
    /// Preset name, as written in a setup file.
    pub name: &'static str,
    /// Young's modulus E, MPa. For printed polymers the XY flexural
    /// modulus: hobby loads are mostly bending (cross-cutting Thread 6).
    pub youngs_modulus: Property,
    /// Stiff bound of the deflection range, MPa: the XY tensile modulus
    /// for printed polymers.
    pub youngs_modulus_stiff_bound: Property,
    /// Poisson's ratio ν.
    pub poisson_ratio: Property,
    /// Strengths.
    pub strength: Strength,
    /// Density, kg/m³.
    pub density: Property,
    /// Thermal conductivity, W/(m·K).
    pub conductivity: Property,
    /// Specific heat capacity c_p, J/(kg·K).
    pub specific_heat: Property,
    /// Default surface emittance ε for the radiation term.
    pub emittance_default: Property,
    /// Heat deflection temperature at 0.455 MPa, °C.
    pub heat_deflection_temp: Property,
    /// Linear thermal expansion coefficient, 1/K.
    pub thermal_expansion: Property,
}

/// Why a material is refused.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum MaterialError {
    /// The material is known but outside the linear-elastic model.
    #[error("material `{name}` is excluded: {reason}")]
    Excluded {
        /// The name asked for.
        name: String,
        /// Why it is excluded.
        reason: &'static str,
    },
    /// No preset has this name.
    #[error("unknown material `{name}`; presets are {known:?}")]
    Unknown {
        /// The name asked for.
        name: String,
        /// The preset names.
        known: Vec<&'static str>,
    },
    /// A required property has no value.
    #[error("material `{material}`: {property} has no value")]
    MissingProperty {
        /// The material.
        material: &'static str,
        /// The property.
        property: &'static str,
    },
    /// A property value is outside its physical range.
    #[error("material `{material}`: {property} = {value} is out of range ({range})")]
    OutOfRange {
        /// The material.
        material: &'static str,
        /// The property.
        property: &'static str,
        /// The value.
        value: f64,
        /// The valid range.
        range: &'static str,
    },
}

impl Material {
    /// Checks the properties a structural solve needs: E > 0,
    /// 0 < ν < 0.5, in-plane strength > 0, and every other present
    /// value finite and positive (emittance in (0, 1]).
    pub fn validate(&self) -> Result<(), MaterialError> {
        let material = self.name;
        let required = |property: &'static str, p: &Property| {
            p.value()
                .ok_or(MaterialError::MissingProperty { material, property })
        };
        let range = |property: &'static str, value: f64, ok: bool, range: &'static str| {
            if ok && value.is_finite() {
                Ok(())
            } else {
                Err(MaterialError::OutOfRange {
                    material,
                    property,
                    value,
                    range,
                })
            }
        };
        let e = required("youngs_modulus", &self.youngs_modulus)?;
        range("youngs_modulus", e, e > 0.0, "> 0")?;
        let nu = required("poisson_ratio", &self.poisson_ratio)?;
        range("poisson_ratio", nu, nu > 0.0 && nu < 0.5, "(0, 0.5)")?;
        let s = required("strength.inplane", &self.strength.inplane)?;
        range("strength.inplane", s, s > 0.0, "> 0")?;
        let optional = [
            (
                "youngs_modulus_stiff_bound",
                self.youngs_modulus_stiff_bound,
            ),
            ("strength.interlayer", self.strength.interlayer),
            ("strength.interlayer_shear", self.strength.interlayer_shear),
            (
                "strength.ultimate_compressive",
                self.strength.ultimate_compressive,
            ),
            ("density", self.density),
            ("conductivity", self.conductivity),
            ("specific_heat", self.specific_heat),
            ("thermal_expansion", self.thermal_expansion),
        ];
        for (property, p) in optional {
            if let Some(v) = p.value() {
                range(property, v, v > 0.0, "> 0")?;
            }
        }
        for (property, p) in [
            ("strength.inplane_scatter", self.strength.inplane_scatter),
            (
                "strength.interlayer_scatter",
                self.strength.interlayer_scatter,
            ),
            (
                "strength.interlayer_shear_scatter",
                self.strength.interlayer_shear_scatter,
            ),
            (
                "strength.ultimate_compressive_scatter",
                self.strength.ultimate_compressive_scatter,
            ),
        ] {
            if let Some(v) = p.value() {
                range(property, v, v >= 0.0, ">= 0")?;
            }
        }
        if let Some(v) = self.emittance_default.value() {
            range("emittance_default", v, v > 0.0 && v <= 1.0, "(0, 1]")?;
        }
        if let Some(v) = self.heat_deflection_temp.value() {
            range("heat_deflection_temp", v, v > -273.15, "> -273.15")?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "material_tests.rs"]
mod material_tests;
