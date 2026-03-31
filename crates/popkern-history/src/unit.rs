//! Unit system for parametric dimensions.
//!
//! All parameter values are stored internally in SI base units (metres for
//! length, radians for angles). [`Unit`] is metadata used for display
//! conversion only.

use serde::{Deserialize, Serialize};

/// Dimensional category for type-checking expressions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum DimensionKind {
    /// Linear distance (meters, millimeters, inches, etc.).
    Length,
    /// Inverse linear distance (curvature: 1/m).
    InverseLength,
    /// Area (m^2).
    Area,
    /// Volume (m^3).
    Volume,
    /// Angular measurement (radians, degrees).
    Angle,
    /// Second moment of area (m^4).
    MomentOfArea,
    /// Unitless ratio, count, or boolean-like value.
    Dimensionless,
}

/// Display/input unit for a parameter value.
///
/// All values are stored internally in SI base units (meters for length,
/// radians for angles). The `Unit` tag is metadata for display conversion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Unit {
    /// SI base unit for length.
    Meters,
    /// One thousandth of a metre.
    Millimeters,
    /// One hundredth of a metre.
    Centimeters,
    /// Imperial inch (25.4 mm exactly).
    Inches,
    /// SI base unit for angle.
    Radians,
    /// Non-SI angular unit (π/180 radians).
    Degrees,
    /// Dimensionless ratio or pure number.
    Ratio,
    /// Square metres (area).
    SquareMeters,
    /// Cubic metres (volume).
    CubicMeters,
    /// Per metre (inverse length, curvature).
    PerMeter,
    /// Metres to the fourth power (second moment of area).
    MetersToFourth,
}

impl Unit {
    /// Return the dimensional category of this unit.
    pub fn kind(&self) -> DimensionKind {
        match self {
            Unit::Meters | Unit::Millimeters | Unit::Centimeters | Unit::Inches => {
                DimensionKind::Length
            }
            Unit::Radians | Unit::Degrees => DimensionKind::Angle,
            Unit::Ratio => DimensionKind::Dimensionless,
            Unit::SquareMeters => DimensionKind::Area,
            Unit::CubicMeters => DimensionKind::Volume,
            Unit::PerMeter => DimensionKind::InverseLength,
            Unit::MetersToFourth => DimensionKind::MomentOfArea,
        }
    }

    /// Convert a value in this display unit to SI base units.
    pub fn to_si(&self, display_value: f64) -> f64 {
        match self {
            Unit::Meters => display_value,
            Unit::Millimeters => display_value / 1000.0,
            Unit::Centimeters => display_value / 100.0,
            Unit::Inches => display_value * 0.0254,
            Unit::Radians => display_value,
            Unit::Degrees => display_value * std::f64::consts::PI / 180.0,
            Unit::Ratio => display_value,
            Unit::SquareMeters | Unit::CubicMeters | Unit::PerMeter | Unit::MetersToFourth => {
                display_value
            }
        }
    }

    /// Convert a value in SI base units to this display unit.
    pub fn from_si(&self, si_value: f64) -> f64 {
        match self {
            Unit::Meters => si_value,
            Unit::Millimeters => si_value * 1000.0,
            Unit::Centimeters => si_value * 100.0,
            Unit::Inches => si_value / 0.0254,
            Unit::Radians => si_value,
            Unit::Degrees => si_value * 180.0 / std::f64::consts::PI,
            Unit::Ratio => si_value,
            Unit::SquareMeters | Unit::CubicMeters | Unit::PerMeter | Unit::MetersToFourth => {
                si_value
            }
        }
    }

    /// Return the short suffix string for this unit (e.g., `"mm"`, `"deg"`).
    pub fn suffix(&self) -> &'static str {
        match self {
            Unit::Meters => "m",
            Unit::Millimeters => "mm",
            Unit::Centimeters => "cm",
            Unit::Inches => "in",
            Unit::Radians => "rad",
            Unit::Degrees => "deg",
            Unit::Ratio => "",
            Unit::SquareMeters => "m²",
            Unit::CubicMeters => "m³",
            Unit::PerMeter => "1/m",
            Unit::MetersToFourth => "m⁴",
        }
    }
}

#[cfg(test)]
#[path = "unit_tests.rs"]
mod tests;
