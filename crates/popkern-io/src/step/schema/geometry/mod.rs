//! Typed geometry entity definitions from the STEP schema (AP203/AP214).
//!
//! Curve and point types live here; surface types are in the [`surfaces`] submodule.

mod surfaces;

pub use surfaces::*;

use super::{
    MappingError, as_bool, as_entity_ref, as_enum, as_int, as_int_list, as_list, as_optional_ref,
    as_real, as_real_list, as_ref_list, as_string, param_at,
};
use crate::step::parser::types::{EntityId, Parameter};

/// CARTESIAN_POINT — a point defined by (x, y, z) coordinates.
///
/// ISO 10303-42 gives a point one to three coordinates, as many as its
/// representation context's dimension. A pcurve's
/// `DEFINITIONAL_REPRESENTATION` is a 2-D parametric context, and
/// writers such as OpenCASCADE emit its points with two coordinates.
/// Missing trailing coordinates read as 0, which embeds a 2-D point in
/// the z = 0 plane that popkern pcurves use.
#[derive(Debug, Clone)]
pub struct CartesianPoint {
    /// Label name.
    pub name: String,
    /// Coordinates `[x, y, z]` (trailing coordinates a lower-dimension
    /// point omits are 0).
    pub coords: [f64; 3],
}

impl CartesianPoint {
    const TYPE: &str = "CARTESIAN_POINT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        let name = as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?;
        let list = as_real_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?;
        Ok(Self {
            name,
            coords: padded(&list, 1, Self::TYPE, "coordinates")?,
        })
    }
}

/// DIRECTION — a unit direction vector.
///
/// Two ratios in a 2-D parametric context (a pcurve's line direction),
/// three in model space; a missing `dz` reads as 0 (see
/// [`CartesianPoint`]).
#[derive(Debug, Clone)]
pub struct Direction {
    /// Label name.
    pub name: String,
    /// Direction ratios `[dx, dy, dz]`.
    pub ratios: [f64; 3],
}

impl Direction {
    const TYPE: &str = "DIRECTION";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        let name = as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?;
        let list = as_real_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?;
        Ok(Self {
            name,
            ratios: padded(&list, 2, Self::TYPE, "ratios")?,
        })
    }
}

/// A 1-to-3 component list as `[x, y, z]`, trailing components 0;
/// fewer than `min` or more than 3 components is malformed.
fn padded(list: &[f64], min: usize, entity: &str, what: &str) -> Result<[f64; 3], MappingError> {
    if list.len() < min || list.len() > 3 {
        return Err(MappingError {
            entity_type: entity.to_string(),
            entity_id: None,
            message: format!("expected {min} to 3 {what}, got {}", list.len()),
        });
    }
    let mut out = [0.0; 3];
    out[..list.len()].copy_from_slice(list);
    Ok(out)
}

/// VECTOR — an oriented magnitude.
#[derive(Debug, Clone)]
pub struct StepVector {
    /// Label name.
    pub name: String,
    /// Reference to the `DIRECTION` entity.
    pub orientation: EntityId,
    /// Scalar magnitude.
    pub magnitude: f64,
}

impl StepVector {
    const TYPE: &str = "VECTOR";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            orientation: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            magnitude: as_real(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// AXIS1_PLACEMENT — a single-axis placement in 3-D space.
#[derive(Debug, Clone)]
pub struct Axis1Placement {
    /// Label name.
    pub name: String,
    /// Reference to the location `CARTESIAN_POINT`.
    pub location: EntityId,
    /// Optional reference to the axis `DIRECTION`.
    pub axis: Option<EntityId>,
}

impl Axis1Placement {
    const TYPE: &str = "AXIS1_PLACEMENT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            location: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            axis: as_optional_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// AXIS2_PLACEMENT_3D — a coordinate system defined by location, axis, and reference direction.
#[derive(Debug, Clone)]
pub struct Axis2Placement3d {
    /// Label name.
    pub name: String,
    /// Reference to the location `CARTESIAN_POINT`.
    pub location: EntityId,
    /// Optional reference to the Z-axis `DIRECTION`.
    pub axis: Option<EntityId>,
    /// Optional reference to the X-axis `DIRECTION`.
    pub ref_direction: Option<EntityId>,
}

impl Axis2Placement3d {
    const TYPE: &str = "AXIS2_PLACEMENT_3D";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            location: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            axis: as_optional_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            ref_direction: as_optional_ref(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// AXIS2_PLACEMENT_2D — a 2-D coordinate system (a pcurve's conic
/// placement in a parametric context): location and x direction; the y
/// direction is the x direction turned a quarter turn counterclockwise.
#[derive(Debug, Clone)]
pub struct Axis2Placement2d {
    /// Label name.
    pub name: String,
    /// Reference to the location `CARTESIAN_POINT`.
    pub location: EntityId,
    /// Optional reference to the x-axis `DIRECTION`.
    pub ref_direction: Option<EntityId>,
}

impl Axis2Placement2d {
    const TYPE: &str = "AXIS2_PLACEMENT_2D";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            location: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            ref_direction: as_optional_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// LINE — an infinite line through a point along a direction.
#[derive(Debug, Clone)]
pub struct StepLine {
    /// Label name.
    pub name: String,
    /// Reference to the `CARTESIAN_POINT`.
    pub point: EntityId,
    /// Reference to the `VECTOR` giving the line direction.
    pub direction: EntityId,
}

impl StepLine {
    const TYPE: &str = "LINE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            point: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            direction: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// CIRCLE — a circle defined by a placement and radius.
#[derive(Debug, Clone)]
pub struct StepCircle {
    /// Label name.
    pub name: String,
    /// Reference to the `AXIS2_PLACEMENT_3D`.
    pub position: EntityId,
    /// Circle radius.
    pub radius: f64,
}

impl StepCircle {
    const TYPE: &str = "CIRCLE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            position: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            radius: as_real(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// ELLIPSE — an ellipse defined by a placement and two semi-axes.
#[derive(Debug, Clone)]
pub struct StepEllipse {
    /// Label name.
    pub name: String,
    /// Reference to the `AXIS2_PLACEMENT_3D`.
    pub position: EntityId,
    /// Length of the first semi-axis.
    pub semi_axis_1: f64,
    /// Length of the second semi-axis.
    pub semi_axis_2: f64,
}

impl StepEllipse {
    const TYPE: &str = "ELLIPSE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            position: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            semi_axis_1: as_real(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            semi_axis_2: as_real(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// B_SPLINE_CURVE_WITH_KNOTS — a B-spline curve with explicit knot data.
#[derive(Debug, Clone)]
pub struct BSplineCurveWithKnots {
    /// Label name.
    pub name: String,
    /// Polynomial degree.
    pub degree: i64,
    /// References to control-point `CARTESIAN_POINT` entities.
    pub control_points: Vec<EntityId>,
    /// Curve form enumeration (e.g. `UNSPECIFIED`).
    pub curve_form: String,
    /// Whether the curve is closed.
    pub closed: bool,
    /// Whether the curve self-intersects.
    pub self_intersect: bool,
    /// Knot multiplicity list.
    pub knot_multiplicities: Vec<i64>,
    /// Knot values.
    pub knots: Vec<f64>,
    /// Knot specification enumeration.
    pub knot_spec: String,
    /// Weights from a `RATIONAL_B_SPLINE_CURVE` complex part, if present.
    pub weights: Option<Vec<f64>>,
}

impl BSplineCurveWithKnots {
    const TYPE: &str = "B_SPLINE_CURVE_WITH_KNOTS";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            degree: as_int(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            control_points: as_ref_list(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            curve_form: as_enum(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
            closed: as_bool(param_at(params, 4, Self::TYPE)?, Self::TYPE)?,
            self_intersect: as_bool(param_at(params, 5, Self::TYPE)?, Self::TYPE)?,
            knot_multiplicities: as_int_list(param_at(params, 6, Self::TYPE)?, Self::TYPE)?,
            knots: as_real_list(param_at(params, 7, Self::TYPE)?, Self::TYPE)?,
            knot_spec: as_enum(param_at(params, 8, Self::TYPE)?, Self::TYPE)?,
            weights: None,
        })
    }
}

/// TRIMMED_CURVE — a bounded portion of a basis curve.
#[derive(Debug, Clone)]
pub struct TrimmedCurve {
    /// Label name.
    pub name: String,
    /// Reference to the basis curve entity.
    pub basis_curve: EntityId,
    /// First trim parameter (may be cartesian point ref or parameter value).
    pub trim_1: Vec<Parameter>,
    /// Second trim parameter.
    pub trim_2: Vec<Parameter>,
    /// Sense agreement with the basis curve.
    pub sense: bool,
    /// Master representation enumeration.
    pub master_representation: String,
}

impl TrimmedCurve {
    const TYPE: &str = "TRIMMED_CURVE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        let trim_1 = as_list(param_at(params, 2, Self::TYPE)?, Self::TYPE)?.to_vec();
        let trim_2 = as_list(param_at(params, 3, Self::TYPE)?, Self::TYPE)?.to_vec();
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            basis_curve: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            trim_1,
            trim_2,
            sense: as_bool(param_at(params, 4, Self::TYPE)?, Self::TYPE)?,
            master_representation: as_enum(param_at(params, 5, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// COMPOSITE_CURVE — a curve composed of multiple segments.
#[derive(Debug, Clone)]
pub struct CompositeCurve {
    /// Label name.
    pub name: String,
    /// References to the `COMPOSITE_CURVE_SEGMENT` entities.
    pub segments: Vec<EntityId>,
    /// Whether the curve self-intersects.
    pub self_intersect: bool,
}

impl CompositeCurve {
    const TYPE: &str = "COMPOSITE_CURVE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            segments: as_ref_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            self_intersect: as_bool(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// SURFACE_CURVE — a curve lying on a surface.
#[derive(Debug, Clone)]
pub struct SurfaceCurve {
    /// Label name.
    pub name: String,
    /// Reference to the 3-D curve entity.
    pub curve_3d: EntityId,
    /// References to associated geometry (pcurves or surfaces).
    pub associated_geometry: Vec<EntityId>,
    /// Master representation enumeration.
    pub master_representation: String,
}

impl SurfaceCurve {
    const TYPE: &str = "SURFACE_CURVE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            curve_3d: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            associated_geometry: as_ref_list(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            master_representation: as_enum(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// SEAM_CURVE — a surface curve at a parametric seam (same structure as `SURFACE_CURVE`).
#[derive(Debug, Clone)]
pub struct SeamCurve {
    /// Label name.
    pub name: String,
    /// Reference to the 3-D curve entity.
    pub curve_3d: EntityId,
    /// References to associated geometry.
    pub associated_geometry: Vec<EntityId>,
    /// Master representation enumeration.
    pub master_representation: String,
}

impl SeamCurve {
    const TYPE: &str = "SEAM_CURVE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            curve_3d: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            associated_geometry: as_ref_list(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            master_representation: as_enum(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// PCURVE — a parametric curve defined in the parameter space of a surface.
#[derive(Debug, Clone)]
pub struct StepPcurve {
    /// Label name.
    pub name: String,
    /// Reference to the basis surface.
    pub basis_surface: EntityId,
    /// Reference to the 2-D definitional representation curve.
    pub reference_curve: EntityId,
}

impl StepPcurve {
    const TYPE: &str = "PCURVE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            basis_surface: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            reference_curve: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
