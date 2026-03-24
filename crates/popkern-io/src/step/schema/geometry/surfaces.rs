//! Surface entity definitions from the STEP schema (AP203/AP214).

use crate::step::parser::types::{EntityId, Parameter};
use crate::step::schema::{
    MappingError, as_bool, as_entity_ref, as_enum, as_int, as_int_list, as_list, as_real,
    as_real_list, as_ref_list, as_string, param_at,
};

/// PLANE — a planar surface.
#[derive(Debug, Clone)]
pub struct StepPlane {
    /// Label name.
    pub name: String,
    /// Reference to the `AXIS2_PLACEMENT_3D`.
    pub position: EntityId,
}

impl StepPlane {
    const TYPE: &str = "PLANE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            position: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// CYLINDRICAL_SURFACE — a cylindrical surface.
#[derive(Debug, Clone)]
pub struct CylindricalSurface {
    /// Label name.
    pub name: String,
    /// Reference to the `AXIS2_PLACEMENT_3D`.
    pub position: EntityId,
    /// Cylinder radius.
    pub radius: f64,
}

impl CylindricalSurface {
    const TYPE: &str = "CYLINDRICAL_SURFACE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            position: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            radius: as_real(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// CONICAL_SURFACE — a conical surface.
#[derive(Debug, Clone)]
pub struct ConicalSurface {
    /// Label name.
    pub name: String,
    /// Reference to the `AXIS2_PLACEMENT_3D`.
    pub position: EntityId,
    /// Radius at the base.
    pub radius: f64,
    /// Semi-angle of the cone.
    pub semi_angle: f64,
}

impl ConicalSurface {
    const TYPE: &str = "CONICAL_SURFACE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            position: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            radius: as_real(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            semi_angle: as_real(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// SPHERICAL_SURFACE — a spherical surface.
#[derive(Debug, Clone)]
pub struct SphericalSurface {
    /// Label name.
    pub name: String,
    /// Reference to the `AXIS2_PLACEMENT_3D`.
    pub position: EntityId,
    /// Sphere radius.
    pub radius: f64,
}

impl SphericalSurface {
    const TYPE: &str = "SPHERICAL_SURFACE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            position: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            radius: as_real(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// TOROIDAL_SURFACE — a toroidal surface.
#[derive(Debug, Clone)]
pub struct ToroidalSurface {
    /// Label name.
    pub name: String,
    /// Reference to the `AXIS2_PLACEMENT_3D`.
    pub position: EntityId,
    /// Major radius.
    pub major_radius: f64,
    /// Minor radius.
    pub minor_radius: f64,
}

impl ToroidalSurface {
    const TYPE: &str = "TOROIDAL_SURFACE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            position: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            major_radius: as_real(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            minor_radius: as_real(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// B_SPLINE_SURFACE_WITH_KNOTS — a B-spline surface with explicit knot data.
#[derive(Debug, Clone)]
pub struct BSplineSurfaceWithKnots {
    /// Label name.
    pub name: String,
    /// Polynomial degree in the U direction.
    pub u_degree: i64,
    /// Polynomial degree in the V direction.
    pub v_degree: i64,
    /// 2-D grid of control-point entity references.
    pub control_points: Vec<Vec<EntityId>>,
    /// Surface form enumeration.
    pub surface_form: String,
    /// Whether the surface is closed in U.
    pub u_closed: bool,
    /// Whether the surface is closed in V.
    pub v_closed: bool,
    /// Whether the surface self-intersects.
    pub self_intersect: bool,
    /// Knot multiplicities in U.
    pub u_multiplicities: Vec<i64>,
    /// Knot multiplicities in V.
    pub v_multiplicities: Vec<i64>,
    /// Knot values in U.
    pub u_knots: Vec<f64>,
    /// Knot values in V.
    pub v_knots: Vec<f64>,
    /// Knot specification enumeration.
    pub knot_spec: String,
    /// Weights from a `RATIONAL_B_SPLINE_SURFACE` complex part, if present.
    pub weights: Option<Vec<Vec<f64>>>,
}

impl BSplineSurfaceWithKnots {
    const TYPE: &str = "B_SPLINE_SURFACE_WITH_KNOTS";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        let cp_outer = as_list(param_at(params, 3, Self::TYPE)?, Self::TYPE)?;
        let mut control_points = Vec::with_capacity(cp_outer.len());
        for row in cp_outer {
            control_points.push(as_ref_list(row, Self::TYPE)?);
        }
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            u_degree: as_int(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            v_degree: as_int(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            control_points,
            surface_form: as_enum(param_at(params, 4, Self::TYPE)?, Self::TYPE)?,
            u_closed: as_bool(param_at(params, 5, Self::TYPE)?, Self::TYPE)?,
            v_closed: as_bool(param_at(params, 6, Self::TYPE)?, Self::TYPE)?,
            self_intersect: as_bool(param_at(params, 7, Self::TYPE)?, Self::TYPE)?,
            u_multiplicities: as_int_list(param_at(params, 8, Self::TYPE)?, Self::TYPE)?,
            v_multiplicities: as_int_list(param_at(params, 9, Self::TYPE)?, Self::TYPE)?,
            u_knots: as_real_list(param_at(params, 10, Self::TYPE)?, Self::TYPE)?,
            v_knots: as_real_list(param_at(params, 11, Self::TYPE)?, Self::TYPE)?,
            knot_spec: as_enum(param_at(params, 12, Self::TYPE)?, Self::TYPE)?,
            weights: None,
        })
    }
}
