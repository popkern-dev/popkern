//! Surface entity definitions from the STEP schema (AP203/AP214).

use crate::step::parser::types::{EntityId, Parameter};
use crate::step::schema::{
    MappingError, as_entity_ref, as_enum, as_int, as_int_list, as_list, as_logical, as_real,
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
            u_closed: as_logical(param_at(params, 5, Self::TYPE)?, Self::TYPE)?,
            v_closed: as_logical(param_at(params, 6, Self::TYPE)?, Self::TYPE)?,
            self_intersect: as_logical(param_at(params, 7, Self::TYPE)?, Self::TYPE)?,
            u_multiplicities: as_int_list(param_at(params, 8, Self::TYPE)?, Self::TYPE)?,
            v_multiplicities: as_int_list(param_at(params, 9, Self::TYPE)?, Self::TYPE)?,
            u_knots: as_real_list(param_at(params, 10, Self::TYPE)?, Self::TYPE)?,
            v_knots: as_real_list(param_at(params, 11, Self::TYPE)?, Self::TYPE)?,
            knot_spec: as_enum(param_at(params, 12, Self::TYPE)?, Self::TYPE)?,
            weights: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step::parser::types::Parameter;

    fn p_str(s: &str) -> Parameter {
        Parameter::String(s.into())
    }
    fn p_real(v: f64) -> Parameter {
        Parameter::Real(v)
    }
    fn p_ref(id: u64) -> Parameter {
        Parameter::EntityRef(id)
    }
    fn p_int(v: i64) -> Parameter {
        Parameter::Int(v)
    }
    fn p_enum(s: &str) -> Parameter {
        Parameter::Enum(s.into())
    }

    #[test]
    fn parse_plane() {
        let params = vec![p_str("base"), p_ref(10)];
        let p = StepPlane::from_params(&params).unwrap();
        assert_eq!(p.name, "base");
        assert_eq!(p.position, 10);
    }

    #[test]
    fn parse_cylindrical_surface() {
        let params = vec![p_str(""), p_ref(42), p_real(5.0)];
        let cs = CylindricalSurface::from_params(&params).unwrap();
        assert_eq!(cs.name, "");
        assert_eq!(cs.position, 42);
        assert!((cs.radius - 5.0).abs() < 1e-12);
    }

    #[test]
    fn parse_conical_surface() {
        let params = vec![p_str("cone"), p_ref(7), p_real(3.0), p_real(0.5)];
        let cs = ConicalSurface::from_params(&params).unwrap();
        assert_eq!(cs.name, "cone");
        assert_eq!(cs.position, 7);
        assert!((cs.radius - 3.0).abs() < 1e-12);
        assert!((cs.semi_angle - 0.5).abs() < 1e-12);
    }

    #[test]
    fn parse_spherical_surface() {
        let params = vec![p_str("ball"), p_ref(20), p_real(12.5)];
        let ss = SphericalSurface::from_params(&params).unwrap();
        assert_eq!(ss.name, "ball");
        assert_eq!(ss.position, 20);
        assert!((ss.radius - 12.5).abs() < 1e-12);
    }

    #[test]
    fn parse_toroidal_surface() {
        let params = vec![p_str("torus"), p_ref(30), p_real(10.0), p_real(2.0)];
        let ts = ToroidalSurface::from_params(&params).unwrap();
        assert_eq!(ts.name, "torus");
        assert_eq!(ts.position, 30);
        assert!((ts.major_radius - 10.0).abs() < 1e-12);
        assert!((ts.minor_radius - 2.0).abs() < 1e-12);
    }

    #[test]
    fn parse_bspline_surface_with_knots() {
        let params = vec![
            p_str("surf"),
            p_int(2),
            p_int(3),
            Parameter::List(vec![
                Parameter::List(vec![p_ref(1), p_ref(2)]),
                Parameter::List(vec![p_ref(3), p_ref(4)]),
            ]),
            p_enum("UNSPECIFIED"),
            p_enum("F"),
            p_enum("F"),
            p_enum("F"),
            Parameter::List(vec![p_int(3), p_int(3)]),
            Parameter::List(vec![p_int(4), p_int(4)]),
            Parameter::List(vec![p_real(0.0), p_real(1.0)]),
            Parameter::List(vec![p_real(0.0), p_real(1.0)]),
            p_enum("UNSPECIFIED"),
        ];
        let bs = BSplineSurfaceWithKnots::from_params(&params).unwrap();
        assert_eq!(bs.name, "surf");
        assert_eq!(bs.u_degree, 2);
        assert_eq!(bs.v_degree, 3);
        assert_eq!(bs.control_points.len(), 2);
        assert_eq!(bs.control_points[0], vec![1, 2]);
        assert_eq!(bs.control_points[1], vec![3, 4]);
        assert_eq!(bs.surface_form, "UNSPECIFIED");
        assert!(!bs.u_closed);
        assert!(!bs.v_closed);
        assert!(!bs.self_intersect);
        assert_eq!(bs.u_multiplicities, vec![3, 3]);
        assert_eq!(bs.v_multiplicities, vec![4, 4]);
        assert_eq!(bs.u_knots, vec![0.0, 1.0]);
        assert_eq!(bs.v_knots, vec![0.0, 1.0]);
        assert_eq!(bs.knot_spec, "UNSPECIFIED");
        assert!(bs.weights.is_none());
    }
}
