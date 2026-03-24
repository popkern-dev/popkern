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
#[derive(Debug, Clone)]
pub struct CartesianPoint {
    /// Label name.
    pub name: String,
    /// Coordinates `[x, y, z]`.
    pub coords: [f64; 3],
}

impl CartesianPoint {
    const TYPE: &str = "CARTESIAN_POINT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        let name = as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?;
        let list = as_real_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?;
        if list.len() < 3 {
            return Err(MappingError {
                entity_type: Self::TYPE.to_string(),
                entity_id: None,
                message: format!("expected 3 coordinates, got {}", list.len()),
            });
        }
        Ok(Self {
            name,
            coords: [list[0], list[1], list[2]],
        })
    }
}

/// DIRECTION — a unit direction vector.
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
        if list.len() < 3 {
            return Err(MappingError {
                entity_type: Self::TYPE.to_string(),
                entity_id: None,
                message: format!("expected 3 ratios, got {}", list.len()),
            });
        }
        Ok(Self {
            name,
            ratios: [list[0], list[1], list[2]],
        })
    }
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
mod tests {
    use super::*;

    fn p_str(s: &str) -> Parameter {
        Parameter::String(s.into())
    }
    fn p_real(v: f64) -> Parameter {
        Parameter::Real(v)
    }
    fn p_ref(id: u64) -> Parameter {
        Parameter::EntityRef(id)
    }
    fn p_enum(s: &str) -> Parameter {
        Parameter::Enum(s.into())
    }
    fn p_int(v: i64) -> Parameter {
        Parameter::Int(v)
    }

    #[test]
    fn cartesian_point_parse() {
        let params = vec![
            p_str("origin"),
            Parameter::List(vec![p_real(1.0), p_real(2.0), p_real(3.0)]),
        ];
        let pt = CartesianPoint::from_params(&params).unwrap();
        assert_eq!(pt.name, "origin");
        assert!((pt.coords[0] - 1.0).abs() < f64::EPSILON);
        assert!((pt.coords[2] - 3.0).abs() < f64::EPSILON);
    }

    #[test]
    fn direction_parse() {
        let params = vec![
            p_str("z"),
            Parameter::List(vec![p_real(0.0), p_real(0.0), p_real(1.0)]),
        ];
        let d = Direction::from_params(&params).unwrap();
        assert_eq!(d.name, "z");
        assert!((d.ratios[2] - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn axis2_placement_3d_parse() {
        let params = vec![p_str(""), p_ref(10), p_ref(11), Parameter::Omitted];
        let a = Axis2Placement3d::from_params(&params).unwrap();
        assert_eq!(a.location, 10);
        assert_eq!(a.axis, Some(11));
        assert_eq!(a.ref_direction, None);
    }

    #[test]
    fn circle_parse() {
        let params = vec![p_str("c1"), p_ref(5), p_real(12.5)];
        let c = StepCircle::from_params(&params).unwrap();
        assert_eq!(c.position, 5);
        assert!((c.radius - 12.5).abs() < f64::EPSILON);
    }

    #[test]
    fn vector_parse() {
        let params = vec![p_str("v1"), p_ref(8), p_real(2.5)];
        let v = StepVector::from_params(&params).unwrap();
        assert_eq!(v.name, "v1");
        assert_eq!(v.orientation, 8);
        assert!((v.magnitude - 2.5).abs() < f64::EPSILON);
    }

    #[test]
    fn axis1_placement_parse() {
        let params = vec![p_str("a1"), p_ref(10), p_ref(11)];
        let a = Axis1Placement::from_params(&params).unwrap();
        assert_eq!(a.name, "a1");
        assert_eq!(a.location, 10);
        assert_eq!(a.axis, Some(11));
    }

    #[test]
    fn axis1_placement_no_axis() {
        let params = vec![p_str("a1"), p_ref(10), Parameter::Omitted];
        let a = Axis1Placement::from_params(&params).unwrap();
        assert_eq!(a.axis, None);
    }

    #[test]
    fn ellipse_parse() {
        let params = vec![p_str("e1"), p_ref(5), p_real(10.0), p_real(6.0)];
        let e = StepEllipse::from_params(&params).unwrap();
        assert_eq!(e.name, "e1");
        assert_eq!(e.position, 5);
        assert!((e.semi_axis_1 - 10.0).abs() < f64::EPSILON);
        assert!((e.semi_axis_2 - 6.0).abs() < f64::EPSILON);
    }

    #[test]
    fn trimmed_curve_parse() {
        let params = vec![
            p_str("tc"),
            p_ref(100),
            Parameter::List(vec![p_ref(50)]),
            Parameter::List(vec![p_real(1.0)]),
            p_enum("T"),
            p_enum("PARAMETER"),
        ];
        let tc = TrimmedCurve::from_params(&params).unwrap();
        assert_eq!(tc.name, "tc");
        assert_eq!(tc.basis_curve, 100);
        assert_eq!(tc.trim_1.len(), 1);
        assert_eq!(tc.trim_2.len(), 1);
        assert!(tc.sense);
        assert_eq!(tc.master_representation, "PARAMETER");
    }

    #[test]
    fn composite_curve_parse() {
        let params = vec![
            p_str("cc"),
            Parameter::List(vec![p_ref(1), p_ref(2), p_ref(3)]),
            p_enum("F"),
        ];
        let cc = CompositeCurve::from_params(&params).unwrap();
        assert_eq!(cc.name, "cc");
        assert_eq!(cc.segments, vec![1, 2, 3]);
        assert!(!cc.self_intersect);
    }

    #[test]
    fn surface_curve_parse() {
        let params = vec![
            p_str("sc"),
            p_ref(10),
            Parameter::List(vec![p_ref(20), p_ref(21)]),
            p_enum("PCURVE_S1"),
        ];
        let sc = SurfaceCurve::from_params(&params).unwrap();
        assert_eq!(sc.name, "sc");
        assert_eq!(sc.curve_3d, 10);
        assert_eq!(sc.associated_geometry, vec![20, 21]);
        assert_eq!(sc.master_representation, "PCURVE_S1");
    }

    #[test]
    fn seam_curve_parse() {
        let params = vec![
            p_str("seam"),
            p_ref(15),
            Parameter::List(vec![p_ref(25)]),
            p_enum("PCURVE_S2"),
        ];
        let sc = SeamCurve::from_params(&params).unwrap();
        assert_eq!(sc.name, "seam");
        assert_eq!(sc.curve_3d, 15);
        assert_eq!(sc.associated_geometry, vec![25]);
        assert_eq!(sc.master_representation, "PCURVE_S2");
    }

    #[test]
    fn pcurve_parse() {
        let params = vec![p_str("pc"), p_ref(40), p_ref(41)];
        let pc = StepPcurve::from_params(&params).unwrap();
        assert_eq!(pc.name, "pc");
        assert_eq!(pc.basis_surface, 40);
        assert_eq!(pc.reference_curve, 41);
    }

    #[test]
    fn bspline_curve_parse() {
        let params = vec![
            p_str("curve"),
            p_int(3),
            Parameter::List(vec![p_ref(1), p_ref(2), p_ref(3), p_ref(4)]),
            p_enum("UNSPECIFIED"),
            p_enum("F"),
            p_enum("F"),
            Parameter::List(vec![p_int(4), p_int(4)]),
            Parameter::List(vec![p_real(0.0), p_real(1.0)]),
            p_enum("UNSPECIFIED"),
        ];
        let bs = BSplineCurveWithKnots::from_params(&params).unwrap();
        assert_eq!(bs.degree, 3);
        assert_eq!(bs.control_points.len(), 4);
        assert!(!bs.closed);
        assert_eq!(bs.knot_multiplicities, vec![4, 4]);
        assert!(bs.weights.is_none());
    }
}
