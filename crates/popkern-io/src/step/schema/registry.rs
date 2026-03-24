//! Entity type registry — dispatches raw [`EntityInstance`] values to typed schema structs.
//!
//! The central entry point is [`resolve_entity`], which matches on the STEP type name string
//! and calls the appropriate `from_params()` constructor.

use super::{
    MappingError,
    geometry::{
        Axis1Placement, Axis2Placement3d, BSplineCurveWithKnots, BSplineSurfaceWithKnots,
        CartesianPoint, CompositeCurve, ConicalSurface, CylindricalSurface, Direction, SeamCurve,
        SphericalSurface, StepCircle, StepEllipse, StepLine, StepPcurve, StepPlane, StepVector,
        SurfaceCurve, ToroidalSurface, TrimmedCurve,
    },
    presentation::{
        ColourRgb, DraughtingPreDefinedColour, FillAreaStyle, FillAreaStyleColour,
        PresentationLayerAssignment, PresentationStyleAssignment, StyledItem, SurfaceSideStyle,
        SurfaceStyleFillArea, SurfaceStyleUsage,
    },
    product::{
        AdvancedBrepShapeRepresentation, ConversionBasedUnit, FacetedBrepShapeRepresentation,
        GeometricRepresentationContext, GlobalUncertaintyAssignedContext,
        ItemDefinedTransformation, LengthMeasureWithUnit, NextAssemblyUsageOccurrence,
        PlaneAngleMeasureWithUnit, ProductDefinition, ProductDefinitionFormation,
        ProductDefinitionShape, RepresentationRelationship,
        RepresentationRelationshipWithTransformation, ShapeDefinitionRepresentation, SiUnit,
        StepProduct, StepShapeRepresentation, UncertaintyMeasureWithUnit,
    },
    topology::{
        AdvancedFace, BrepWithVoids, ClosedShell, EdgeCurve, EdgeLoop, FaceBound, FaceOuterBound,
        ManifoldSolidBrep, OpenShell, OrientedEdge, VertexPoint,
    },
};
use crate::step::parser::types::EntityInstance;

/// A type-erased resolved STEP entity.
#[derive(Debug)]
pub enum ResolvedEntity {
    // Geometry - points, directions, placements
    /// CARTESIAN_POINT
    CartesianPoint(CartesianPoint),
    /// DIRECTION
    Direction(Direction),
    /// VECTOR
    Vector(StepVector),
    /// AXIS1_PLACEMENT
    Axis1Placement(Axis1Placement),
    /// AXIS2_PLACEMENT_3D
    Axis2Placement3d(Axis2Placement3d),
    // Geometry - curves
    /// LINE
    Line(StepLine),
    /// CIRCLE
    Circle(StepCircle),
    /// ELLIPSE
    Ellipse(StepEllipse),
    /// B_SPLINE_CURVE_WITH_KNOTS
    BSplineCurveWithKnots(BSplineCurveWithKnots),
    /// TRIMMED_CURVE
    TrimmedCurve(TrimmedCurve),
    /// COMPOSITE_CURVE
    CompositeCurve(CompositeCurve),
    /// SURFACE_CURVE
    SurfaceCurve(SurfaceCurve),
    /// SEAM_CURVE
    SeamCurve(SeamCurve),
    /// PCURVE
    Pcurve(StepPcurve),
    // Geometry - surfaces
    /// PLANE
    Plane(StepPlane),
    /// CYLINDRICAL_SURFACE
    CylindricalSurface(CylindricalSurface),
    /// CONICAL_SURFACE
    ConicalSurface(ConicalSurface),
    /// SPHERICAL_SURFACE
    SphericalSurface(SphericalSurface),
    /// TOROIDAL_SURFACE
    ToroidalSurface(ToroidalSurface),
    /// B_SPLINE_SURFACE_WITH_KNOTS
    BSplineSurfaceWithKnots(BSplineSurfaceWithKnots),
    // Topology
    /// VERTEX_POINT
    VertexPoint(VertexPoint),
    /// EDGE_CURVE
    EdgeCurve(EdgeCurve),
    /// ORIENTED_EDGE
    OrientedEdge(OrientedEdge),
    /// EDGE_LOOP
    EdgeLoop(EdgeLoop),
    /// FACE_BOUND
    FaceBound(FaceBound),
    /// FACE_OUTER_BOUND
    FaceOuterBound(FaceOuterBound),
    /// ADVANCED_FACE
    AdvancedFace(AdvancedFace),
    /// CLOSED_SHELL
    ClosedShell(ClosedShell),
    /// OPEN_SHELL
    OpenShell(OpenShell),
    /// MANIFOLD_SOLID_BREP
    ManifoldSolidBrep(ManifoldSolidBrep),
    /// BREP_WITH_VOIDS
    BrepWithVoids(BrepWithVoids),
    // Product
    /// PRODUCT
    Product(StepProduct),
    /// PRODUCT_DEFINITION_FORMATION
    ProductDefinitionFormation(ProductDefinitionFormation),
    /// PRODUCT_DEFINITION
    ProductDefinition(ProductDefinition),
    /// PRODUCT_DEFINITION_SHAPE
    ProductDefinitionShape(ProductDefinitionShape),
    /// SHAPE_DEFINITION_REPRESENTATION
    ShapeDefinitionRepresentation(ShapeDefinitionRepresentation),
    /// SHAPE_REPRESENTATION
    ShapeRepresentation(StepShapeRepresentation),
    /// ADVANCED_BREP_SHAPE_REPRESENTATION
    AdvancedBrepShapeRepresentation(AdvancedBrepShapeRepresentation),
    /// FACETED_BREP_SHAPE_REPRESENTATION
    FacetedBrepShapeRepresentation(FacetedBrepShapeRepresentation),
    /// NEXT_ASSEMBLY_USAGE_OCCURRENCE
    NextAssemblyUsageOccurrence(NextAssemblyUsageOccurrence),
    /// REPRESENTATION_RELATIONSHIP
    RepresentationRelationship(RepresentationRelationship),
    /// REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION
    RepresentationRelationshipWithTransformation(RepresentationRelationshipWithTransformation),
    /// ITEM_DEFINED_TRANSFORMATION
    ItemDefinedTransformation(ItemDefinedTransformation),
    // Units/Context
    /// GEOMETRIC_REPRESENTATION_CONTEXT
    GeometricRepresentationContext(GeometricRepresentationContext),
    /// GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT
    GlobalUncertaintyAssignedContext(GlobalUncertaintyAssignedContext),
    /// UNCERTAINTY_MEASURE_WITH_UNIT
    UncertaintyMeasureWithUnit(UncertaintyMeasureWithUnit),
    /// SI_UNIT
    SiUnit(SiUnit),
    /// CONVERSION_BASED_UNIT
    ConversionBasedUnit(ConversionBasedUnit),
    /// LENGTH_MEASURE_WITH_UNIT
    LengthMeasureWithUnit(LengthMeasureWithUnit),
    /// PLANE_ANGLE_MEASURE_WITH_UNIT
    PlaneAngleMeasureWithUnit(PlaneAngleMeasureWithUnit),
    // Presentation
    /// STYLED_ITEM
    StyledItem(StyledItem),
    /// PRESENTATION_STYLE_ASSIGNMENT
    PresentationStyleAssignment(PresentationStyleAssignment),
    /// SURFACE_STYLE_USAGE
    SurfaceStyleUsage(SurfaceStyleUsage),
    /// SURFACE_SIDE_STYLE
    SurfaceSideStyle(SurfaceSideStyle),
    /// SURFACE_STYLE_FILL_AREA
    SurfaceStyleFillArea(SurfaceStyleFillArea),
    /// FILL_AREA_STYLE
    FillAreaStyle(FillAreaStyle),
    /// FILL_AREA_STYLE_COLOUR
    FillAreaStyleColour(FillAreaStyleColour),
    /// COLOUR_RGB
    ColourRgb(ColourRgb),
    /// DRAUGHTING_PRE_DEFINED_COLOUR
    DraughtingPreDefinedColour(DraughtingPreDefinedColour),
    /// PRESENTATION_LAYER_ASSIGNMENT
    PresentationLayerAssignment(PresentationLayerAssignment),
    /// Unrecognized entity type — preserved for round-trip export.
    Unknown(EntityInstance),
}

/// Resolve a raw [`EntityInstance`] into a typed [`ResolvedEntity`].
///
/// Returns [`ResolvedEntity::Unknown`] for unrecognized type names.
/// Returns `Err(MappingError)` if the type is recognized but parameters are malformed.
pub fn resolve_entity(entity: &EntityInstance) -> Result<ResolvedEntity, MappingError> {
    match entity.type_name.as_str() {
        // Geometry - points, directions, placements
        "CARTESIAN_POINT" => Ok(ResolvedEntity::CartesianPoint(CartesianPoint::from_params(
            &entity.params,
        )?)),
        "DIRECTION" => Ok(ResolvedEntity::Direction(Direction::from_params(
            &entity.params,
        )?)),
        "VECTOR" => Ok(ResolvedEntity::Vector(StepVector::from_params(
            &entity.params,
        )?)),
        "AXIS1_PLACEMENT" => Ok(ResolvedEntity::Axis1Placement(Axis1Placement::from_params(
            &entity.params,
        )?)),
        "AXIS2_PLACEMENT_3D" => Ok(ResolvedEntity::Axis2Placement3d(
            Axis2Placement3d::from_params(&entity.params)?,
        )),
        // Geometry - curves
        "LINE" => Ok(ResolvedEntity::Line(StepLine::from_params(&entity.params)?)),
        "CIRCLE" => Ok(ResolvedEntity::Circle(StepCircle::from_params(
            &entity.params,
        )?)),
        "ELLIPSE" => Ok(ResolvedEntity::Ellipse(StepEllipse::from_params(
            &entity.params,
        )?)),
        "B_SPLINE_CURVE_WITH_KNOTS" => Ok(ResolvedEntity::BSplineCurveWithKnots(
            BSplineCurveWithKnots::from_params(&entity.params)?,
        )),
        "TRIMMED_CURVE" => Ok(ResolvedEntity::TrimmedCurve(TrimmedCurve::from_params(
            &entity.params,
        )?)),
        "COMPOSITE_CURVE" => Ok(ResolvedEntity::CompositeCurve(CompositeCurve::from_params(
            &entity.params,
        )?)),
        "SURFACE_CURVE" => Ok(ResolvedEntity::SurfaceCurve(SurfaceCurve::from_params(
            &entity.params,
        )?)),
        "SEAM_CURVE" => Ok(ResolvedEntity::SeamCurve(SeamCurve::from_params(
            &entity.params,
        )?)),
        "PCURVE" => Ok(ResolvedEntity::Pcurve(StepPcurve::from_params(
            &entity.params,
        )?)),
        // Geometry - surfaces
        "PLANE" => Ok(ResolvedEntity::Plane(StepPlane::from_params(
            &entity.params,
        )?)),
        "CYLINDRICAL_SURFACE" => Ok(ResolvedEntity::CylindricalSurface(
            CylindricalSurface::from_params(&entity.params)?,
        )),
        "CONICAL_SURFACE" => Ok(ResolvedEntity::ConicalSurface(ConicalSurface::from_params(
            &entity.params,
        )?)),
        "SPHERICAL_SURFACE" => Ok(ResolvedEntity::SphericalSurface(
            SphericalSurface::from_params(&entity.params)?,
        )),
        "TOROIDAL_SURFACE" => Ok(ResolvedEntity::ToroidalSurface(
            ToroidalSurface::from_params(&entity.params)?,
        )),
        "B_SPLINE_SURFACE_WITH_KNOTS" => Ok(ResolvedEntity::BSplineSurfaceWithKnots(
            BSplineSurfaceWithKnots::from_params(&entity.params)?,
        )),
        // Topology
        "VERTEX_POINT" => Ok(ResolvedEntity::VertexPoint(VertexPoint::from_params(
            &entity.params,
        )?)),
        "EDGE_CURVE" => Ok(ResolvedEntity::EdgeCurve(EdgeCurve::from_params(
            &entity.params,
        )?)),
        "ORIENTED_EDGE" => Ok(ResolvedEntity::OrientedEdge(OrientedEdge::from_params(
            &entity.params,
        )?)),
        "EDGE_LOOP" => Ok(ResolvedEntity::EdgeLoop(EdgeLoop::from_params(
            &entity.params,
        )?)),
        "FACE_BOUND" => Ok(ResolvedEntity::FaceBound(FaceBound::from_params(
            &entity.params,
        )?)),
        "FACE_OUTER_BOUND" => Ok(ResolvedEntity::FaceOuterBound(FaceOuterBound::from_params(
            &entity.params,
        )?)),
        "ADVANCED_FACE" => Ok(ResolvedEntity::AdvancedFace(AdvancedFace::from_params(
            &entity.params,
        )?)),
        "CLOSED_SHELL" => Ok(ResolvedEntity::ClosedShell(ClosedShell::from_params(
            &entity.params,
        )?)),
        "OPEN_SHELL" => Ok(ResolvedEntity::OpenShell(OpenShell::from_params(
            &entity.params,
        )?)),
        "MANIFOLD_SOLID_BREP" => Ok(ResolvedEntity::ManifoldSolidBrep(
            ManifoldSolidBrep::from_params(&entity.params)?,
        )),
        "BREP_WITH_VOIDS" => Ok(ResolvedEntity::BrepWithVoids(BrepWithVoids::from_params(
            &entity.params,
        )?)),
        // Product
        "PRODUCT" => Ok(ResolvedEntity::Product(StepProduct::from_params(
            &entity.params,
        )?)),
        "PRODUCT_DEFINITION_FORMATION" => Ok(ResolvedEntity::ProductDefinitionFormation(
            ProductDefinitionFormation::from_params(&entity.params)?,
        )),
        "PRODUCT_DEFINITION" => Ok(ResolvedEntity::ProductDefinition(
            ProductDefinition::from_params(&entity.params)?,
        )),
        "PRODUCT_DEFINITION_SHAPE" => Ok(ResolvedEntity::ProductDefinitionShape(
            ProductDefinitionShape::from_params(&entity.params)?,
        )),
        "SHAPE_DEFINITION_REPRESENTATION" => Ok(ResolvedEntity::ShapeDefinitionRepresentation(
            ShapeDefinitionRepresentation::from_params(&entity.params)?,
        )),
        "SHAPE_REPRESENTATION" => Ok(ResolvedEntity::ShapeRepresentation(
            StepShapeRepresentation::from_params(&entity.params)?,
        )),
        "ADVANCED_BREP_SHAPE_REPRESENTATION" => {
            Ok(ResolvedEntity::AdvancedBrepShapeRepresentation(
                AdvancedBrepShapeRepresentation::from_params(&entity.params)?,
            ))
        }
        "FACETED_BREP_SHAPE_REPRESENTATION" => Ok(ResolvedEntity::FacetedBrepShapeRepresentation(
            FacetedBrepShapeRepresentation::from_params(&entity.params)?,
        )),
        "NEXT_ASSEMBLY_USAGE_OCCURRENCE" => Ok(ResolvedEntity::NextAssemblyUsageOccurrence(
            NextAssemblyUsageOccurrence::from_params(&entity.params)?,
        )),
        "REPRESENTATION_RELATIONSHIP" => Ok(ResolvedEntity::RepresentationRelationship(
            RepresentationRelationship::from_params(&entity.params)?,
        )),
        "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION" => Ok(
            ResolvedEntity::RepresentationRelationshipWithTransformation(
                RepresentationRelationshipWithTransformation::from_params(&entity.params)?,
            ),
        ),
        "ITEM_DEFINED_TRANSFORMATION" => Ok(ResolvedEntity::ItemDefinedTransformation(
            ItemDefinedTransformation::from_params(&entity.params)?,
        )),
        // Units/Context
        "GEOMETRIC_REPRESENTATION_CONTEXT" => Ok(ResolvedEntity::GeometricRepresentationContext(
            GeometricRepresentationContext::from_params(&entity.params)?,
        )),
        "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT" => {
            Ok(ResolvedEntity::GlobalUncertaintyAssignedContext(
                GlobalUncertaintyAssignedContext::from_params(&entity.params)?,
            ))
        }
        "UNCERTAINTY_MEASURE_WITH_UNIT" => Ok(ResolvedEntity::UncertaintyMeasureWithUnit(
            UncertaintyMeasureWithUnit::from_params(&entity.params)?,
        )),
        "SI_UNIT" => Ok(ResolvedEntity::SiUnit(SiUnit::from_params(&entity.params)?)),
        "CONVERSION_BASED_UNIT" => Ok(ResolvedEntity::ConversionBasedUnit(
            ConversionBasedUnit::from_params(&entity.params)?,
        )),
        "LENGTH_MEASURE_WITH_UNIT" => Ok(ResolvedEntity::LengthMeasureWithUnit(
            LengthMeasureWithUnit::from_params(&entity.params)?,
        )),
        "PLANE_ANGLE_MEASURE_WITH_UNIT" => Ok(ResolvedEntity::PlaneAngleMeasureWithUnit(
            PlaneAngleMeasureWithUnit::from_params(&entity.params)?,
        )),
        // Presentation
        "STYLED_ITEM" => Ok(ResolvedEntity::StyledItem(StyledItem::from_params(
            &entity.params,
        )?)),
        "PRESENTATION_STYLE_ASSIGNMENT" => Ok(ResolvedEntity::PresentationStyleAssignment(
            PresentationStyleAssignment::from_params(&entity.params)?,
        )),
        "SURFACE_STYLE_USAGE" => Ok(ResolvedEntity::SurfaceStyleUsage(
            SurfaceStyleUsage::from_params(&entity.params)?,
        )),
        "SURFACE_SIDE_STYLE" => Ok(ResolvedEntity::SurfaceSideStyle(
            SurfaceSideStyle::from_params(&entity.params)?,
        )),
        "SURFACE_STYLE_FILL_AREA" => Ok(ResolvedEntity::SurfaceStyleFillArea(
            SurfaceStyleFillArea::from_params(&entity.params)?,
        )),
        "FILL_AREA_STYLE" => Ok(ResolvedEntity::FillAreaStyle(FillAreaStyle::from_params(
            &entity.params,
        )?)),
        "FILL_AREA_STYLE_COLOUR" => Ok(ResolvedEntity::FillAreaStyleColour(
            FillAreaStyleColour::from_params(&entity.params)?,
        )),
        "COLOUR_RGB" => Ok(ResolvedEntity::ColourRgb(ColourRgb::from_params(
            &entity.params,
        )?)),
        "DRAUGHTING_PRE_DEFINED_COLOUR" => Ok(ResolvedEntity::DraughtingPreDefinedColour(
            DraughtingPreDefinedColour::from_params(&entity.params)?,
        )),
        "PRESENTATION_LAYER_ASSIGNMENT" => Ok(ResolvedEntity::PresentationLayerAssignment(
            PresentationLayerAssignment::from_params(&entity.params)?,
        )),
        _ => Ok(ResolvedEntity::Unknown(entity.clone())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step::parser::types::{EntityInstance, Parameter};

    fn make_entity(type_name: &str, params: Vec<Parameter>) -> EntityInstance {
        EntityInstance {
            id: 1,
            type_name: type_name.to_string(),
            params,
            complex_parts: vec![],
        }
    }

    #[test]
    fn resolve_cartesian_point() {
        let entity = make_entity(
            "CARTESIAN_POINT",
            vec![
                Parameter::String("origin".into()),
                Parameter::List(vec![
                    Parameter::Real(1.0),
                    Parameter::Real(2.0),
                    Parameter::Real(3.0),
                ]),
            ],
        );
        let resolved = resolve_entity(&entity).expect("should resolve");
        assert!(
            matches!(resolved, ResolvedEntity::CartesianPoint(p) if p.coords == [1.0, 2.0, 3.0])
        );
    }

    #[test]
    fn resolve_manifold_solid_brep() {
        let entity = make_entity(
            "MANIFOLD_SOLID_BREP",
            vec![Parameter::String("body".into()), Parameter::EntityRef(42)],
        );
        let resolved = resolve_entity(&entity).expect("should resolve");
        assert!(matches!(resolved, ResolvedEntity::ManifoldSolidBrep(b) if b.outer == 42));
    }

    #[test]
    fn resolve_unknown_type() {
        let entity = make_entity("FOOBAR", vec![]);
        let resolved = resolve_entity(&entity).expect("should not error for unknown types");
        assert!(matches!(resolved, ResolvedEntity::Unknown(_)));
    }

    // --- Geometry: points, directions, placements ---

    #[test]
    fn resolve_direction() {
        let entity = make_entity(
            "DIRECTION",
            vec![
                Parameter::String("z".into()),
                Parameter::List(vec![
                    Parameter::Real(0.0),
                    Parameter::Real(0.0),
                    Parameter::Real(1.0),
                ]),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::Direction(_)));
    }

    #[test]
    fn resolve_vector() {
        let entity = make_entity(
            "VECTOR",
            vec![
                Parameter::String("v".into()),
                Parameter::EntityRef(2),
                Parameter::Real(1.0),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::Vector(_)));
    }

    #[test]
    fn resolve_axis1_placement() {
        let entity = make_entity(
            "AXIS1_PLACEMENT",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::EntityRef(2),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::Axis1Placement(_)));
    }

    #[test]
    fn resolve_axis2_placement_3d() {
        let entity = make_entity(
            "AXIS2_PLACEMENT_3D",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::Omitted,
                Parameter::Omitted,
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::Axis2Placement3d(_)));
    }

    // --- Geometry: curves ---

    #[test]
    fn resolve_line() {
        let entity = make_entity(
            "LINE",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::EntityRef(2),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::Line(_)));
    }

    #[test]
    fn resolve_circle() {
        let entity = make_entity(
            "CIRCLE",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::Real(5.0),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::Circle(_)));
    }

    #[test]
    fn resolve_ellipse() {
        let entity = make_entity(
            "ELLIPSE",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::Real(10.0),
                Parameter::Real(5.0),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::Ellipse(_)));
    }

    #[test]
    fn resolve_bspline_curve() {
        let entity = make_entity(
            "B_SPLINE_CURVE_WITH_KNOTS",
            vec![
                Parameter::String("".into()),
                Parameter::Int(3),
                Parameter::List(vec![Parameter::EntityRef(1), Parameter::EntityRef(2)]),
                Parameter::Enum("UNSPECIFIED".into()),
                Parameter::Enum("F".into()),
                Parameter::Enum("F".into()),
                Parameter::List(vec![Parameter::Int(2), Parameter::Int(2)]),
                Parameter::List(vec![Parameter::Real(0.0), Parameter::Real(1.0)]),
                Parameter::Enum("UNSPECIFIED".into()),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::BSplineCurveWithKnots(_)));
    }

    #[test]
    fn resolve_trimmed_curve() {
        let entity = make_entity(
            "TRIMMED_CURVE",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::List(vec![Parameter::Real(0.0)]),
                Parameter::List(vec![Parameter::Real(1.0)]),
                Parameter::Enum("T".into()),
                Parameter::Enum("PARAMETER".into()),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::TrimmedCurve(_)));
    }

    #[test]
    fn resolve_composite_curve() {
        let entity = make_entity(
            "COMPOSITE_CURVE",
            vec![
                Parameter::String("".into()),
                Parameter::List(vec![Parameter::EntityRef(1)]),
                Parameter::Enum("F".into()),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::CompositeCurve(_)));
    }

    #[test]
    fn resolve_surface_curve() {
        let entity = make_entity(
            "SURFACE_CURVE",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::List(vec![Parameter::EntityRef(2)]),
                Parameter::Enum("PCURVE_S1".into()),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::SurfaceCurve(_)));
    }

    #[test]
    fn resolve_seam_curve() {
        let entity = make_entity(
            "SEAM_CURVE",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::List(vec![Parameter::EntityRef(2)]),
                Parameter::Enum("PCURVE_S1".into()),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::SeamCurve(_)));
    }

    #[test]
    fn resolve_pcurve() {
        let entity = make_entity(
            "PCURVE",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::EntityRef(2),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::Pcurve(_)));
    }

    // --- Geometry: surfaces ---

    #[test]
    fn resolve_plane() {
        let entity = make_entity(
            "PLANE",
            vec![Parameter::String("".into()), Parameter::EntityRef(1)],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::Plane(_)));
    }

    #[test]
    fn resolve_cylindrical_surface() {
        let entity = make_entity(
            "CYLINDRICAL_SURFACE",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::Real(5.0),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::CylindricalSurface(_)));
    }

    #[test]
    fn resolve_conical_surface() {
        let entity = make_entity(
            "CONICAL_SURFACE",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::Real(3.0),
                Parameter::Real(0.5),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::ConicalSurface(_)));
    }

    #[test]
    fn resolve_spherical_surface() {
        let entity = make_entity(
            "SPHERICAL_SURFACE",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::Real(10.0),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::SphericalSurface(_)));
    }

    #[test]
    fn resolve_toroidal_surface() {
        let entity = make_entity(
            "TOROIDAL_SURFACE",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::Real(10.0),
                Parameter::Real(2.0),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::ToroidalSurface(_)));
    }

    #[test]
    fn resolve_bspline_surface() {
        let entity = make_entity(
            "B_SPLINE_SURFACE_WITH_KNOTS",
            vec![
                Parameter::String("".into()),
                Parameter::Int(2),
                Parameter::Int(2),
                Parameter::List(vec![Parameter::List(vec![
                    Parameter::EntityRef(1),
                    Parameter::EntityRef(2),
                ])]),
                Parameter::Enum("UNSPECIFIED".into()),
                Parameter::Enum("F".into()),
                Parameter::Enum("F".into()),
                Parameter::Enum("F".into()),
                Parameter::List(vec![Parameter::Int(3), Parameter::Int(3)]),
                Parameter::List(vec![Parameter::Int(3), Parameter::Int(3)]),
                Parameter::List(vec![Parameter::Real(0.0), Parameter::Real(1.0)]),
                Parameter::List(vec![Parameter::Real(0.0), Parameter::Real(1.0)]),
                Parameter::Enum("UNSPECIFIED".into()),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::BSplineSurfaceWithKnots(_)
        ));
    }

    // --- Topology ---

    #[test]
    fn resolve_vertex_point() {
        let entity = make_entity(
            "VERTEX_POINT",
            vec![Parameter::String("".into()), Parameter::EntityRef(1)],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::VertexPoint(_)));
    }

    #[test]
    fn resolve_edge_curve() {
        let entity = make_entity(
            "EDGE_CURVE",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::EntityRef(2),
                Parameter::EntityRef(3),
                Parameter::Enum("T".into()),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::EdgeCurve(_)));
    }

    #[test]
    fn resolve_oriented_edge() {
        let entity = make_entity(
            "ORIENTED_EDGE",
            vec![
                Parameter::String("".into()),
                Parameter::Omitted,
                Parameter::Omitted,
                Parameter::EntityRef(1),
                Parameter::Enum("T".into()),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::OrientedEdge(_)));
    }

    #[test]
    fn resolve_edge_loop() {
        let entity = make_entity(
            "EDGE_LOOP",
            vec![
                Parameter::String("".into()),
                Parameter::List(vec![Parameter::EntityRef(1), Parameter::EntityRef(2)]),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::EdgeLoop(_)));
    }

    #[test]
    fn resolve_face_bound() {
        let entity = make_entity(
            "FACE_BOUND",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::Enum("T".into()),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::FaceBound(_)));
    }

    #[test]
    fn resolve_face_outer_bound() {
        let entity = make_entity(
            "FACE_OUTER_BOUND",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::Enum("T".into()),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::FaceOuterBound(_)));
    }

    #[test]
    fn resolve_advanced_face() {
        let entity = make_entity(
            "ADVANCED_FACE",
            vec![
                Parameter::String("".into()),
                Parameter::List(vec![Parameter::EntityRef(1)]),
                Parameter::EntityRef(2),
                Parameter::Enum("T".into()),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::AdvancedFace(_)));
    }

    #[test]
    fn resolve_closed_shell() {
        let entity = make_entity(
            "CLOSED_SHELL",
            vec![
                Parameter::String("".into()),
                Parameter::List(vec![Parameter::EntityRef(1)]),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::ClosedShell(_)));
    }

    #[test]
    fn resolve_open_shell() {
        let entity = make_entity(
            "OPEN_SHELL",
            vec![
                Parameter::String("".into()),
                Parameter::List(vec![Parameter::EntityRef(1)]),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::OpenShell(_)));
    }

    #[test]
    fn resolve_brep_with_voids() {
        let entity = make_entity(
            "BREP_WITH_VOIDS",
            vec![
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::List(vec![Parameter::EntityRef(2)]),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::BrepWithVoids(_)));
    }

    // --- Product ---

    #[test]
    fn resolve_product() {
        let entity = make_entity(
            "PRODUCT",
            vec![
                Parameter::String("id".into()),
                Parameter::String("name".into()),
                Parameter::String("desc".into()),
                Parameter::List(vec![Parameter::EntityRef(1)]),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::Product(_)));
    }

    #[test]
    fn resolve_product_definition_formation() {
        let entity = make_entity(
            "PRODUCT_DEFINITION_FORMATION",
            vec![
                Parameter::String("v1".into()),
                Parameter::String("desc".into()),
                Parameter::EntityRef(1),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::ProductDefinitionFormation(_)
        ));
    }

    #[test]
    fn resolve_product_definition() {
        let entity = make_entity(
            "PRODUCT_DEFINITION",
            vec![
                Parameter::String("".into()),
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::EntityRef(2),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::ProductDefinition(_)));
    }

    #[test]
    fn resolve_product_definition_shape() {
        let entity = make_entity(
            "PRODUCT_DEFINITION_SHAPE",
            vec![
                Parameter::String("".into()),
                Parameter::String("".into()),
                Parameter::EntityRef(1),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::ProductDefinitionShape(_)
        ));
    }

    #[test]
    fn resolve_shape_definition_representation() {
        let entity = make_entity(
            "SHAPE_DEFINITION_REPRESENTATION",
            vec![Parameter::EntityRef(1), Parameter::EntityRef(2)],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::ShapeDefinitionRepresentation(_)
        ));
    }

    #[test]
    fn resolve_shape_representation() {
        let entity = make_entity(
            "SHAPE_REPRESENTATION",
            vec![
                Parameter::String("".into()),
                Parameter::List(vec![Parameter::EntityRef(1)]),
                Parameter::EntityRef(2),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::ShapeRepresentation(_)));
    }

    #[test]
    fn resolve_advanced_brep_shape_representation() {
        let entity = make_entity(
            "ADVANCED_BREP_SHAPE_REPRESENTATION",
            vec![
                Parameter::String("".into()),
                Parameter::List(vec![Parameter::EntityRef(1)]),
                Parameter::EntityRef(2),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::AdvancedBrepShapeRepresentation(_)
        ));
    }

    #[test]
    fn resolve_faceted_brep_shape_representation() {
        let entity = make_entity(
            "FACETED_BREP_SHAPE_REPRESENTATION",
            vec![
                Parameter::String("".into()),
                Parameter::List(vec![Parameter::EntityRef(1)]),
                Parameter::EntityRef(2),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::FacetedBrepShapeRepresentation(_)
        ));
    }

    #[test]
    fn resolve_next_assembly_usage_occurrence() {
        let entity = make_entity(
            "NEXT_ASSEMBLY_USAGE_OCCURRENCE",
            vec![
                Parameter::String("".into()),
                Parameter::String("".into()),
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::EntityRef(2),
                Parameter::Omitted,
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::NextAssemblyUsageOccurrence(_)
        ));
    }

    #[test]
    fn resolve_representation_relationship() {
        let entity = make_entity(
            "REPRESENTATION_RELATIONSHIP",
            vec![
                Parameter::String("".into()),
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::EntityRef(2),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::RepresentationRelationship(_)
        ));
    }

    #[test]
    fn resolve_representation_relationship_with_transformation() {
        let entity = make_entity(
            "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION",
            vec![
                Parameter::String("".into()),
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::EntityRef(2),
                Parameter::EntityRef(3),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::RepresentationRelationshipWithTransformation(_)
        ));
    }

    #[test]
    fn resolve_item_defined_transformation() {
        let entity = make_entity(
            "ITEM_DEFINED_TRANSFORMATION",
            vec![
                Parameter::String("".into()),
                Parameter::String("".into()),
                Parameter::EntityRef(1),
                Parameter::EntityRef(2),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::ItemDefinedTransformation(_)
        ));
    }

    // --- Units/Context ---

    #[test]
    fn resolve_geometric_representation_context() {
        let entity = make_entity(
            "GEOMETRIC_REPRESENTATION_CONTEXT",
            vec![Parameter::String("".into()), Parameter::Int(3)],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::GeometricRepresentationContext(_)
        ));
    }

    #[test]
    fn resolve_global_uncertainty_assigned_context() {
        let entity = make_entity(
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT",
            vec![
                Parameter::String("".into()),
                Parameter::Int(3),
                Parameter::List(vec![Parameter::EntityRef(1)]),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::GlobalUncertaintyAssignedContext(_)
        ));
    }

    #[test]
    fn resolve_uncertainty_measure_with_unit() {
        let entity = make_entity(
            "UNCERTAINTY_MEASURE_WITH_UNIT",
            vec![
                Parameter::Real(0.001),
                Parameter::EntityRef(1),
                Parameter::String("tol".into()),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::UncertaintyMeasureWithUnit(_)
        ));
    }

    #[test]
    fn resolve_si_unit() {
        let entity = make_entity(
            "SI_UNIT",
            vec![
                Parameter::Derived,
                Parameter::Omitted,
                Parameter::Enum("METRE".into()),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::SiUnit(_)));
    }

    #[test]
    fn resolve_conversion_based_unit() {
        let entity = make_entity(
            "CONVERSION_BASED_UNIT",
            vec![
                Parameter::Omitted,
                Parameter::String("INCH".into()),
                Parameter::EntityRef(1),
            ],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::ConversionBasedUnit(_)));
    }

    #[test]
    fn resolve_length_measure_with_unit() {
        let entity = make_entity(
            "LENGTH_MEASURE_WITH_UNIT",
            vec![Parameter::Real(25.4), Parameter::EntityRef(1)],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(resolved, ResolvedEntity::LengthMeasureWithUnit(_)));
    }

    #[test]
    fn resolve_plane_angle_measure_with_unit() {
        let entity = make_entity(
            "PLANE_ANGLE_MEASURE_WITH_UNIT",
            vec![Parameter::Real(0.01745), Parameter::EntityRef(1)],
        );
        let resolved = resolve_entity(&entity).unwrap();
        assert!(matches!(
            resolved,
            ResolvedEntity::PlaneAngleMeasureWithUnit(_)
        ));
    }

    #[test]
    fn resolve_malformed_params() {
        // CARTESIAN_POINT with no params — missing required name string
        let entity = make_entity("CARTESIAN_POINT", vec![]);
        let result = resolve_entity(&entity);
        assert!(result.is_err(), "malformed params should return Err");
    }
}
