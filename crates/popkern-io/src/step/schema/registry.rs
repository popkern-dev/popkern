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
#[path = "registry_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "registry_tests_2.rs"]
mod tests_2;
