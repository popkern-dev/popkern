//! Every entity name a feature references, as one visitor
//! ([`for_each_entity_name_mut`]): pattern instances remap operation
//! ids through it, and replay rebinds references by provenance through
//! it (PK-677).

use popkern_core::naming::EntityName;

use crate::feature::Feature;

/// Visit every [`EntityName`] a feature references, mutably: sketch
/// planes, faces, edges, axes, end conditions and measure queries.
/// Pattern features reference seed features, not names.
pub(crate) fn for_each_entity_name_mut(feature: &mut Feature, f: &mut dyn FnMut(&mut EntityName)) {
    match feature {
        Feature::Sketch { plane, .. } => {
            visit_sketch_plane(plane, f);
        }
        Feature::Extrude {
            face,
            direction,
            end_condition,
            second_direction,
            ..
        } => {
            f(face);
            if let Some(d) = direction {
                visit_vec3_desc(d, f);
            }
            visit_extrude_end(end_condition, f);
            if let Some(sd) = second_direction {
                visit_vec3_desc(&mut sd.direction, f);
                visit_extrude_end(&mut sd.end_condition, f);
            }
        }
        Feature::Fillet { edges, .. } => {
            for name in edges.iter_mut() {
                f(name);
            }
        }
        Feature::Chamfer { edges, kind, .. } => {
            for name in edges.iter_mut() {
                f(name);
            }
            match kind {
                crate::feature_blend::ChamferKindDesc::Equal => {}
                crate::feature_blend::ChamferKindDesc::TwoDistance { face, .. }
                | crate::feature_blend::ChamferKindDesc::DistanceAngle { face, .. } => {
                    f(face);
                }
            }
        }
        Feature::Boolean { .. }
        | Feature::Thicken { .. }
        | Feature::OffsetSolid { .. }
        | Feature::RawBody { .. } => {}
        Feature::Transform { motion } => {
            if let Some(name) = motion.entity_name_mut() {
                f(name);
            }
        }
        Feature::Shell {
            faces_to_remove, ..
        } => {
            for name in faces_to_remove.iter_mut() {
                f(name);
            }
        }
        Feature::Generative { spec } => {
            for name in spec.entity_names_mut() {
                f(name);
            }
        }
        Feature::Draft {
            faces,
            neutral_plane,
            ..
        } => {
            for name in faces.iter_mut() {
                f(name);
            }
            if let Some(name) = neutral_plane.entity_name_mut() {
                f(name);
            }
        }
        Feature::Sweep {
            face, path_edge, ..
        } => {
            f(face);
            f(path_edge);
        }
        Feature::Loft { sections, .. } => {
            for name in sections.iter_mut() {
                f(name);
            }
        }
        Feature::Pipe { path_edge, .. } => {
            f(path_edge);
        }
        Feature::Fastener { .. }
        | Feature::Hinge { .. }
        | Feature::SeamHinge { .. }
        | Feature::LivingHinge { .. } => {}
        Feature::TappedHole { axis, .. } | Feature::ClearanceHole { axis, .. } => {
            visit_axis_desc(axis, f);
        }
        Feature::HelicalSweep { face, axis, .. } => {
            f(face);
            visit_axis_desc(axis, f);
        }
        Feature::Revolve {
            face,
            axis,
            end_condition,
            ..
        } => {
            f(face);
            visit_axis_desc(axis, f);
            visit_revolve_end(end_condition, f);
        }
        Feature::Reference { query, .. } => {
            visit_measure_query(query, f);
        }
        Feature::LinearPattern { .. }
        | Feature::CircularPattern { .. }
        | Feature::CurvePattern { .. }
        | Feature::Mirror { .. } => {
            // Pattern features reference seed FeatureIds, not entity names
            // that need remapping.
        }
    }
}

/// Visit the entity names inside a [`crate::feature::SketchPlaneDesc`].
fn visit_sketch_plane(
    plane: &mut crate::feature::SketchPlaneDesc,
    f: &mut dyn FnMut(&mut EntityName),
) {
    match plane {
        crate::feature::SketchPlaneDesc::OnFace(name) => f(name),
        crate::feature::SketchPlaneDesc::OnFaceOffset(name, _) => f(name),
        crate::feature::SketchPlaneDesc::XY
        | crate::feature::SketchPlaneDesc::XZ
        | crate::feature::SketchPlaneDesc::YZ
        | crate::feature::SketchPlaneDesc::Custom { .. }
        | crate::feature::SketchPlaneDesc::BaseOffset { .. } => {}
    }
}

/// Visit the entity names inside a [`crate::feature::Vec3Desc`].
fn visit_vec3_desc(desc: &mut crate::feature::Vec3Desc, f: &mut dyn FnMut(&mut EntityName)) {
    match desc {
        crate::feature::Vec3Desc::Literal(_) => {}
        crate::feature::Vec3Desc::FaceNormal(name)
        | crate::feature::Vec3Desc::EdgeDirection(name) => f(name),
        crate::feature::Vec3Desc::Negated(inner) => visit_vec3_desc(inner, f),
    }
}

/// Visit the entity names inside an [`crate::feature::ExtrudeEnd`].
fn visit_extrude_end(end: &mut crate::feature::ExtrudeEnd, f: &mut dyn FnMut(&mut EntityName)) {
    match end {
        crate::feature::ExtrudeEnd::UpToSurface { face } => f(face),
        crate::feature::ExtrudeEnd::OffsetFromSurface { face, .. } => {
            f(face);
        }
        crate::feature::ExtrudeEnd::UpToVertex { vertex } => f(vertex),
        crate::feature::ExtrudeEnd::Blind { .. }
        | crate::feature::ExtrudeEnd::ThroughAll
        | crate::feature::ExtrudeEnd::MidPlane { .. } => {}
    }
}

/// Visit the entity names inside an [`crate::feature::AxisDesc`].
fn visit_axis_desc(axis: &mut crate::feature::AxisDesc, f: &mut dyn FnMut(&mut EntityName)) {
    match axis {
        crate::feature::AxisDesc::Edge(name) => f(name),
        crate::feature::AxisDesc::Literal { .. } => {}
    }
}

/// Visit the entity names inside a [`crate::feature::RevolveEnd`].
fn visit_revolve_end(end: &mut crate::feature::RevolveEnd, f: &mut dyn FnMut(&mut EntityName)) {
    match end {
        crate::feature::RevolveEnd::UpToSurface { face } => f(face),
        crate::feature::RevolveEnd::Blind { .. }
        | crate::feature::RevolveEnd::FullCircle
        | crate::feature::RevolveEnd::MidPlane { .. } => {}
    }
}

/// Visit the entity names inside a [`crate::measure_query::MeasureQuery`].
fn visit_measure_query(
    query: &mut crate::measure_query::MeasureQuery,
    f: &mut dyn FnMut(&mut EntityName),
) {
    use crate::measure_query::MeasureQuery;
    match query {
        MeasureQuery::EdgeLength(n)
        | MeasureQuery::EdgeTangent(n, _)
        | MeasureQuery::FaceArea(n)
        | MeasureQuery::FacePerimeter(n)
        | MeasureQuery::FaceNormal(n, _)
        | MeasureQuery::FaceCentroid(n, _)
        | MeasureQuery::FaceRadius(n)
        | MeasureQuery::FaceMinCurvature(n)
        | MeasureQuery::FaceMaxCurvature(n)
        | MeasureQuery::VertexPosition(n, _) => {
            f(n);
        }
        MeasureQuery::MinDistance(a, b)
        | MeasureQuery::DihedralAngle(a, b)
        | MeasureQuery::EdgeAngle(a, b) => {
            f(a);
            f(b);
        }
        MeasureQuery::BodyVolume
        | MeasureQuery::BodySurfaceArea
        | MeasureQuery::BodyCentroid(_)
        | MeasureQuery::BodyBound(_)
        | MeasureQuery::BodyMomentOfInertia(_) => {}
    }
}
