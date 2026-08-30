//! Parameter-id listing per [`Feature`] variant (scope assignment +
//! qualified-name indexing). Split from `feature.rs` (500-line
//! ratchet, PK-328).

use super::feature::{ExtrudeEnd, Feature, RevolveEnd, SketchPlaneDesc};
use crate::parameter::ParameterId;

/// Extract all [`ParameterId`] fields from a [`Feature`] with their field names.
///
/// Used for automatic scope assignment and qualified-name indexing when a
/// feature is inserted into the timeline.
pub(crate) fn feature_param_ids(feature: &Feature) -> Vec<(&str, ParameterId)> {
    let mut ids: Vec<(&str, ParameterId)> = match feature {
        Feature::Extrude {
            end_condition,
            draft_angle,
            second_direction,
            ..
        } => {
            let mut v = Vec::new();
            match end_condition {
                ExtrudeEnd::Blind { distance } => v.push(("distance", *distance)),
                ExtrudeEnd::OffsetFromSurface { offset, .. } => v.push(("offset", *offset)),
                ExtrudeEnd::MidPlane { half_distance } => v.push(("half_distance", *half_distance)),
                ExtrudeEnd::ThroughAll
                | ExtrudeEnd::UpToSurface { .. }
                | ExtrudeEnd::UpToVertex { .. } => {}
            }
            if let Some(id) = draft_angle {
                v.push(("draft_angle", *id));
            }
            if let Some(sd) = second_direction {
                match &sd.end_condition {
                    ExtrudeEnd::Blind { distance } => v.push(("second_distance", *distance)),
                    ExtrudeEnd::OffsetFromSurface { offset, .. } => {
                        v.push(("second_offset", *offset))
                    }
                    ExtrudeEnd::MidPlane { half_distance } => {
                        v.push(("second_half_distance", *half_distance))
                    }
                    ExtrudeEnd::ThroughAll
                    | ExtrudeEnd::UpToSurface { .. }
                    | ExtrudeEnd::UpToVertex { .. } => {}
                }
                if let Some(id) = sd.draft_angle {
                    v.push(("second_draft_angle", id));
                }
            }
            v
        }
        Feature::Fillet { radius, recipe, .. } => {
            let mut v = vec![("radius", *radius)];
            v.extend(recipe.param_ids());
            v
        }
        Feature::Chamfer { distance, kind, .. } => {
            let mut v = vec![("distance", *distance)];
            v.extend(kind.param_ids());
            v
        }
        Feature::Transform { motion } => motion.param_ids(),
        Feature::Shell { thickness, .. } => vec![("thickness", *thickness)],
        Feature::Thicken { thickness, .. } => vec![("thickness", *thickness)],
        Feature::OffsetSolid { distance, .. } => vec![("distance", *distance)],
        Feature::Draft { angle, .. } => vec![("angle", *angle)],
        Feature::Sweep { twist, .. } => vec![("twist", *twist)],
        Feature::Pipe { radius, .. } => vec![("radius", *radius)],
        Feature::Revolve { end_condition, .. } => match end_condition {
            RevolveEnd::Blind { angle } => vec![("angle", *angle)],
            RevolveEnd::MidPlane { half_angle } => vec![("half_angle", *half_angle)],
            RevolveEnd::FullCircle | RevolveEnd::UpToSurface { .. } => vec![],
        },
        Feature::LinearPattern {
            primary, secondary, ..
        } => {
            let mut v = vec![
                ("primary_count", primary.count),
                ("primary_total_distance", primary.total_distance),
            ];
            if let Some(sec) = secondary {
                v.push(("secondary_count", sec.count));
                v.push(("secondary_total_distance", sec.total_distance));
            }
            v
        }
        Feature::CircularPattern {
            count, total_angle, ..
        } => vec![("count", *count), ("total_angle", *total_angle)],
        Feature::CurvePattern { count, .. } => vec![("count", *count)],
        Feature::Sketch { sketch_data, .. } => sketch_data.param_ids(),
        Feature::Boolean { .. }
        | Feature::Loft { .. }
        | Feature::Reference { .. }
        | Feature::Mirror { .. } => vec![],
    };
    // Handle SketchPlaneDesc::OnFaceOffset parameter.
    if let Feature::Sketch {
        plane: SketchPlaneDesc::OnFaceOffset(_, offset_id),
        ..
    } = feature
    {
        ids.push(("offset", *offset_id));
    }
    ids
}
