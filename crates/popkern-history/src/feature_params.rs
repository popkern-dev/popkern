//! Parameter-id listing per [`Feature`] variant (scope assignment +
//! qualified-name indexing). Split from `feature.rs` (500-line
//! ratchet, PK-328).

use super::feature::{ExtrudeEnd, Feature, RevolveEnd, SketchPlaneDesc};
use crate::parameter::ParameterId;

/// Extract all [`ParameterId`] fields from a [`Feature`] with their field names.
///
/// Used for automatic scope assignment and qualified-name indexing when a
/// feature is inserted into the timeline.
pub(crate) fn feature_param_ids(feature: &Feature) -> Vec<(String, ParameterId)> {
    let mut ids: Vec<(String, ParameterId)> = static_param_ids(feature)
        .into_iter()
        .map(|(role, id)| (role.to_string(), id))
        .collect();
    if let Feature::Sketch { sketch_data, .. } = feature {
        ids.extend(sketch_data.param_ids());
    }
    ids
}

/// The statically named parameter roles of a feature (everything but
/// the analytic sketch's shape scalars, whose roles depend on the shape).
fn static_param_ids(feature: &Feature) -> Vec<(&'static str, ParameterId)> {
    let mut ids: Vec<(&'static str, ParameterId)> = match feature {
        Feature::RawBody { .. } => vec![],
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
        Feature::Fastener {
            thread,
            length,
            thread_length,
            clearance,
            ..
        } => {
            let mut v = vec![
                ("major_diameter", thread.major_diameter),
                ("pitch", thread.pitch),
                ("length", *length),
            ];
            if let Some(id) = thread_length {
                v.push(("thread_length", *id));
            }
            if let Some(id) = clearance {
                v.push(("clearance", *id));
            }
            v
        }
        Feature::TappedHole {
            thread,
            depth,
            clearance,
            ..
        } => {
            let mut v = vec![
                ("major_diameter", thread.major_diameter),
                ("pitch", thread.pitch),
                ("clearance", *clearance),
            ];
            if let Some(id) = depth {
                v.push(("depth", *id));
            }
            v
        }
        Feature::ClearanceHole {
            thread,
            fit,
            depth,
            countersink,
            counterbore,
            ..
        } => {
            let mut v = vec![
                ("major_diameter", thread.major_diameter),
                ("pitch", thread.pitch),
            ];
            if let crate::feature_fastener::ClearanceFitDesc::Custom(id) = fit {
                v.push(("hole_diameter", *id));
            }
            if let Some(id) = depth {
                v.push(("depth", *id));
            }
            if let Some(id) = countersink {
                v.push(("countersink", *id));
            }
            if let Some(id) = counterbore {
                v.push(("counterbore", *id));
            }
            v
        }
        Feature::HelicalSweep {
            start_radius,
            end_radius,
            pitch,
            turns,
            start_angle,
            ..
        } => vec![
            ("start_radius", *start_radius),
            ("end_radius", *end_radius),
            ("pitch", *pitch),
            ("turns", *turns),
            ("start_angle", *start_angle),
        ],
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
        Feature::Sketch { .. } => vec![],
        Feature::Boolean { .. }
        | Feature::Loft { .. }
        | Feature::Reference { .. }
        | Feature::Mirror { .. } => vec![],
    };
    // Handle SketchPlaneDesc::OnFaceOffset parameter.
    match feature {
        Feature::Sketch {
            plane: SketchPlaneDesc::OnFaceOffset(_, offset_id),
            ..
        } => ids.push(("offset", *offset_id)),
        Feature::Sketch {
            plane: SketchPlaneDesc::BaseOffset { distance, .. },
            ..
        } => ids.push(("plane_offset", *distance)),
        _ => {}
    }
    ids
}

impl Feature {
    /// The operation label of a feature — the same label its structured
    /// errors carry (`ReplayError::FeatureFailed { label }`) and the stem
    /// of its auto-label (`fillet-3`).
    pub fn operation_label(&self) -> &'static str {
        match self {
            Feature::Sketch { .. } => "sketch",
            Feature::RawBody { .. } => "raw_body",
            Feature::Extrude { .. } => "extrude",
            Feature::Fillet { .. } => "fillet",
            Feature::Chamfer { .. } => "chamfer",
            Feature::Boolean { .. } => "boolean",
            Feature::Shell { .. } => "shell",
            Feature::Transform { .. } => "transform",
            Feature::Thicken { .. } => "thicken",
            Feature::OffsetSolid { .. } => "offset_solid",
            Feature::Draft { .. } => "draft",
            Feature::Sweep { .. } => "sweep",
            Feature::Loft { .. } => "loft",
            Feature::Pipe { .. } => "pipe",
            Feature::HelicalSweep { .. } => "helical_sweep",
            Feature::Fastener { kind, .. } => match kind {
                crate::feature_fastener::FastenerKind::HexBolt => "hex_bolt",
                crate::feature_fastener::FastenerKind::SocketCapScrew => "socket_cap_screw",
                crate::feature_fastener::FastenerKind::CountersunkScrew => "countersunk_screw",
                crate::feature_fastener::FastenerKind::ThreadedRod => "threaded_rod",
                crate::feature_fastener::FastenerKind::HexNut => "hex_nut",
            },
            Feature::TappedHole { .. } => "tapped_hole",
            Feature::ClearanceHole { .. } => "clearance_hole",
            Feature::Revolve { .. } => "revolve",
            Feature::Reference { .. } => "reference",
            Feature::LinearPattern { .. } => "linear_pattern",
            Feature::CircularPattern { .. } => "circular_pattern",
            Feature::CurvePattern { .. } => "curve_pattern",
            Feature::Mirror { .. } => "mirror",
        }
    }
}
