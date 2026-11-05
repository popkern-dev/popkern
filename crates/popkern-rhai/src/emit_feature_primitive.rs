//! Emission of `Feature::Primitive` (PK-1083) as the model-level
//! `m.box(..)` / `m.cylinder(..)` / `m.cone(..)` / `m.sphere(..)` /
//! `m.prism(..)` / `m.pyramid(..)` / `m.wedge(..)` / `m.torus(..)`
//! statement; the placement is spelled only when it is
//! not the world frame, its origin in the model's unit.

use popkern_history::FeatureTree;
use popkern_history::feature_primitive::{PrimitiveDesc, PrimitiveFrame, PrimitiveShape};
use popkern_history::unit::Unit;

use super::{Statement, label_option, options};
use crate::emit::EmitError;
use crate::emit_spelling::{scalar, vec3};

/// The primitive's statement.
pub(super) fn primitive_statement(
    tree: &FeatureTree,
    unit: Unit,
    context: &str,
    label: &str,
    desc: &PrimitiveDesc,
) -> Result<Statement, EmitError> {
    let s = |id| scalar(tree, id, context);
    // Wire points are recorded in metres and spelled in the model unit.
    let point = |p: [f64; 3]| vec3(p.map(|c| unit.from_si(c)));
    let mut pairs = vec![label_option(label)];
    let (verb, args) = match &desc.shape {
        PrimitiveShape::Box { x, y, z, centred } => {
            if *centred {
                pairs.push(("centred", "true".to_string()));
            }
            ("box", vec![s(*x)?, s(*y)?, s(*z)?])
        }
        PrimitiveShape::Cylinder { radius, height, .. } => {
            ("cylinder", vec![s(*radius)?, s(*height)?])
        }
        PrimitiveShape::Cone { r1, r2, height, .. } => {
            ("cone", vec![s(*r1)?, s(*r2)?, s(*height)?])
        }
        PrimitiveShape::Sphere {
            radius, latitudes, ..
        } => {
            if let Some([lo, hi]) = latitudes {
                pairs.push(("lat", format!("[{}, {}]", s(*lo)?, s(*hi)?)));
            }
            ("sphere", vec![s(*radius)?])
        }
        PrimitiveShape::Prism {
            sides,
            radius,
            height,
            flats,
        } => {
            if *flats {
                pairs.push(("flats", "true".to_string()));
            }
            ("prism", vec![sides.to_string(), s(*radius)?, s(*height)?])
        }
        PrimitiveShape::Pyramid {
            sides,
            r1,
            r2,
            height,
            flats,
        } => {
            if *flats {
                pairs.push(("flats", "true".to_string()));
            }
            (
                "pyramid",
                vec![sides.to_string(), s(*r1)?, s(*r2)?, s(*height)?],
            )
        }
        PrimitiveShape::Wedge { x, y, z, top_x } => {
            ("wedge", vec![s(*x)?, s(*y)?, s(*z)?, s(*top_x)?])
        }
        PrimitiveShape::Torus { major, minor, .. } => ("torus", vec![s(*major)?, s(*minor)?]),
        PrimitiveShape::WirePolyline { points, closed } => {
            let spelled: Vec<String> = points.iter().map(|p| point(*p)).collect();
            if points.len() == 2 && !closed {
                ("wire_line", spelled)
            } else {
                if *closed {
                    pairs.push(("closed", "true".to_string()));
                }
                ("wire_polyline", vec![format!("[{}]", spelled.join(", "))])
            }
        }
        PrimitiveShape::WireArc { start, via, end } => {
            ("wire_arc", vec![point(*start), point(*via), point(*end)])
        }
        PrimitiveShape::WireCircle { radius } => ("wire_circle", vec![s(*radius)?]),
        PrimitiveShape::WireHelix {
            radius,
            pitch,
            turns,
            left_handed,
        } => {
            if *left_handed {
                pairs.push(("left_hand", "true".to_string()));
            }
            ("wire_helix", vec![s(*radius)?, s(*pitch)?, s(*turns)?])
        }
        PrimitiveShape::SheetRect { x, y, centred } => {
            if *centred {
                pairs.push(("centred", "true".to_string()));
            }
            ("sheet_rect", vec![s(*x)?, s(*y)?])
        }
        PrimitiveShape::SheetDisc { radius } => ("sheet_disc", vec![s(*radius)?]),
        PrimitiveShape::SheetPolygon { points } => {
            let spelled: Vec<String> = points
                .iter()
                .map(|p| {
                    format!(
                        "[{}, {}]",
                        crate::emit_spelling::float(unit.from_si(p[0])),
                        crate::emit_spelling::float(unit.from_si(p[1]))
                    )
                })
                .collect();
            ("sheet_polygon", vec![format!("[{}]", spelled.join(", "))])
        }
        PrimitiveShape::SheetFill { .. } => {
            return Err(EmitError::not_applicable(
                label,
                "sheet_fill",
                "a fill is spelled with its wire (emit_feature)",
            ));
        }
        PrimitiveShape::WireSpline {
            points,
            degree,
            control,
        } => {
            if *degree != 3 {
                pairs.push(("degree", degree.to_string()));
            }
            if *control {
                pairs.push(("control", "true".to_string()));
            }
            let spelled: Vec<String> = points.iter().map(|p| point(*p)).collect();
            ("wire_spline", vec![format!("[{}]", spelled.join(", "))])
        }
    };
    if let Some(a) = desc.shape.angle() {
        pairs.push(("angle", s(a)?));
    }
    let world = PrimitiveFrame::default();
    if desc.frame.origin != world.origin {
        pairs.push(("at", vec3(desc.frame.origin.map(|c| unit.from_si(c)))));
    }
    if desc.frame.axis != world.axis {
        pairs.push(("axis", vec3(desc.frame.axis)));
    }
    if desc.frame.reference != world.reference {
        pairs.push(("ref", vec3(desc.frame.reference)));
    }
    Ok(Statement {
        call: format!("m.{verb}({}{})", args.join(", "), options(pairs)),
        starts_body: true,
    })
}
