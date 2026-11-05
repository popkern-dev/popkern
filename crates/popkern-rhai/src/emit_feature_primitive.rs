//! Emission of `Feature::Primitive` (PK-1083) as the model-level
//! `m.box(..)` / `m.cylinder(..)` / `m.cone(..)` / `m.sphere(..)` /
//! `m.torus(..)` statement; the placement is spelled only when it is
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
    let mut pairs = vec![label_option(label)];
    let (verb, args) = match &desc.shape {
        PrimitiveShape::Box { x, y, z, centred } => {
            if *centred {
                pairs.push(("centred", "true".to_string()));
            }
            ("box", vec![s(*x)?, s(*y)?, s(*z)?])
        }
        PrimitiveShape::Cylinder { radius, height } => ("cylinder", vec![s(*radius)?, s(*height)?]),
        PrimitiveShape::Cone { r1, r2, height } => ("cone", vec![s(*r1)?, s(*r2)?, s(*height)?]),
        PrimitiveShape::Sphere { radius, latitudes } => {
            if let Some([lo, hi]) = latitudes {
                pairs.push(("lat", format!("[{}, {}]", s(*lo)?, s(*hi)?)));
            }
            ("sphere", vec![s(*radius)?])
        }
        PrimitiveShape::Torus {
            major,
            minor,
            angle,
        } => {
            if let Some(a) = angle {
                pairs.push(("angle", s(*a)?));
            }
            ("torus", vec![s(*major)?, s(*minor)?])
        }
    };
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
