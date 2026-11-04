//! Spellings of the offset family and draft (split from `emit_feature`
//! under the file budget): `.offset_solid(d, #{..})`,
//! `.thicken(t, #{..})` and `.draft(faces, angle, neutral, #{..})`.

use popkern_history::feature::{Feature, JoinStyleKind, TrimPolicyKind};
use popkern_history::feature_draft::NeutralPlaneDesc;
use popkern_script::Model;

use crate::emit::EmitError;
use crate::emit_spelling::{Speller, label_option, options, quoted, scalar, vec3};

/// The `join` and `trim` options an offset carries when not default.
fn offset_pairs(
    label: &str,
    join: &JoinStyleKind,
    trim: &TrimPolicyKind,
) -> Vec<(&'static str, String)> {
    let mut pairs = vec![label_option(label)];
    if matches!(join, JoinStyleKind::Intersection) {
        pairs.push(("join", quoted("intersection")));
    }
    if matches!(trim, TrimPolicyKind::Decline) {
        pairs.push(("trim", quoted("decline")));
    }
    pairs
}

/// The call for an offset or draft `feature`; `None` for any other.
pub(super) fn call(
    model: &Model,
    feature: &Feature,
    label: &str,
    spell: Speller<'_>,
) -> Option<Result<String, EmitError>> {
    if !matches!(
        feature,
        Feature::OffsetSolid { .. } | Feature::Thicken { .. } | Feature::Draft { .. }
    ) {
        return None;
    }
    Some(spelling(model, feature, label, spell))
}

fn spelling(
    model: &Model,
    feature: &Feature,
    label: &str,
    spell: Speller<'_>,
) -> Result<String, EmitError> {
    let tree = model.tree();
    let context = label;
    Ok(match feature {
        Feature::OffsetSolid {
            distance,
            join_style,
            trim_policy,
        } => format!(
            ".offset_solid({}{})",
            scalar(tree, *distance, context)?,
            options(offset_pairs(label, join_style, trim_policy))
        ),
        Feature::Thicken {
            thickness,
            join_style,
            trim_policy,
        } => format!(
            ".thicken({}{})",
            scalar(tree, *thickness, context)?,
            options(offset_pairs(label, join_style, trim_policy))
        ),
        Feature::Draft {
            faces,
            angle,
            neutral_plane,
        } => {
            let neutral = match neutral_plane {
                NeutralPlaneDesc::Face(face) => spell("neutral", std::slice::from_ref(face), true)?,
                // Stored in metres; the script reads the model's unit.
                NeutralPlaneDesc::Literal { origin, normal } => {
                    let unit = model.length_unit();
                    format!(
                        "plane_at({}, {})",
                        vec3(origin.map(|c| unit.from_si(c))),
                        vec3(*normal)
                    )
                }
            };
            format!(
                ".draft({}, {}, {}{})",
                spell("faces", faces, true)?,
                scalar(tree, *angle, context)?,
                neutral,
                options(vec![label_option(label)])
            )
        }
        _ => return Err(EmitError::not_applicable(label, "shaping", "not an offset")),
    })
}
