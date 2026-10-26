//! The Rhai spelling of a generative design (PK-956):
//! `.generative(#{ supports: [..], loads: [..], keep: [..], ... })`.

use popkern_history::FeatureTree;
use popkern_history::feature_generative::{GenerativeDesc, GenerativeSupportKindDesc};

use crate::emit::EmitError;
use crate::emit_spelling::{Speller, label_option, options, quoted, scalar, vec3};

/// The call that records `spec`.
pub(super) fn generative_call(
    tree: &FeatureTree,
    context: &str,
    label: &str,
    spell: Speller<'_>,
    spec: &GenerativeDesc,
) -> Result<String, EmitError> {
    let supports = spec
        .supports
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let kind = match s.kind {
                GenerativeSupportKindDesc::Fixed => "fixed",
                GenerativeSupportKindDesc::SurfaceRoller => "surface_roller",
            };
            Ok(format!(
                "#{{ faces: {}, kind: {} }}",
                spell(&format!("supports[{i}]"), &s.faces, true)?,
                quoted(kind)
            ))
        })
        .collect::<Result<Vec<_>, EmitError>>()?;
    let loads = spec
        .loads
        .iter()
        .enumerate()
        .map(|(i, l)| {
            Ok(format!(
                "#{{ faces: {}, force: {}, direction: {} }}",
                spell(&format!("loads[{i}]"), &l.faces, true)?,
                scalar(tree, l.force, context)?,
                vec3(l.direction)
            ))
        })
        .collect::<Result<Vec<_>, EmitError>>()?;
    let mut pairs = vec![
        ("material", quoted(&spec.material)),
        ("supports", format!("[{}]", supports.join(", "))),
        ("loads", format!("[{}]", loads.join(", "))),
    ];
    if !spec.keep.is_empty() {
        let keep = spec
            .keep
            .iter()
            .enumerate()
            .map(|(i, names)| spell(&format!("keep[{i}]"), names, true))
            .collect::<Result<Vec<_>, EmitError>>()?;
        pairs.push(("keep", format!("[{}]", keep.join(", "))));
    }
    pairs.push((
        "volume_fraction",
        scalar(tree, spec.volume_fraction, context)?,
    ));
    pairs.push(("filter", scalar(tree, spec.filter, context)?));
    if let Some(id) = spec.keep_depth {
        pairs.push(("keep_depth", scalar(tree, id, context)?));
    }
    if let Some(id) = spec.cell {
        pairs.push(("cell", scalar(tree, id, context)?));
    }
    if let Some(b) = spec.build_direction {
        pairs.push(("build_direction", quoted(&b.spelling())));
    }
    pairs.push(("max_iterations", spec.max_iterations.to_string()));
    if spec.robust {
        pairs.push(("robust", "true".to_string()));
    }
    if let Some(c) = spec.max_cells {
        pairs.push(("max_cells", c.to_string()));
    }
    pairs.push(label_option(label));
    // `options` renders `, #{ .. }` for a trailing argument; the map is
    // the only one.
    let map = options(pairs);
    Ok(format!(".generative({})", map.trim_start_matches(", ")))
}
