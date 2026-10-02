//! Spelling of a fillet with radius groups (PK-602):
//! `.fillet_each([[edges, r0], [edges_1, r1], …], opts)`. The base group
//! is every edge no later group names. Also the fillet setback options
//! (PK-623/624), shared by both fillet spellings.

use popkern_core::naming::EntityName;
use popkern_history::feature_blend::FilletRecipeDesc;
use popkern_history::parameter::ParameterId;
use popkern_script::Model;

use crate::emit::EmitError;
use crate::emit_spelling::{Speller, float, scalar};

/// The `.fillet_each(...)` call for a grouped fillet feature.
pub(super) fn fillet_each_call(
    model: &Model,
    context: &str,
    spell: Speller<'_>,
    label: &str,
    edges: &[EntityName],
    radius: ParameterId,
    recipe: &FilletRecipeDesc,
) -> Result<String, EmitError> {
    let tree = model.tree();
    let grouped: Vec<&EntityName> = recipe
        .edge_radii
        .iter()
        .flat_map(|g| g.edges.iter())
        .collect();
    let base: Vec<EntityName> = edges
        .iter()
        .filter(|e| !grouped.contains(e))
        .cloned()
        .collect();
    let mut groups = vec![format!(
        "[{}, {}]",
        spell("edges", &base, false)?,
        scalar(tree, radius, context)?
    )];
    for group in &recipe.edge_radii {
        groups.push(format!(
            "[{}, {}]",
            spell("edges", &group.edges, false)?,
            scalar(tree, group.radius, context)?
        ));
    }
    Ok(format!(
        ".fillet_each([{}]{})",
        groups.join(", "),
        super::fillet_options(model, recipe, label, context, spell)?
    ))
}

/// The fillet's setback options (PK-623 `setback`, PK-624
/// `setback_ends`), spelled only when set.
pub(super) fn setback_pairs(
    recipe: &FilletRecipeDesc,
    spell: Speller<'_>,
) -> Result<Vec<(&'static str, String)>, EmitError> {
    let mut pairs = Vec::new();
    if let Some(setback) = recipe.setback {
        pairs.push(("setback", float(setback)));
    }
    if !recipe.end_setbacks.is_empty() {
        let groups = recipe
            .end_setbacks
            .iter()
            .map(|g| {
                Ok(format!(
                    "[{}, {}, {}]",
                    spell("edges", &g.edges, false)?,
                    spell("edges", &g.at, false)?,
                    float(g.distance)
                ))
            })
            .collect::<Result<Vec<_>, EmitError>>()?;
        pairs.push(("setback_ends", format!("[{}]", groups.join(", "))));
    }
    Ok(pairs)
}
