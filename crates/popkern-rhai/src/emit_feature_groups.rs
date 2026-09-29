//! Spelling of a fillet with radius groups (PK-602):
//! `.fillet_each([[edges, r0], [edges_1, r1], …], opts)`. The base group
//! is every edge no later group names.

use popkern_core::naming::EntityName;
use popkern_history::feature_blend::FilletRecipeDesc;
use popkern_history::parameter::ParameterId;
use popkern_script::Model;

use crate::emit::EmitError;
use crate::emit_spelling::{Speller, scalar};

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
        super::fillet_options(model, recipe, label, context)?
    ))
}
