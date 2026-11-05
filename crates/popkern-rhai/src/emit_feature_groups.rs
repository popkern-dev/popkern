//! Spelling of a fillet with radius groups (PK-602):
//! `.fillet_each([[edges, r0], [edges_1, r1], …], opts)`. The base group
//! is every edge no later group names. Also the fillet setback options
//! (PK-623/624), shared by both fillet spellings.

use popkern_core::naming::EntityName;
use popkern_history::feature_blend::{
    BlendProfileKind, CapKindDesc, ContinuityKind, CornerConstructionDesc, CornerModeDesc,
    FilletRadiusLaw, FilletRecipeDesc, InterferencePolicyKind,
};
use popkern_history::parameter::ParameterId;
use popkern_script::Model;

use crate::emit::EmitError;
use crate::emit_spelling::{Speller, float, label_option, options, quoted, scalar};

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
        fillet_options(model, recipe, label, context, spell)?
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

/// The options map of a fillet: its label and every recipe choice
/// that differs from the default.
pub(super) fn fillet_options(
    model: &Model,
    recipe: &FilletRecipeDesc,
    label: &str,
    context: &str,
    spell: Speller<'_>,
) -> Result<String, EmitError> {
    let tree = model.tree();
    let mut pairs = vec![label_option(label)];
    match recipe.continuity {
        ContinuityKind::G1 => {}
        ContinuityKind::G0 => pairs.push(("continuity", quoted("g0"))),
        ContinuityKind::G2 => pairs.push(("continuity", quoted("g2"))),
    }
    match recipe.cap {
        CapKindDesc::Open => {}
        CapKindDesc::Flat => pairs.push(("cap", quoted("flat"))),
        CapKindDesc::Rounded => pairs.push(("cap", quoted("rounded"))),
        CapKindDesc::TangentMatch => pairs.push(("cap", quoted("tangent_match"))),
    }
    match recipe.profile {
        BlendProfileKind::Circular => {}
        BlendProfileKind::Conic { rho } => {
            pairs.push(("profile", format!("conic({})", float(rho))))
        }
        BlendProfileKind::CurvatureContinuous => {
            pairs.push(("profile", quoted("curvature_continuous")))
        }
    }
    if recipe.interference == InterferencePolicyKind::ServeNonInterfering {
        pairs.push(("interference", quoted("serve_non_interfering")));
    }
    if recipe.corner == CornerModeDesc::Mitre {
        pairs.push(("corner", quoted("mitre")));
    }
    match recipe.corner_cap.construction {
        CornerConstructionDesc::Auto => {}
        CornerConstructionDesc::Setback => pairs.push(("corner_cap", quoted("setback"))),
        CornerConstructionDesc::Fitted => pairs.push(("corner_cap", quoted("fitted"))),
        CornerConstructionDesc::Classic => pairs.push(("corner_cap", quoted("classic"))),
        CornerConstructionDesc::RollingBall => pairs.push(("corner_cap", quoted("rolling_ball"))),
    }
    if recipe.corner_cap.strict {
        pairs.push(("corner_cap_strict", "true".to_string()));
    }
    pairs.extend(setback_pairs(recipe, spell)?);
    match &recipe.radius_law {
        FilletRadiusLaw::Constant => {}
        FilletRadiusLaw::Linear { end } => {
            pairs.push(("law", format!("linear({})", scalar(tree, *end, context)?)))
        }
        FilletRadiusLaw::ByStations(stations) => {
            let items: Vec<String> = stations
                .iter()
                .map(|s| {
                    Ok(format!(
                        "[{}, {}]",
                        float(s.t),
                        scalar(tree, s.radius, context)?
                    ))
                })
                .collect::<Result<_, EmitError>>()?;
            pairs.push(("law", format!("by_stations([{}])", items.join(", "))));
        }
    }
    Ok(options(pairs))
}
