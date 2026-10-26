//! Spellings the canonical emitter needs: scalar bindings, selection
//! queries, planes, and profile shapes, each in the Rhai dialect of
//! the surface spec.

use popkern_history::FeatureTree;
use popkern_history::expression::Expression;
use popkern_history::feature::{BasePlane, SketchPlaneDesc};
use popkern_history::feature_sketch::ScalarBinding;
use popkern_history::parameter::ParameterId;
use popkern_script::Model;
use popkern_sketch::profile::{AnalyticProfile, PathSegment, ProfileShape};

use crate::emit::EmitError;

/// A float literal that Rhai reads back as the same `f64`.
pub fn float(v: f64) -> String {
    let text = format!("{v:?}");
    if text.contains(['.', 'e', 'n', 'i']) {
        text
    } else {
        format!("{text}.0")
    }
}

/// A quoted Rhai string literal.
pub fn quoted(text: &str) -> String {
    format!("{text:?}")
}

/// An options map argument (`, #{ k: v, .. }`), or nothing when empty.
pub fn options(pairs: Vec<(&str, String)>) -> String {
    if pairs.is_empty() {
        String::new()
    } else {
        let inner: Vec<String> = pairs.iter().map(|(k, v)| format!("{k}: {v}")).collect();
        format!(", #{{ {} }}", inner.join(", "))
    }
}

/// The `label` option pair.
pub fn label_option(label: &str) -> (&'static str, String) {
    ("label", quoted(label))
}

/// The scalar argument that re-binds parameter `id`: `p("name")` for
/// a parameter reference, `expr("…")` for a driven parameter, else the
/// literal in the parameter's display unit.
pub fn scalar(tree: &FeatureTree, id: ParameterId, context: &str) -> Result<String, EmitError> {
    let table = tree.parameters();
    let param = table.get(id).ok_or_else(|| {
        EmitError::not_applicable(
            context,
            "parameter",
            "a feature parameter is missing from the table",
        )
    })?;
    match &param.expression {
        Some(Expression::Param(target)) => {
            let name = table.qualified_name(*target).ok_or_else(|| {
                EmitError::not_applicable(context, "parameter", "unnamed parameter reference")
            })?;
            Ok(format!("p({})", quoted(&name)))
        }
        Some(expression) => Ok(format!(
            "expr({})",
            quoted(&popkern_history::format_expression(expression, table))
        )),
        None => Ok(float(param.unit.from_si(param.value))),
    }
}

/// A scalar binding of a profile.
pub fn binding(tree: &FeatureTree, b: &ScalarBinding, context: &str) -> Result<String, EmitError> {
    match b {
        ScalarBinding::Literal(v) => Ok(float(*v)),
        ScalarBinding::Parameter(id) => scalar(tree, *id, context),
    }
}

/// A profile coordinate that the skin only accepts as a literal
/// (polygon points, arc via/to).
fn literal_only(tree: &FeatureTree, b: &ScalarBinding, context: &str) -> Result<String, EmitError> {
    match b {
        ScalarBinding::Literal(v) => Ok(float(*v)),
        ScalarBinding::Parameter(id) => {
            let param = tree.parameters().get(*id).ok_or_else(|| {
                EmitError::not_applicable(context, "sketch", "missing profile parameter")
            })?;
            if param.expression.is_some() {
                return Err(EmitError::not_applicable(
                    context,
                    "sketch",
                    "polygon and arc coordinates are literals in the surface; a driven \
                     coordinate has no spelling",
                ));
            }
            Ok(float(param.unit.from_si(param.value)))
        }
    }
}

/// The Rhai form of a recorded selection query (`describe` text):
/// enum words become strings. Opaque and Debug-form `named(...)`
/// queries have no spelling.
pub fn selection(describe: &str, context: &str) -> Result<String, EmitError> {
    if describe.contains("named(EntityName") || describe.contains("opaque") {
        return Err(EmitError::not_applicable(
            context,
            "selection",
            format!("selection '{describe}' was recorded without a re-spellable query"),
        ));
    }
    const ENUM_ARGS: [&str; 5] = ["between(", "direction(", "kind(", "facing(", "facing_neg("];
    let mut out = String::with_capacity(describe.len() + 16);
    let mut rest = describe;
    while !rest.is_empty() {
        let hit = ENUM_ARGS
            .iter()
            .filter_map(|name| rest.find(name).map(|i| (i, *name)))
            .filter(|(i, name)| {
                // `facing(` must not match inside `facing_neg(`; require a
                // '.' before the name.
                *i == 0 || rest.as_bytes()[*i - 1] == b'.' || name.starts_with('.')
            })
            .min_by_key(|(i, _)| *i);
        let Some((i, name)) = hit else {
            out.push_str(rest);
            break;
        };
        let start = i + name.len();
        let close = rest[start..].find(')').map(|c| start + c).ok_or_else(|| {
            EmitError::not_applicable(context, "selection", "unbalanced selection text")
        })?;
        out.push_str(&rest[..start]);
        let args: Vec<String> = rest[start..close]
            .split(',')
            .map(|a| {
                let a = a.trim();
                if a.starts_with('"') || a.parse::<f64>().is_ok() {
                    a.to_string()
                } else {
                    quoted(a)
                }
            })
            .collect();
        out.push_str(&args.join(", "));
        out.push(')');
        rest = &rest[close + 1..];
    }
    Ok(out)
}

/// A selection spelled as its name snapshot: `edges().named("a")
/// .or(edges().named("b"))`. Face selections have no disjunction in
/// the grammar, so only a single face can be spelled this way.
pub fn named_chain(
    model: &Model,
    names: &[popkern_core::naming::EntityName],
    faces: bool,
    context: &str,
) -> Result<String, EmitError> {
    if faces && names.len() != 1 {
        return Err(EmitError::not_applicable(
            context,
            "selection",
            format!(
                "a stale face query over {} faces cannot be spelled: the grammar has no \
                 face disjunction",
                names.len()
            ),
        ));
    }
    let kind = if faces { "faces" } else { "edges" };
    let mut spelled = Vec::with_capacity(names.len());
    for name in names {
        let text = model
            .spell_entity(name)
            .map_err(|e| EmitError::not_applicable(context, "selection", e.to_string()))?;
        spelled.push(format!("{kind}().named({})", quoted(&text)));
    }
    let mut out = spelled
        .first()
        .cloned()
        .ok_or_else(|| EmitError::not_applicable(context, "selection", "empty selection"))?;
    for part in spelled.iter().skip(1) {
        out.push_str(&format!(".or({part})"));
    }
    Ok(out)
}

fn base_at(base: BasePlane) -> &'static str {
    match base {
        BasePlane::XY => "xy_at",
        BasePlane::XZ => "xz_at",
        BasePlane::YZ => "yz_at",
    }
}

/// Spells a recorded selection: `(role, names, faces)` → query text.
pub type Speller<'a> =
    &'a dyn Fn(&str, &[popkern_core::naming::EntityName], bool) -> Result<String, EmitError>;

/// A sketch plane.
pub fn plane(
    model: &Model,
    plane: &SketchPlaneDesc,
    spell: Speller<'_>,
    context: &str,
) -> Result<String, EmitError> {
    let tree = model.tree();
    Ok(match plane {
        SketchPlaneDesc::XY => "xy()".into(),
        SketchPlaneDesc::XZ => "xz()".into(),
        SketchPlaneDesc::YZ => "yz()".into(),
        SketchPlaneDesc::BaseOffset { base, distance } => {
            format!("{}({})", base_at(*base), scalar(tree, *distance, context)?)
        }
        SketchPlaneDesc::OnFace(face) => {
            format!(
                "on_face({})",
                spell("plane", std::slice::from_ref(face), true)?
            )
        }
        SketchPlaneDesc::OnFaceOffset(face, distance) => format!(
            "on_face_offset({}, {})",
            spell("plane", std::slice::from_ref(face), true)?,
            scalar(tree, *distance, context)?
        ),
        SketchPlaneDesc::Custom { origin, normal } => {
            format!("plane_at({}, {})", vec3(*origin), vec3(*normal))
        }
    })
}

/// `[x, y, z]`.
pub fn vec3(v: [f64; 3]) -> String {
    format!("[{}, {}, {}]", float(v[0]), float(v[1]), float(v[2]))
}

/// The profile builder chain after `m.profile(plane)`.
pub fn shape(
    tree: &FeatureTree,
    profile: &AnalyticProfile<ScalarBinding>,
    context: &str,
) -> Result<String, EmitError> {
    let mut out = String::new();
    match &profile.shape {
        ProfileShape::Rect { width, height } => out.push_str(&format!(
            ".rect({}, {})",
            binding(tree, width, context)?,
            binding(tree, height, context)?
        )),
        ProfileShape::Circle { radius } => {
            out.push_str(&format!(".circle({})", binding(tree, radius, context)?))
        }
        ProfileShape::Polygon { points } => {
            let pts: Vec<String> = points
                .iter()
                .map(|[x, y]| {
                    Ok(format!(
                        "[{}, {}]",
                        literal_only(tree, x, context)?,
                        literal_only(tree, y, context)?
                    ))
                })
                .collect::<Result<_, EmitError>>()?;
            out.push_str(&format!(".polygon([{}])", pts.join(", ")));
        }
        ProfileShape::Path { start, segments } => {
            out.push_str(&format!(
                ".path({}, {})",
                binding(tree, &start[0], context)?,
                binding(tree, &start[1], context)?
            ));
            for segment in segments {
                match segment {
                    PathSegment::Line { to } => out.push_str(&format!(
                        ".line_to({}, {})",
                        binding(tree, &to[0], context)?,
                        binding(tree, &to[1], context)?
                    )),
                    PathSegment::Arc { via, to } => out.push_str(&format!(
                        ".arc_to([{}, {}], [{}, {}])",
                        literal_only(tree, &via[0], context)?,
                        literal_only(tree, &via[1], context)?,
                        literal_only(tree, &to[0], context)?,
                        literal_only(tree, &to[1], context)?
                    )),
                }
            }
        }
    }
    let centred = profile.center.iter().all(|c| match c {
        ScalarBinding::Literal(v) => *v == 0.0,
        ScalarBinding::Parameter(id) => tree
            .parameters()
            .get(*id)
            .is_some_and(|p| p.expression.is_none() && p.value == 0.0),
    });
    if !centred {
        out.push_str(&format!(
            ".at({}, {})",
            binding(tree, &profile.center[0], context)?,
            binding(tree, &profile.center[1], context)?
        ));
    }
    Ok(out)
}
