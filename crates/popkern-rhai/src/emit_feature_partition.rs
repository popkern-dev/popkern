//! Spellings of the partition's features (split from `emit_feature`
//! under the file budget): `.split([t..], keep)`,
//! `.nary_boolean([t..], expr, #{ non_manifold })`,
//! `.symmetric_difference(t)` and `.imprint(t)`, with set expressions
//! and piece rules spelled by their builders.

use popkern_history::feature::{Feature, NonManifoldPolicy, PieceRule, SetExpr};
use popkern_history::node::{BodyId, FeatureNode};

use crate::emit_spelling::{float, label_option, options, quoted};

/// A set expression as its builder calls.
fn set_expr(expr: &SetExpr) -> String {
    match expr {
        SetExpr::Body(i) => format!("body({i})"),
        SetExpr::Any => "any_body()".into(),
        SetExpr::All => "all_bodies()".into(),
        SetExpr::ExactlyOne => "exactly_one()".into(),
        SetExpr::Odd => "odd()".into(),
        SetExpr::AtLeast(k) => format!("at_least({k})"),
        SetExpr::And(a, b) => format!("{}.and({})", set_expr(a), set_expr(b)),
        SetExpr::Or(a, b) => format!("{}.or({})", set_expr(a), set_expr(b)),
        SetExpr::Minus(a, b) => format!("{}.minus({})", set_expr(a), set_expr(b)),
    }
}

/// A piece rule as its builder call.
fn piece_rule(rule: &PieceRule) -> String {
    match rule {
        PieceRule::All => "all_pieces()".into(),
        PieceRule::InsideTools(expr) => format!("inside_tools({})", set_expr(expr)),
        PieceRule::Containing([x, y, z]) => {
            format!("containing([{}, {}, {}])", float(*x), float(*y), float(*z))
        }
        PieceRule::Largest(n) => format!("largest({n})"),
    }
}

/// The call for a split, N-ary boolean or imprint `feature`; `None`
/// for any other feature.
pub(super) fn call(
    feature: &Feature,
    label: &str,
    tool_var: &dyn Fn(BodyId) -> String,
) -> Option<String> {
    let list = |tools: &[BodyId]| {
        let names: Vec<String> = tools.iter().map(|t| tool_var(*t)).collect();
        format!("[{}]", names.join(", "))
    };
    Some(match feature {
        Feature::Split { tools, keep } => format!(
            ".split({}, {}{})",
            list(tools),
            piece_rule(keep),
            options(vec![label_option(label)])
        ),
        Feature::NaryBoolean {
            tools,
            expr: SetExpr::ExactlyOne,
            non_manifold: NonManifoldPolicy::SplitToManifold,
        } if tools.len() == 1 => format!(
            ".symmetric_difference({}{})",
            tool_var(tools[0]),
            options(vec![label_option(label)])
        ),
        Feature::NaryBoolean {
            tools,
            expr,
            non_manifold,
        } => {
            let mut pairs = Vec::new();
            if let Some(spelling) = crate::partition_api::non_manifold_spelling(*non_manifold) {
                pairs.push(("non_manifold", quoted(spelling)));
            }
            pairs.push(label_option(label));
            format!(
                ".nary_boolean({}, {}{})",
                list(tools),
                set_expr(expr),
                options(pairs)
            )
        }
        Feature::Imprint { tool } => format!(
            ".imprint({}{})",
            tool_var(*tool),
            options(vec![label_option(label)])
        ),
        Feature::KeepPieces { rule } => format!(
            ".keep_pieces({}{})",
            piece_rule(rule),
            options(vec![label_option(label)])
        ),
        // A new body, spelled with its receiver: `let b2 = b0.separate_pieces(..)`.
        Feature::SeparatePieces { source, rule } => format!(
            "{}.separate_pieces({}{})",
            tool_var(*source),
            piece_rule(rule),
            options(vec![label_option(label)])
        ),
        _ => return None,
    })
}

/// The array of `node`'s served bodies, for a feature serving several
/// (PK-1120): a split's call returns its pieces (bound to `{body}_pieces`),
/// a boolean's are read back with `served_bodies()`.
pub(crate) fn pieces_var(node: &FeatureNode, body: &str) -> String {
    if matches!(node.feature, Feature::Split { .. }) {
        format!("{body}_pieces")
    } else {
        format!("{body}.served_bodies()")
    }
}

/// The lines binding the further bodies `node` served, each by its
/// place in the served array (the own body first, then the others in
/// the order they were recorded): `let b7 = b0_pieces[1];`.
pub(crate) fn served_body_lines(
    node: &FeatureNode,
    var: &dyn Fn(BodyId) -> String,
) -> Vec<(BodyId, String)> {
    let pieces = pieces_var(node, &var(node.body));
    node.output_bodies
        .iter()
        .map(|o| o.body)
        .filter(|&b| b != node.body)
        .enumerate()
        .map(|(k, b)| (b, format!("let {} = {pieces}[{}];", var(b), k + 1)))
        .collect()
}
