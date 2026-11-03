//! Spellings of the partition's features (split from `emit_feature`
//! under the file budget): `.split([t..], keep)`,
//! `.nary_boolean([t..], expr, #{ non_manifold })`,
//! `.symmetric_difference(t)` and `.imprint(t)`, with set expressions
//! and piece rules spelled by their builders.

use popkern_history::feature::{Feature, NonManifoldPolicy, PieceRule, SetExpr};
use popkern_history::node::BodyId;

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
            if *non_manifold == NonManifoldPolicy::SplitToManifold {
                pairs.push(("non_manifold", quoted("split_to_manifold")));
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
        _ => return None,
    })
}
