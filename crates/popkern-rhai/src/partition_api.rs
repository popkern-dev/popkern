//! The partition's operations on a script body (design: operations on
//! the cell complex, 5.4): `b.split([t1, t2], keep?)`,
//! `b.nary_boolean([t1, t2], expr, options?)`,
//! `b.symmetric_difference(t)`, `b.imprint(t)`, `b.keep_pieces(rule)`,
//! `b.separate_pieces(rule)` (returns the new body), with set expressions
//! built by `body(0).and(body(1))`, `at_least(2)`, `exactly_one()`, and
//! piece rules by `all_pieces()`, `inside_tools(expr)`,
//! `containing([x, y, z])`, `largest(n)`.

use popkern_history::feature::{NonManifoldPolicy, PieceRule, SetExpr};
use popkern_script::Body;
use popkern_script::prelude::PartitionOp;
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args::{self, Options};
use crate::body_api::ScriptBody;
use crate::error::{IntoScriptResult, argument_error};

const NARY_OPTIONS: &[&str] = &["non_manifold", "label"];
const LABEL_ONLY: &[&str] = &["label"];

/// The tool bodies of a `[t1, t2, ..]` array.
fn tools(argument: &str, value: rhai::Array) -> Result<Vec<Body>, Box<EvalAltResult>> {
    value
        .into_iter()
        .enumerate()
        .map(|(i, v)| {
            v.try_cast::<ScriptBody>()
                .map(|b| b.body)
                .ok_or_else(|| argument_error(argument, format!("item {i} is not a body")))
        })
        .collect()
}

/// A body index, a non-negative integer.
fn index(argument: &str, value: rhai::INT) -> Result<usize, Box<EvalAltResult>> {
    usize::try_from(value)
        .map_err(|_| argument_error(argument, format!("expected a count, got {value}")))
}

/// A non-manifold policy spelled `decline` or `split_to_manifold`.
pub(crate) fn non_manifold(
    argument: &str,
    value: &Dynamic,
) -> Result<NonManifoldPolicy, Box<EvalAltResult>> {
    args::spelling(
        argument,
        value,
        &[
            ("decline", NonManifoldPolicy::Decline),
            ("split_to_manifold", NonManifoldPolicy::SplitToManifold),
            ("general", NonManifoldPolicy::General),
        ],
    )
}

/// The script spelling of a non-default policy (`None` for the
/// default, which the emitted call leaves out).
pub(crate) fn non_manifold_spelling(policy: NonManifoldPolicy) -> Option<&'static str> {
    match policy {
        NonManifoldPolicy::Decline => None,
        NonManifoldPolicy::SplitToManifold => Some("split_to_manifold"),
        NonManifoldPolicy::General => Some("general"),
    }
}

impl ScriptBody {
    /// Record `op` with the options map's label.
    fn partition(
        &self,
        op: PartitionOp,
        options: &Options,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let body = self
            .body
            .partition_with(&mut self.model.borrow_mut(), op, &options.label()?)
            .script()?;
        Ok(ScriptBody::new(self.model.clone(), body))
    }

    fn split_by(
        &self,
        tools: Vec<Body>,
        keep: PieceRule,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let options = Options::parse("split", options, LABEL_ONLY)?;
        self.partition(PartitionOp::Split { tools, keep }, &options)
    }

    fn nary_by(
        &self,
        tools: Vec<Body>,
        expr: SetExpr,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let options = Options::parse("nary_boolean", options, NARY_OPTIONS)?;
        let non_manifold = options
            .with("non_manifold", non_manifold)?
            .unwrap_or_default();
        let op = PartitionOp::Nary {
            tools,
            expr,
            non_manifold,
        };
        self.partition(op, &options)
    }

    fn symmetric_difference_with(
        &self,
        tool: ScriptBody,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let options = Options::parse("symmetric_difference", options, LABEL_ONLY)?;
        let op = PartitionOp::Nary {
            tools: vec![tool.body],
            expr: SetExpr::ExactlyOne,
            non_manifold: NonManifoldPolicy::SplitToManifold,
        };
        self.partition(op, &options)
    }

    fn pieces_by(
        &self,
        name: &'static str,
        op: PartitionOp,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let options = Options::parse(name, options, LABEL_ONLY)?;
        self.partition(op, &options)
    }

    fn imprint_with(
        &self,
        tool: ScriptBody,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let options = Options::parse("imprint", options, LABEL_ONLY)?;
        self.partition(PartitionOp::Imprint { tool: tool.body }, &options)
    }
}

/// Register the set expressions, piece rules and partition operations.
pub fn register(engine: &mut Engine) {
    engine
        .register_type_with_name::<SetExpr>("SetExpr")
        .register_type_with_name::<PieceRule>("PieceRule")
        .register_fn("body", |i: rhai::INT| index("body", i).map(SetExpr::Body))
        .register_fn("any_body", || SetExpr::Any)
        .register_fn("all_bodies", || SetExpr::All)
        .register_fn("exactly_one", || SetExpr::ExactlyOne)
        .register_fn("odd", || SetExpr::Odd)
        .register_fn("at_least", |k: rhai::INT| {
            index("at_least", k).map(SetExpr::AtLeast)
        })
        .register_fn("and", |a: &mut SetExpr, b: SetExpr| {
            SetExpr::And(Box::new(a.clone()), Box::new(b))
        })
        .register_fn("or", |a: &mut SetExpr, b: SetExpr| {
            SetExpr::Or(Box::new(a.clone()), Box::new(b))
        })
        .register_fn("minus", |a: &mut SetExpr, b: SetExpr| {
            SetExpr::Minus(Box::new(a.clone()), Box::new(b))
        })
        .register_fn("all_pieces", || PieceRule::All)
        .register_fn("inside_tools", PieceRule::InsideTools)
        .register_fn("containing", |p: Dynamic| {
            args::vec3("point", &p).map(PieceRule::Containing)
        })
        .register_fn("largest", |n: rhai::INT| {
            index("largest", n).map(PieceRule::Largest)
        })
        .register_fn(
            "split",
            |b: &mut ScriptBody, t: rhai::Array| -> Result<ScriptBody, Box<EvalAltResult>> {
                b.split_by(tools("tools", t)?, PieceRule::All, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "split",
            |b: &mut ScriptBody,
             t: rhai::Array,
             keep: PieceRule|
             -> Result<ScriptBody, Box<EvalAltResult>> {
                b.split_by(tools("tools", t)?, keep, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "split",
            |b: &mut ScriptBody,
             t: rhai::Array,
             keep: PieceRule,
             o: Dynamic|
             -> Result<ScriptBody, Box<EvalAltResult>> {
                b.split_by(tools("tools", t)?, keep, &o)
            },
        )
        .register_fn(
            "nary_boolean",
            |b: &mut ScriptBody,
             t: rhai::Array,
             expr: SetExpr|
             -> Result<ScriptBody, Box<EvalAltResult>> {
                b.nary_by(tools("tools", t)?, expr, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "nary_boolean",
            |b: &mut ScriptBody,
             t: rhai::Array,
             expr: SetExpr,
             o: Dynamic|
             -> Result<ScriptBody, Box<EvalAltResult>> {
                b.nary_by(tools("tools", t)?, expr, &o)
            },
        )
        .register_fn(
            "symmetric_difference",
            |b: &mut ScriptBody, tool: ScriptBody| {
                b.symmetric_difference_with(tool, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "symmetric_difference",
            |b: &mut ScriptBody, tool: ScriptBody, o: Dynamic| {
                b.symmetric_difference_with(tool, &o)
            },
        )
        .register_fn("imprint", |b: &mut ScriptBody, tool: ScriptBody| {
            b.imprint_with(tool, &Dynamic::UNIT)
        })
        .register_fn(
            "imprint",
            |b: &mut ScriptBody, tool: ScriptBody, o: Dynamic| b.imprint_with(tool, &o),
        )
        .register_fn("keep_pieces", |b: &mut ScriptBody, rule: PieceRule| {
            b.pieces_by(
                "keep_pieces",
                PartitionOp::KeepPieces { rule },
                &Dynamic::UNIT,
            )
        })
        .register_fn(
            "keep_pieces",
            |b: &mut ScriptBody, rule: PieceRule, o: Dynamic| {
                b.pieces_by("keep_pieces", PartitionOp::KeepPieces { rule }, &o)
            },
        )
        .register_fn("separate_pieces", |b: &mut ScriptBody, rule: PieceRule| {
            let op = PartitionOp::SeparatePieces { rule };
            b.pieces_by("separate_pieces", op, &Dynamic::UNIT)
        })
        .register_fn(
            "separate_pieces",
            |b: &mut ScriptBody, rule: PieceRule, o: Dynamic| {
                let op = PartitionOp::SeparatePieces { rule };
                b.pieces_by("separate_pieces", op, &o)
            },
        );
}
