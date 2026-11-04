//! The partition's operations on a body (design: operations on the cell
//! complex, 5.4): a split keeping the pieces a rule selects, an N-ary
//! boolean by set expression, the symmetric difference, an imprint,
//! and taking a body's pieces apart. Each consumes its tools, as a
//! boolean consumes its tool.

use popkern_history::feature::{BoolOpKind, Feature, NonManifoldPolicy, PieceRule, SetExpr};
use popkern_history::node::BodyId;

use crate::body::Body;
use crate::error::ScriptError;
use crate::model::Model;

/// One of the partition's operations, with its tool bodies.
#[derive(Debug, Clone)]
pub enum PartitionOp {
    /// Cut by `tools`, keeping the pieces `keep` selects.
    Split {
        /// The cutting bodies (tool k is `InsideTools` body k).
        tools: Vec<Body>,
        /// Which pieces the body keeps.
        keep: PieceRule,
    },
    /// The regions of the body (body 0) and `tools` (bodies 1..N) that
    /// `expr` keeps.
    Nary {
        /// The tool bodies.
        tools: Vec<Body>,
        /// The regions kept.
        expr: SetExpr,
        /// Kept regions touching only along edges or at points.
        non_manifold: NonManifoldPolicy,
    },
    /// The body's faces split where `tool` crosses them.
    Imprint {
        /// The imprinting body.
        tool: Body,
    },
    /// Keep the body's pieces (solids) `rule` selects, removing the
    /// rest (a geometric rule: `All`, `Containing`, `Largest`).
    KeepPieces {
        /// Which pieces stay.
        rule: PieceRule,
    },
    /// Move the body's pieces `rule` selects into a new body (the
    /// result); the body keeps the rest.
    SeparatePieces {
        /// Which pieces move.
        rule: PieceRule,
    },
}

/// A piece rule a body's own pieces can answer: tool membership exists
/// only on the split that made them.
fn geometric(rule: PieceRule) -> Result<PieceRule, ScriptError> {
    match rule {
        PieceRule::InsideTools(_) => Err(ScriptError::Argument {
            argument: "rule".into(),
            reason: "a body's pieces have no tools to be inside: select by tool on the split, \
                     or by `containing` / `largest` here"
                .into(),
        }),
        rule => Ok(rule),
    }
}

impl Body {
    /// Check `tools` are alive, distinct and not this body; their ids.
    fn tool_ids(self, model: &mut Model, tools: &[Body]) -> Result<Vec<BodyId>, ScriptError> {
        self.require_alive(model)?;
        if tools.is_empty() {
            return Err(ScriptError::Argument {
                argument: "tools".into(),
                reason: "at least one tool body is needed".into(),
            });
        }
        let mut ids: Vec<BodyId> = Vec::with_capacity(tools.len());
        for tool in tools {
            tool.require_alive(model)?;
            if tool.id == self.id {
                return Err(ScriptError::Argument {
                    argument: "tools".into(),
                    reason: "a body cannot be its own tool".into(),
                });
            }
            if ids.contains(&tool.id) {
                return Err(ScriptError::Argument {
                    argument: "tools".into(),
                    reason: "a tool body is given twice".into(),
                });
            }
            ids.push(tool.id);
        }
        Ok(ids)
    }

    /// Record `op` on this body under `label` (empty: generated). Its
    /// tools are consumed.
    ///
    /// # Errors
    ///
    /// [`ScriptError::Argument`] for no tools, a repeated tool, the body
    /// itself, or a set expression naming a body past the tools; a
    /// replay failure (a split rule selecting no piece included).
    pub fn partition_with(
        self,
        model: &mut Model,
        op: PartitionOp,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let feature = match op {
            PartitionOp::Split { tools, keep } => Feature::Split {
                tools: self.tool_ids(model, &tools)?,
                keep,
            },
            PartitionOp::Nary {
                tools,
                expr,
                non_manifold,
            } => {
                let tools = self.tool_ids(model, &tools)?;
                if let Some(i) = expr.max_body()
                    && i > tools.len()
                {
                    return Err(ScriptError::Argument {
                        argument: "expr".into(),
                        reason: format!(
                            "the expression names body {i}, but there are {} bodies",
                            tools.len() + 1
                        ),
                    });
                }
                Feature::NaryBoolean {
                    tools,
                    expr,
                    non_manifold,
                }
            }
            PartitionOp::Imprint { tool } => Feature::Imprint {
                tool: self.tool_ids(model, &[tool])?[0],
            },
            PartitionOp::KeepPieces { rule } => {
                self.require_alive(model)?;
                Feature::KeepPieces {
                    rule: geometric(rule)?,
                }
            }
            PartitionOp::SeparatePieces { rule } => {
                self.require_alive(model)?;
                Feature::SeparatePieces {
                    source: self.id,
                    rule: geometric(rule)?,
                }
            }
        };
        self.record(model, feature, label, Vec::new())
    }

    /// This body cut into its pieces by `tools` (consumed), keeping the
    /// pieces `keep` selects as one body of several solids.
    ///
    /// # Errors
    ///
    /// As [`Body::partition_with`].
    pub fn split(
        self,
        model: &mut Model,
        tools: &[Body],
        keep: PieceRule,
    ) -> Result<Body, ScriptError> {
        let tools = tools.to_vec();
        self.partition_with(model, PartitionOp::Split { tools, keep }, "")
    }

    /// The regions of this body (body 0) and `tools` (bodies 1..N,
    /// consumed) that `expr` keeps, as one body; kept regions touching
    /// only along edges or at points follow `non_manifold`.
    ///
    /// # Errors
    ///
    /// As [`Body::partition_with`].
    pub fn nary_boolean(
        self,
        model: &mut Model,
        tools: &[Body],
        expr: SetExpr,
        non_manifold: NonManifoldPolicy,
    ) -> Result<Body, ScriptError> {
        let op = PartitionOp::Nary {
            tools: tools.to_vec(),
            expr,
            non_manifold,
        };
        self.partition_with(model, op, "")
    }

    /// The material in exactly one of this body and `tool` (consumed).
    /// Where the two boundaries cross, the two sides' pieces meet along
    /// the crossing curves with zero thickness: they are served as
    /// manifold pieces with recorded contacts.
    ///
    /// # Errors
    ///
    /// As [`Body::partition_with`].
    pub fn symmetric_difference(self, model: &mut Model, tool: Body) -> Result<Body, ScriptError> {
        self.nary_boolean(
            model,
            &[tool],
            SetExpr::ExactlyOne,
            NonManifoldPolicy::SplitToManifold,
        )
    }

    /// Boolean `operation` with `tool` (consumed), labelled; a result
    /// touching only along edges or at points follows `non_manifold`.
    ///
    /// # Errors
    ///
    /// As [`Body::boolean_with`].
    pub fn boolean_under(
        self,
        model: &mut Model,
        operation: BoolOpKind,
        tool: Body,
        non_manifold: NonManifoldPolicy,
        label: &str,
    ) -> Result<Body, ScriptError> {
        self.boolean_labeled(model, operation, tool, non_manifold, label)
    }

    /// This body with only the pieces (solids) `rule` selects.
    ///
    /// # Errors
    ///
    /// As [`Body::partition_with`]; `InsideTools` is an argument error.
    pub fn keep_pieces(self, model: &mut Model, rule: PieceRule) -> Result<Body, ScriptError> {
        self.partition_with(model, PartitionOp::KeepPieces { rule }, "")
    }

    /// A new body of this body's pieces `rule` selects; this body keeps
    /// the rest. Returns the new body.
    ///
    /// # Errors
    ///
    /// As [`Body::keep_pieces`]; a rule selecting every piece fails
    /// (this body would be empty).
    pub fn separate_pieces(self, model: &mut Model, rule: PieceRule) -> Result<Body, ScriptError> {
        self.partition_with(model, PartitionOp::SeparatePieces { rule }, "")
    }

    /// This body with its faces split where `tool` (consumed) crosses
    /// them.
    ///
    /// # Errors
    ///
    /// As [`Body::partition_with`].
    pub fn imprint(self, model: &mut Model, tool: Body) -> Result<Body, ScriptError> {
        self.partition_with(model, PartitionOp::Imprint { tool }, "")
    }
}
