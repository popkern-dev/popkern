//! The partition's operations on a body (design: operations on the cell
//! complex, 5.4): a split keeping the pieces a rule selects, an N-ary
//! boolean by set expression, the symmetric difference, an imprint.
//! Each consumes its tools, as a boolean consumes its tool.

use popkern_history::feature::{Feature, NonManifoldPolicy, PieceRule, SetExpr};
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
