//! Migration verbs (spec §7): the canonical code form of every
//! history-tree edit. Each verb is transactional — the tree is
//! snapshotted, edited, replayed, and restored if the replay fails.

use popkern_history::feature::Feature;
use popkern_history::labels::qualified_prefix;
use popkern_history::node::FeatureId;
use popkern_history::replay_report::{ReplayOutcome, ReplayReport};
use popkern_history::{DeleteBehavior, FeatureFailure, FeatureTree};

use crate::body::Body;
use crate::error::{ScriptError, SelectionKind};
use crate::feature_edit::Edit;
use crate::model::{FeatureHandle, Model};
use crate::selection::{EntityRef, resolve_edges, resolve_faces};

/// A parameter as the surface lists it.
#[derive(Debug, Clone, PartialEq)]
pub struct ParamInfo {
    /// Qualified name (`wall`, `fillet-3.radius`).
    pub name: String,
    /// Current SI value.
    pub value: f64,
    /// Display unit.
    pub unit: popkern_history::unit::Unit,
    /// Expression text, when driven.
    pub expression: Option<String>,
    /// Declared bounds.
    pub bounds: Option<(f64, f64)>,
    /// Read-only (published by a `reference`).
    pub read_only: bool,
}

impl Model {
    fn index_of(&self, feature: &FeatureHandle) -> Result<usize, ScriptError> {
        self.tree
            .timeline()
            .iter()
            .position(|&id| id == feature.id)
            .ok_or_else(|| popkern_history::HistoryError::FeatureNotFound.into())
    }

    /// Run `edit` on the tree transactionally: replay must succeed with
    /// no failed feature, else the tree is restored and the error raised.
    fn transact<T>(
        &mut self,
        edit: impl FnOnce(&mut Model) -> Result<T, ScriptError>,
    ) -> Result<T, ScriptError> {
        let snapshot: FeatureTree = self.tree.clone();
        let active = self.active;
        let outcome = edit(self).and_then(|value| {
            let report = self.tree.replay(self.tolerance)?.report;
            match report.failures().next() {
                Some(row) => Err(ScriptError::FeatureFailed {
                    label: row.label.clone(),
                    operation: row.operation.clone(),
                    error: Box::new(match &row.outcome {
                        ReplayOutcome::Failed { error } => error.clone(),
                        other => FeatureFailure::other(format!("{other:?}")),
                    }),
                    report: Box::new(report.clone()),
                }),
                None => Ok(value),
            }
        });
        if outcome.is_err() {
            self.tree = snapshot;
            self.active = active;
            let _ = self.tree.replay(self.tolerance);
        }
        self.note_full_replay();
        outcome
    }

    /// Plain-language caption of a feature (spec §7 `describe`).
    pub fn describe(&self, feature: &FeatureHandle) -> Result<String, ScriptError> {
        let node = self
            .tree
            .feature(feature.id)
            .ok_or(popkern_history::HistoryError::FeatureNotFound)?;
        let prefix = qualified_prefix(&node.label);
        let params = self.tree.parameters();
        let scalar = |role: &str| {
            params
                .resolve_name(&format!("{prefix}.{role}"))
                .and_then(|id| params.get(id))
                .map(|p| p.value)
        };
        let selection = |role: &str| {
            node.selections
                .iter()
                .find(|s| s.role == role)
                .map(|s| s.describe.clone())
        };
        let text = match &node.feature {
            Feature::Sketch { .. } => "sketch".to_string(),
            Feature::Extrude { .. } => match scalar("distance") {
                Some(d) => format!("extrude {d}"),
                None => "extrude".to_string(),
            },
            Feature::Fillet { edges, recipe, .. } if !recipe.edge_radii.is_empty() => format!(
                "fillet {} radius groups ({}) on {} edges{}",
                recipe.edge_radii.len() + 1,
                std::iter::once(scalar("radius"))
                    .chain((1..=recipe.edge_radii.len()).map(|k| scalar(&format!("radius_{k}"))))
                    .map(|r| r.map_or_else(|| "?".to_string(), |v| v.to_string()))
                    .collect::<Vec<_>>()
                    .join(", "),
                edges.len(),
                match recipe.continuity {
                    popkern_history::feature_blend::ContinuityKind::G2 => " (G2)",
                    _ => "",
                },
            ),
            Feature::Fillet { edges, recipe, .. } => format!(
                "fillet {} on {} edges{}: {}",
                scalar("radius").unwrap_or(f64::NAN),
                edges.len(),
                match recipe.continuity {
                    popkern_history::feature_blend::ContinuityKind::G2 => " (G2)",
                    _ => "",
                },
                selection("edges").unwrap_or_else(|| "(unrecorded selection)".into())
            ),
            Feature::Chamfer { edges, .. } => format!(
                "chamfer {} on {} edges: {}",
                scalar("distance").unwrap_or(f64::NAN),
                edges.len(),
                selection("edges").unwrap_or_else(|| "(unrecorded selection)".into())
            ),
            Feature::Shell {
                faces_to_remove, ..
            } => format!(
                "shell {} removing {} faces: {}",
                scalar("thickness").unwrap_or(f64::NAN),
                faces_to_remove.len(),
                selection("faces_to_remove").unwrap_or_else(|| "(unrecorded selection)".into())
            ),
            Feature::Boolean { operation, .. } => format!("boolean {operation:?}"),
            other => other.operation_label().to_string(),
        };
        Ok(format!("{} — {text}", node.label))
    }

    /// Apply a typed edit to a feature and replay.
    pub fn modify(&mut self, feature: &FeatureHandle, edit: Edit) -> Result<(), ScriptError> {
        let id = feature.id;
        let node = self
            .tree
            .feature(id)
            .ok_or(popkern_history::HistoryError::FeatureNotFound)?
            .clone();
        edit.applies_to(&node.feature)
            .map_err(|reason| ScriptError::Argument {
                argument: "edit".into(),
                reason,
            })?;
        let index = self.index_of(feature)?;
        self.transact(|m| {
            let prefix = qualified_prefix(&node.label);
            for (role, value) in &edit.scalars {
                let name = format!("{prefix}.{role}");
                let pid = m.tree.parameters().resolve_name(&name).ok_or_else(|| {
                    ScriptError::Argument {
                        argument: role.clone(),
                        reason: format!("feature '{}' has no parameter '{role}'", node.label),
                    }
                })?;
                match value {
                    crate::scalar::Scalar::Literal(v) => {
                        let unit = m
                            .tree
                            .parameters()
                            .get(pid)
                            .map(|p| p.unit)
                            .unwrap_or(m.unit);
                        m.tree.set_parameter_value_display(pid, *v, unit)?;
                    }
                    crate::scalar::Scalar::Param(target) => {
                        let target_id =
                            m.tree.parameters().resolve_name(target).ok_or_else(|| {
                                ScriptError::Argument {
                                    argument: role.clone(),
                                    reason: format!("unknown parameter '{target}'"),
                                }
                            })?;
                        m.tree.set_parameter_expression(
                            pid,
                            popkern_history::Expression::Param(target_id),
                        )?;
                    }
                    crate::scalar::Scalar::Expr(text) => {
                        let expression =
                            popkern_history::parse::parse_expression(text, m.tree.parameters())
                                .map_err(|e| ScriptError::Argument {
                                    argument: role.clone(),
                                    reason: format!("expression '{text}': {e}"),
                                })?;
                        m.tree.set_parameter_expression(pid, expression)?;
                    }
                }
            }
            let mut new_feature = node.feature.clone();
            let mut records = node.selections.clone();
            if edit.continuity.is_some() || edit.edges.is_some() || edit.faces.is_some() {
                // Selections resolve against the body the feature runs
                // on, as it stands just before the feature.
                let before = m
                    .tree
                    .state_at(index, m.tolerance)?
                    .body(node.body)
                    .cloned()
                    .ok_or_else(|| ScriptError::Argument {
                        argument: "edit".into(),
                        reason: format!(
                            "feature '{}' has no body to select on before it runs",
                            node.label
                        ),
                    })?;
                let labels = |r: EntityRef| m.resolve_entity_ref(r);
                match (&mut new_feature, &edit.edges, &edit.faces) {
                    (Feature::Fillet { edges, recipe, .. }, sel, _) => {
                        if let Some(kind) = edit.continuity {
                            recipe.continuity = kind;
                        }
                        if sel.is_some() && !recipe.edge_radii.is_empty() {
                            return Err(ScriptError::Argument {
                                argument: "edges".into(),
                                reason: "a fillet with radius groups selects its edges per group: \
                                         re-record it with fillet_each"
                                    .into(),
                            });
                        }
                        if let Some(sel) = sel {
                            let snap = resolve_edges(&before, sel, &labels)?;
                            *edges = snap.names.clone();
                            records.retain(|r| r.role != "edges");
                            records.push(popkern_history::SelectionRecord {
                                role: "edges".into(),
                                describe: snap.describe,
                                labels_used: snap.labels_used,
                                seams_skipped: snap.seams_skipped,
                            });
                        }
                    }
                    (Feature::Chamfer { edges, .. }, Some(sel), _) => {
                        let snap = resolve_edges(&before, sel, &labels)?;
                        *edges = snap.names.clone();
                        records.retain(|r| r.role != "edges");
                        records.push(popkern_history::SelectionRecord {
                            role: "edges".into(),
                            describe: snap.describe,
                            labels_used: snap.labels_used,
                            seams_skipped: snap.seams_skipped,
                        });
                    }
                    (
                        Feature::Shell {
                            faces_to_remove, ..
                        },
                        _,
                        Some(sel),
                    ) => {
                        let snap = resolve_faces(&before, sel, &labels)?;
                        *faces_to_remove = snap.names.clone();
                        records.retain(|r| r.role != "faces_to_remove");
                        records.push(popkern_history::SelectionRecord {
                            role: "faces_to_remove".into(),
                            describe: snap.describe,
                            labels_used: snap.labels_used,
                            seams_skipped: 0,
                        });
                    }
                    _ => {}
                }
                m.tree.modify_feature(id, new_feature);
                let _ = m.tree.set_selection_records(id, records);
            }
            Ok(())
        })
    }

    /// Record one feature immediately before `anchor`.
    pub fn insert_before(
        &mut self,
        anchor: &FeatureHandle,
        build: impl FnOnce(&mut Model, Body) -> Result<Body, ScriptError>,
    ) -> Result<FeatureHandle, ScriptError> {
        let index = self.index_of(anchor)?;
        self.insert_at(index, build)
    }

    /// Record one feature immediately after `anchor`.
    pub fn insert_after(
        &mut self,
        anchor: &FeatureHandle,
        build: impl FnOnce(&mut Model, Body) -> Result<Body, ScriptError>,
    ) -> Result<FeatureHandle, ScriptError> {
        let index = self.index_of(anchor)? + 1;
        self.insert_at(index, build)
    }

    fn insert_at(
        &mut self,
        index: usize,
        build: impl FnOnce(&mut Model, Body) -> Result<Body, ScriptError>,
    ) -> Result<FeatureHandle, ScriptError> {
        let before = self.tree.timeline().len();
        self.transact(|m| {
            m.tree.rollback_to(index);
            m.tree.replay(m.tolerance)?;
            let body = m
                .tree
                .timeline()
                .get(index.saturating_sub(1))
                .and_then(|&id| m.tree.feature(id))
                .map(|n| Body::new(n.body))
                .ok_or(ScriptError::NoActiveBody)?;
            let outcome = build(m, body);
            m.tree.rollback_to_end();
            outcome?;
            if m.tree.timeline().len() != before + 1 {
                return Err(ScriptError::Argument {
                    argument: "insert".into(),
                    reason: "the closure must record exactly one feature".into(),
                });
            }
            let id = m.tree.timeline()[index];
            m.handle(id)
        })
    }

    /// Move a feature to `index` in the timeline.
    pub fn reorder(&mut self, feature: &FeatureHandle, index: usize) -> Result<(), ScriptError> {
        let id = feature.id;
        self.transact(|m| Ok(m.tree.reorder(id, index)?))
    }

    /// Suppress a feature.
    pub fn suppress(&mut self, feature: &FeatureHandle) -> Result<(), ScriptError> {
        let id = feature.id;
        self.transact(|m| {
            m.tree.suppress(id);
            Ok(())
        })
    }

    /// Unsuppress a feature.
    pub fn unsuppress(&mut self, feature: &FeatureHandle) -> Result<(), ScriptError> {
        let id = feature.id;
        self.transact(|m| {
            m.tree.unsuppress(id);
            Ok(())
        })
    }

    /// Delete a feature and everything that depends on it.
    pub fn delete(&mut self, feature: &FeatureHandle) -> Result<Vec<FeatureId>, ScriptError> {
        let id = feature.id;
        self.transact(|m| Ok(m.tree.delete_feature(id, DeleteBehavior::Cascade)?.deleted))
    }

    /// Replay up to (not including) `feature`; later recording inserts there.
    pub fn rollback_to(&mut self, feature: &FeatureHandle) -> Result<ReplayReport, ScriptError> {
        let index = self.index_of(feature)?;
        self.tree.rollback_to(index);
        let report = self.tree.replay(self.tolerance)?.report;
        self.note_full_replay();
        Ok(report)
    }

    /// Clear the rollback bar.
    pub fn rollback_to_end(&mut self) -> Result<ReplayReport, ScriptError> {
        self.tree.rollback_to_end();
        let report = self.tree.replay(self.tolerance)?.report;
        self.note_full_replay();
        Ok(report)
    }

    /// Features that transitively depend on `feature`, in timeline order.
    pub fn dependents(&self, feature: &FeatureHandle) -> Vec<FeatureHandle> {
        self.tree
            .dependents(feature.id)
            .into_iter()
            .filter_map(|id| self.handle(id).ok())
            .collect()
    }

    /// Undo the last document edit; `false` when nothing to undo.
    pub fn undo(&mut self) -> Result<bool, ScriptError> {
        let done = self.tree.undo();
        if done {
            self.tree.replay(self.tolerance)?;
            self.note_full_replay();
        }
        Ok(done)
    }

    /// Redo the last undone edit.
    pub fn redo(&mut self) -> Result<bool, ScriptError> {
        let done = self.tree.redo();
        if done {
            self.tree.replay(self.tolerance)?;
            self.note_full_replay();
        }
        Ok(done)
    }

    /// Every parameter as the surface lists it.
    pub fn params(&self) -> Vec<ParamInfo> {
        let table = self.tree.parameters();
        table
            .all_params()
            .filter_map(|(id, p)| {
                Some(ParamInfo {
                    name: table.qualified_name(id)?,
                    value: p.value,
                    unit: p.unit,
                    expression: p
                        .expression
                        .as_ref()
                        .map(|e| popkern_history::format_expression(e, table)),
                    bounds: p.bounds,
                    read_only: p.read_only,
                })
            })
            .collect()
    }

    /// Selection kinds a feature's arguments use (for skins' editors).
    pub fn selection_kinds(&self, feature: &FeatureHandle) -> Vec<(String, SelectionKind)> {
        self.tree
            .feature(feature.id)
            .map(|n| {
                n.selections
                    .iter()
                    .map(|s| {
                        let kind = if s.role.contains("face") {
                            SelectionKind::Faces
                        } else {
                            SelectionKind::Edges
                        };
                        (s.role.clone(), kind)
                    })
                    .collect()
            })
            .unwrap_or_default()
    }
}
