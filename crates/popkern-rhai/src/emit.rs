//! The canonical emitter (design D5/D7 corollary, spec §12 rule): a
//! document shown back as Rhai — only §§2–9 spellings, timeline order,
//! one statement per feature, every scalar as its current parameter
//! binding. Running the emitted script re-records an identical
//! document, and emitting that document yields the same text (fixed
//! point). Anything the skin cannot spell is a `NotApplicable` error
//! naming the feature, never a silent omission.

use std::collections::BTreeSet;

use popkern_core::tolerance::Tolerance;
use popkern_history::config::ParamOverride;
use popkern_history::feature::{Feature, OperationMode};
use popkern_history::node::BodyId;
use popkern_history::parameter::ParamScope;
use popkern_script::Model;

use crate::args::unit_name;
use crate::emit_feature::statement;
use crate::emit_spelling::{float, named_chain, quoted, selection};
use popkern_core::naming::EntityName;
use popkern_script::{parse_edges, parse_faces};

/// Why a document could not be shown as code.
#[derive(Debug, Clone, thiserror::Error)]
pub enum EmitError {
    /// A feature, parameter, or configuration has no spelling.
    #[error("cannot emit '{feature}' ({operation}): {detail}")]
    NotApplicable {
        /// The feature label (or parameter/config name).
        feature: String,
        /// Its operation identifier.
        operation: String,
        /// What is missing.
        detail: String,
    },
}

impl EmitError {
    pub(crate) fn not_applicable(
        feature: &str,
        operation: &str,
        detail: impl Into<String>,
    ) -> EmitError {
        EmitError::NotApplicable {
            feature: feature.to_string(),
            operation: operation.to_string(),
            detail: detail.into(),
        }
    }
}

fn var(body: BodyId) -> String {
    format!("b{}", body.0)
}

struct Emitter {
    out: String,
    declared: BTreeSet<u32>,
    active: Option<BodyId>,
}

impl Emitter {
    fn line(&mut self, text: &str) {
        self.out.push_str(text);
        self.out.push('\n');
    }

    fn header(&mut self, model: &Model) {
        let tolerance = model.tolerance();
        if tolerance.linear() == Tolerance::default().linear() {
            self.line("let m = model(#{ unit: \"m\" });");
        } else {
            self.line(&format!(
                "let m = model(#{{ unit: \"m\", tolerance: {} }});",
                float(tolerance.linear())
            ));
        }
    }

    fn globals(&mut self, model: &Model) -> Result<(), EmitError> {
        let table = model.tree().parameters();
        let mut literals: Vec<(String, f64, &'static str)> = Vec::new();
        let mut driven: Vec<(String, String, &'static str)> = Vec::new();
        for (id, param) in table.all_params() {
            if param.scope != ParamScope::Global || param.read_only {
                continue;
            }
            let Some(name) = table.qualified_name(id) else {
                continue;
            };
            if param.bounds.is_some() {
                return Err(EmitError::not_applicable(
                    &name,
                    "param",
                    "parameter bounds have no spelling in the skin",
                ));
            }
            let unit = unit_name(param.unit);
            match &param.expression {
                Some(expression) => driven.push((
                    name,
                    popkern_history::format_expression(expression, table),
                    unit,
                )),
                None => literals.push((name, param.unit.from_si(param.value), unit)),
            }
        }
        literals.sort_by(|a, b| a.0.cmp(&b.0));
        for (name, value, unit) in literals {
            self.line(&format!(
                "m.param({}, {}, {});",
                quoted(&name),
                float(value),
                quoted(unit)
            ));
        }
        for (name, text, unit) in driven {
            self.line(&format!(
                "m.param_expr({}, {}, {});",
                quoted(&name),
                quoted(&text),
                quoted(unit)
            ));
        }
        Ok(())
    }

    /// The body that owns the face a sketch is placed on, at the
    /// point of the timeline where the sketch is recorded: the body
    /// of the face's operation, followed through the booleans that
    /// consumed it before `index`.
    fn owner_before(
        model: &Model,
        operation: popkern_core::naming::OperationId,
        index: usize,
    ) -> Option<BodyId> {
        let tree = model.tree();
        let mut owner = tree
            .timeline()
            .iter()
            .filter_map(|&id| tree.feature(id))
            .find(|n| n.operation_id == operation)
            .map(|n| n.body)?;
        loop {
            let consumed = tree
                .timeline()
                .iter()
                .take(index)
                .filter_map(|&id| tree.feature(id))
                .find(|n| matches!(n.feature, Feature::Boolean { tool_body, .. } if tool_body == owner))
                .map(|n| n.body);
            match consumed {
                Some(host) if host != owner => owner = host,
                _ => return Some(owner),
            }
        }
    }

    fn features(&mut self, model: &mut Model) -> Result<(), EmitError> {
        let timeline: Vec<_> = model.tree().timeline().to_vec();
        for (index, id) in timeline.iter().enumerate() {
            let Some(node) = model.tree().feature(*id) else {
                continue;
            };
            let needs_active = match &node.feature {
                Feature::Sketch {
                    plane:
                        popkern_history::feature::SketchPlaneDesc::OnFace(face)
                        | popkern_history::feature::SketchPlaneDesc::OnFaceOffset(face, _),
                    ..
                } => Self::owner_before(model, face.operation, index),
                _ => None,
            };
            if let Some(owner) = needs_active
                && self.active != Some(owner)
            {
                if !self.declared.contains(&owner.0) {
                    return Err(EmitError::not_applicable(
                        &node.label,
                        node.feature.operation_label(),
                        "the sketch's face belongs to a body the script has not created",
                    ));
                }
                self.line(&format!("m.set_active({});", var(owner)));
                self.active = Some(owner);
            }
            let (label, op, body, has_selections) = (
                node.label.clone(),
                node.feature.operation_label(),
                node.body,
                !node.selections.is_empty(),
            );
            let before = if has_selections {
                model
                    .body_before(index, body)
                    .map_err(|e| EmitError::not_applicable(&label, op, e.to_string()))?
            } else {
                None
            };
            let model: &Model = &*model;
            let tree = model.tree();
            let node = tree.feature(*id).ok_or_else(|| {
                EmitError::not_applicable("timeline", "feature", "feature vanished during emission")
            })?;
            let spell = |role: &str, names: &[EntityName], faces: bool| {
                verified_selection(model, node, before.as_ref(), role, names, faces)
            };
            let stmt = statement(model, node, &var, &spell)?;
            let served = match &node.feature {
                Feature::Extrude {
                    mode: OperationMode::Boss | OperationMode::Cut,
                    ..
                } => tree.merge_host(*id).unwrap_or(node.body),
                _ => node.body,
            };
            let target = var(served);
            let text = if stmt.starts_body {
                format!("let {target} = {};", stmt.call)
            } else {
                if !self.declared.contains(&node.body.0) {
                    return Err(EmitError::not_applicable(
                        &node.label,
                        node.feature.operation_label(),
                        "the feature's body was not started by a sketch the skin can spell",
                    ));
                }
                let keyword = if self.declared.contains(&served.0) {
                    ""
                } else {
                    "let "
                };
                format!("{keyword}{target} = {}{};", var(node.body), stmt.call)
            };
            self.line(&text);
            self.declared.insert(served.0);
            self.active = Some(served);
            if node.suppressed {
                self.line(&format!("m.feature({}).suppress();", quoted(&node.label)));
            }
        }
        Ok(())
    }

    fn labels(&mut self, model: &Model) -> Result<(), EmitError> {
        let tree = model.tree();
        for (label, body) in &tree.labels().bodies {
            if !self.declared.contains(&body.0) {
                return Err(EmitError::not_applicable(
                    label,
                    "label_body",
                    "labelled body was never started in the timeline",
                ));
            }
            self.line(&format!("m.label_body({}, {});", quoted(label), var(*body)));
        }
        for (label, entity) in &tree.labels().entities {
            let spelled = model
                .spell_entity(entity)
                .map_err(|e| EmitError::not_applicable(label, "label", e.to_string()))?;
            self.line(&format!(
                "m.label_entity({}, {});",
                quoted(label),
                quoted(&spelled)
            ));
        }
        Ok(())
    }

    fn configs(&mut self, model: &Model) -> Result<(), EmitError> {
        let tree = model.tree();
        let table = tree.parameters();
        for cid in tree.configs() {
            let Some(config) = tree.config(cid) else {
                continue;
            };
            if config.parent.is_some() {
                return Err(EmitError::not_applicable(
                    &config.name,
                    "config",
                    "derived configurations have no spelling in the skin",
                ));
            }
            self.line(&format!(
                "m.add_config({}, {});",
                quoted(&config.name),
                quoted(&config.description)
            ));
            let mut overrides: Vec<(String, String)> = Vec::new();
            for (pid, value) in &config.overrides {
                let Some(name) = table.qualified_name(*pid) else {
                    continue;
                };
                let unit = table.get(*pid).map(|p| p.unit).unwrap_or_default_unit();
                match value {
                    ParamOverride::Literal(si) => {
                        overrides.push((name, float(unit.from_si(*si))));
                    }
                    ParamOverride::Expression(_) => {
                        return Err(EmitError::not_applicable(
                            &config.name,
                            "config",
                            format!("expression override of '{name}' has no spelling in the skin"),
                        ));
                    }
                }
            }
            overrides.sort();
            for (name, value) in overrides {
                self.line(&format!(
                    "m.set_override({}, {}, {});",
                    quoted(&config.name),
                    quoted(&name),
                    value
                ));
            }
            let mut suppressions: Vec<(String, bool)> = tree
                .config_suppressions(cid)
                .map(|s| {
                    s.iter()
                        .filter_map(|(fid, on)| tree.feature(*fid).map(|n| (n.label.clone(), *on)))
                        .collect()
                })
                .unwrap_or_default();
            suppressions.sort();
            for (label, on) in suppressions {
                self.line(&format!(
                    "m.set_suppression({}, m.feature({}), {on});",
                    quoted(&config.name),
                    quoted(&label)
                ));
            }
        }
        if let Some(active) = tree.active_config().and_then(|id| tree.config(id)) {
            self.line(&format!("m.set_active_config({});", quoted(&active.name)));
        }
        Ok(())
    }
}

trait UnitOrDefault {
    fn unwrap_or_default_unit(self) -> popkern_history::unit::Unit;
}

impl UnitOrDefault for Option<popkern_history::unit::Unit> {
    fn unwrap_or_default_unit(self) -> popkern_history::unit::Unit {
        self.unwrap_or(popkern_history::unit::Unit::Meters)
    }
}

/// The recorded query for `role` when it still reproduces the name
/// snapshot on the body the feature saw; otherwise the snapshot
/// itself, spelled name by name (decision D3: names are the truth,
/// the query text is diagnostics).
fn verified_selection(
    model: &Model,
    node: &popkern_history::node::FeatureNode,
    before: Option<&popkern_core::topology::BRep>,
    role: &str,
    names: &[EntityName],
    faces: bool,
) -> Result<String, EmitError> {
    let context = node.label.as_str();
    let record = node.selections.iter().find(|r| r.role == role);
    if let (Some(record), Some(brep)) = (record, before) {
        let resolved: Option<Vec<EntityName>> = if faces {
            parse_faces(&record.describe)
                .ok()
                .and_then(|sel| model.resolve_faces_in(brep, &sel).ok())
                .map(|s| s.names)
        } else {
            parse_edges(&record.describe)
                .ok()
                .and_then(|sel| model.resolve_edges_in(brep, &sel).ok())
                .map(|s| s.names)
        };
        if let Some(mut resolved) = resolved {
            let mut expected = names.to_vec();
            let key = |n: &EntityName| format!("{n:?}");
            resolved.sort_by_key(key);
            expected.sort_by_key(key);
            if resolved == expected {
                return selection(&record.describe, context);
            }
        }
    }
    named_chain(model, names, faces, context)
}

/// Show a document as the canonical Rhai script. Takes the model
/// mutably because verifying each recorded query replays the timeline
/// up to that feature.
pub fn to_rhai(model: &mut Model) -> Result<String, EmitError> {
    let mut emitter = Emitter {
        out: String::new(),
        declared: BTreeSet::new(),
        active: None,
    };
    emitter.header(model);
    emitter.globals(model)?;
    emitter.features(model)?;
    emitter.labels(model)?;
    emitter.configs(model)?;
    emitter.line("m");
    Ok(emitter.out)
}
