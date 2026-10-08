//! The history timeline of a session's document
//! (live-projects-and-history.md §5.6): one row per timeline feature of
//! every part, and the meshes of the bodies alive after any step — each
//! keyed by the body and the step that last changed it, so a page that
//! already holds a key never receives (or re-meshes) that body again,
//! and each carrying the faces the step itself created (its
//! highlight).

use std::collections::HashMap;

use serde::Serialize;

use popkern_core::topology::BRep;
use popkern_history::node::{BodyId, FeatureId};
use popkern_history::replay_report::ReplayOutcome;
use popkern_rhai::Outcome;
use popkern_script::Model;
use popkern_script::source::SourceSpan;
use slotmap::Key;

use crate::envelope::BodyMesh;
use crate::mesh_buffers::{MeshPurpose, body_mesh};
use crate::session::Session;

/// One timeline step.
#[derive(Debug, Clone, Serialize)]
pub struct TimelineStep {
    /// 1-based: the document after this many features.
    pub index: usize,
    /// The feature's stable id.
    pub feature: u64,
    /// Its label.
    pub label: String,
    /// The operation (`extrude`, `fillet`, …).
    pub operation: String,
    /// The project migration that recorded it, in a project.
    pub migration: Option<String>,
    /// Whether the last replay skipped it as suppressed.
    pub suppressed: bool,
    /// Whether the last replay recorded it failed.
    pub failed: bool,
    /// Bodies this step built, modified or merged into.
    pub bodies_changed: Vec<u64>,
    /// Where in a script it was recorded, when tracked.
    pub source: Option<SourceSpan>,
}

/// The timeline of one part.
#[derive(Debug, Clone, Serialize)]
pub struct PartTimeline {
    /// The part definition (assemblies); `None` for a part document.
    pub part: Option<String>,
    /// Steps in timeline order.
    pub steps: Vec<TimelineStep>,
}

/// One body after a step.
#[derive(Debug, Clone, Serialize)]
pub struct StepBody {
    /// `body:feature` — the body and the step that last changed it.
    /// Equal keys are identical meshes.
    pub key: String,
    /// The body's tree id.
    pub body: u64,
    /// Indices into the mesh's `face_labels` of the faces this step
    /// created (empty when the step left the body untouched).
    pub created_faces: Vec<u32>,
    /// Whether the session's step buffers hold this body's mesh (false
    /// when the page said it already has the key).
    pub buffers: bool,
    /// The mesh description (buffers travel separately).
    pub mesh: Option<BodyMesh>,
}

/// The session's cache of step meshes, keyed by part, body and the
/// step that last changed the body.
pub(crate) type StepMeshCache = HashMap<(Option<String>, u32, FeatureId), BodyMesh>;

fn feature_key(id: FeatureId) -> u64 {
    id.data().as_ffi()
}

fn part_timeline(
    part: Option<String>,
    model: &mut Model,
    migrations: &HashMap<(Option<String>, FeatureId), String>,
) -> PartTimeline {
    let ids: Vec<FeatureId> = model.tree().timeline().to_vec();
    let mut steps = Vec::with_capacity(ids.len());
    for (i, id) in ids.iter().enumerate() {
        let Some(node) = model.tree().feature(*id) else {
            continue;
        };
        let row = model
            .tree()
            .last_report()
            .features
            .iter()
            .find(|r| r.feature == *id);
        let bodies_changed = model
            .tree()
            .cached_state_at(i + 1)
            .map(|state| {
                state
                    .body_ids()
                    .into_iter()
                    .filter(|b| state.changed_by(*b) == Some(*id))
                    .map(|b| u64::from(b.0))
                    .collect()
            })
            .unwrap_or_default();
        steps.push(TimelineStep {
            index: i + 1,
            feature: feature_key(*id),
            label: node.label.clone(),
            operation: node.feature.operation_label().to_string(),
            migration: migrations.get(&(part.clone(), *id)).cloned(),
            suppressed: row.is_some_and(|r| matches!(r.outcome, ReplayOutcome::Suppressed { .. })),
            failed: row.is_some_and(|r| matches!(r.outcome, ReplayOutcome::Failed { .. })),
            bodies_changed,
            source: model.source_of(*id).cloned(),
        });
    }
    PartTimeline { part, steps }
}

/// The faces of `brep` named under one of `operations`.
fn created_faces(brep: &BRep, operations: &[u64]) -> Vec<u32> {
    brep.faces
        .iter()
        .enumerate()
        .filter(|(_, (_, face))| {
            face.name
                .and_then(|n| brep.name_pool.get(n))
                .is_some_and(|name| operations.contains(&name.operation.0))
        })
        .map(|(i, _)| i as u32)
        .collect()
}

impl Session {
    /// The timeline of every part of the document.
    pub fn timeline(&mut self) -> Result<Vec<PartTimeline>, Outcome> {
        let migrations = self.migrations.clone();
        let document = self
            .document
            .as_mut()
            .ok_or_else(crate::session::no_model)?;
        Ok(document
            .models_mut()
            .into_iter()
            .map(|(part, model)| part_timeline(part, model, &migrations))
            .collect())
    }

    /// The bodies alive after the first `index` features of `part`
    /// (`None` for a part document). Bodies whose key is in `known`
    /// carry no mesh (the page holds it); the others are meshed once
    /// per key and kept.
    pub fn step_bodies(
        &mut self,
        part: Option<&str>,
        index: usize,
        known: &[String],
    ) -> Result<Vec<StepBody>, Outcome> {
        let document = self
            .document
            .as_mut()
            .ok_or_else(crate::session::no_model)?;
        let gate = |detail: String| Outcome::GateFailure {
            step: "step_meshes".into(),
            detail,
        };
        let (part, model) = document
            .models_mut()
            .into_iter()
            .find(|(p, _)| p.as_deref() == part)
            .ok_or_else(|| gate(format!("no part {part:?} in the document")))?;
        let state = model.step_state(index).map_err(|e| gate(e.to_string()))?;
        let step = index
            .checked_sub(1)
            .and_then(|i| model.tree().timeline().get(i).copied());
        let operations: Vec<u64> = step
            .and_then(|id| model.tree().feature(id))
            .map(|node| {
                let mut ops = vec![node.operation_id.0];
                ops.extend(node.instance_op_ids.values().map(|op| op.0));
                ops
            })
            .unwrap_or_default();
        let mut out = Vec::new();
        for (body, brep) in state.bodies() {
            let Some(changed_by) = state.changed_by(body) else {
                continue;
            };
            let key = format!("{}:{}", body.0, feature_key(changed_by));
            let created = if Some(changed_by) == step {
                created_faces(brep, &operations)
            } else {
                Vec::new()
            };
            if known.contains(&key) {
                out.push(StepBody {
                    key,
                    body: u64::from(body.0),
                    created_faces: created,
                    buffers: false,
                    mesh: None,
                });
                continue;
            }
            let cache_key = (part.clone(), body.0, changed_by);
            if !self.step_meshes.contains_key(&cache_key) {
                let mesh = mesh_step_body(model, body, brep).map_err(gate)?;
                self.step_meshes.insert(cache_key.clone(), mesh);
            }
            out.push(StepBody {
                key,
                body: u64::from(body.0),
                created_faces: created,
                buffers: true,
                mesh: self.step_meshes.get(&cache_key).cloned(),
            });
        }
        Ok(out)
    }
}

/// The display mesh of `brep` (body `body` of `model`).
fn mesh_step_body(model: &Model, body: BodyId, brep: &BRep) -> Result<BodyMesh, String> {
    let label_face = |handle| {
        let name = brep.faces.get(handle)?.name?;
        let entity = brep.name_pool.get(name)?;
        model.spell_entity(entity).ok()
    };
    body_mesh(
        brep,
        u64::from(body.0),
        model.tree().body_label(body).map(str::to_string),
        label_face,
        MeshPurpose::Display,
    )
    .map(|mut mesh| {
        mesh.appearance = model.appearance();
        mesh
    })
    .map_err(|e| e.to_string())
}
