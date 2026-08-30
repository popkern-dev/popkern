//! Sketch dialects for [`Feature::Sketch`](crate::Feature::Sketch)
//! (scripted-modeling design D2, PK-364).
//!
//! A sketch feature's payload is one of two arms:
//!
//! - **Constraint** — the `popkern_constraint::Sketch` JSON the tree has
//!   always carried, re-solved on replay.
//! - **Analytic** — an [`AnalyticProfile`] whose scalars are
//!   [`ScalarBinding`]s (literal numbers or parameter references),
//!   resolved against the parameter table and compiled on replay to
//!   exact geometry by `popkern-sketch` — never through the solver
//!   (design invariant 6), never converted into a constraint sketch
//!   (D2).
//!
//! **Wire form.** [`SketchPayload`] is `#[serde(untagged)]`: the
//! constraint arm serialises as the bare JSON string v3 documents
//! carry in `sketch_data`, so every existing document loads unchanged
//! and no envelope bump is needed; the analytic arm is the profile
//! object. [`ScalarBinding`] is untagged too: a bare number is a
//! literal, a parameter id is the `{idx, version}` key object.

use popkern_sketch::AnalyticProfile;
use serde::{Deserialize, Serialize};

use crate::error::ReplayError;
use crate::parameter::{ParameterId, ParameterTable};

/// A dimension of an analytic sketch: a literal value or a reference
/// into the parameter table.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum ScalarBinding {
    /// A literal value in document units.
    Literal(f64),
    /// The value of a parameter (driven by expressions and
    /// configurations like any other feature parameter).
    Parameter(ParameterId),
}

impl From<f64> for ScalarBinding {
    fn from(v: f64) -> Self {
        ScalarBinding::Literal(v)
    }
}

impl From<ParameterId> for ScalarBinding {
    fn from(id: ParameterId) -> Self {
        ScalarBinding::Parameter(id)
    }
}

impl ScalarBinding {
    /// The bound parameter, if any.
    pub fn parameter(&self) -> Option<ParameterId> {
        match self {
            ScalarBinding::Literal(_) => None,
            ScalarBinding::Parameter(id) => Some(*id),
        }
    }

    /// The scalar's value under `table`.
    pub(crate) fn resolve(&self, table: &ParameterTable) -> Result<f64, ReplayError> {
        match self {
            ScalarBinding::Literal(v) => Ok(*v),
            ScalarBinding::Parameter(id) => {
                table
                    .get(*id)
                    .map(|p| p.value)
                    .ok_or_else(|| ReplayError::FeatureFailed {
                        label: "sketch".into(),
                        error: format!("parameter {id:?} not found").into(),
                    })
            }
        }
    }
}

/// An analytic sketch as stored in a document: shape and placement
/// with every dimension a [`ScalarBinding`].
pub type AnalyticSketchDesc = AnalyticProfile<ScalarBinding>;

/// The payload of a sketch feature — one of the two sketch dialects.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum SketchPayload {
    /// JSON-encoded `popkern_constraint::Sketch` (geometry +
    /// constraints), re-solved on replay.
    Constraint(String),
    /// An analytic profile, compiled exactly on replay.
    Analytic(AnalyticSketchDesc),
}

impl From<String> for SketchPayload {
    fn from(json: String) -> Self {
        crate::SketchPayload::Constraint(json)
    }
}

impl From<&str> for SketchPayload {
    fn from(json: &str) -> Self {
        crate::SketchPayload::Constraint(json.to_string())
    }
}

impl From<AnalyticSketchDesc> for SketchPayload {
    fn from(profile: AnalyticSketchDesc) -> Self {
        SketchPayload::Analytic(profile)
    }
}

impl SketchPayload {
    /// The constraint-sketch JSON, when this is the constraint arm.
    pub fn constraint_json(&self) -> Option<&str> {
        match self {
            crate::SketchPayload::Constraint(json) => Some(json),
            SketchPayload::Analytic(_) => None,
        }
    }

    /// The analytic profile, when this is the analytic arm.
    pub fn analytic(&self) -> Option<&AnalyticSketchDesc> {
        match self {
            crate::SketchPayload::Constraint(_) => None,
            SketchPayload::Analytic(profile) => Some(profile),
        }
    }

    /// Every parameter the payload binds, with a field label for
    /// dependency reporting (the constraint arm binds none).
    pub(crate) fn param_ids(&self) -> Vec<(String, ParameterId)> {
        match self {
            crate::SketchPayload::Constraint(_) => vec![],
            SketchPayload::Analytic(profile) => profile
                .scalars()
                .into_iter()
                .zip(profile.scalar_roles())
                .filter_map(|(s, role)| s.parameter().map(|id| (role, id)))
                .collect(),
        }
    }
}

/// Resolve every binding of an analytic sketch against `table`.
pub(crate) fn resolve_profile(
    profile: &AnalyticSketchDesc,
    table: &ParameterTable,
) -> Result<AnalyticProfile<f64>, ReplayError> {
    profile.map_scalars(|s| s.resolve(table))
}

#[cfg(test)]
#[path = "feature_sketch_tests.rs"]
mod tests;
