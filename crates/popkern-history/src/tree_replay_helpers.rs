//! Helper functions and methods used by the tree replay loop.
//!
//! Extracted from `tree_replay.rs` to keep it within the 500-line structural
//! limit.

use std::collections::HashMap;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::FeatureTree;
use crate::config::ConfigId;
use crate::config_eval;
use crate::error::ReplayError;
use crate::feature::Feature;
use crate::node::{BodyId, FeatureId};
use crate::parameter::ParameterId;
use crate::unit::{DimensionKind, Unit};

// -- Reference parameter helpers ----------------------------------------------

/// Map a [`DimensionKind`] to its default SI [`Unit`].
pub(super) fn dimension_kind_to_unit(kind: DimensionKind) -> Unit {
    match kind {
        DimensionKind::Length => Unit::Meters,
        DimensionKind::Area => Unit::SquareMeters,
        DimensionKind::Volume => Unit::CubicMeters,
        DimensionKind::Angle => Unit::Radians,
        DimensionKind::InverseLength => Unit::PerMeter,
        DimensionKind::MomentOfArea => Unit::MetersToFourth,
        DimensionKind::Dimensionless => Unit::Ratio,
    }
}

impl FeatureTree {
    /// Create or update the read-only parameter for a Reference feature.
    ///
    /// If the parameter already exists (by qualified name), its value and
    /// `last_known` are updated in place.  Otherwise a new feature-scoped
    /// parameter is created.
    pub(super) fn upsert_reference_param(
        &mut self,
        feature_id: FeatureId,
        qualified_name: &str,
        value: f64,
        unit: Unit,
    ) {
        if let Some(pid) = self.parameter_table.resolve_name(qualified_name) {
            if let Some(param) = self.parameter_table.get_mut(pid) {
                param.value = value;
                param.last_known = Some(value);
            }
        } else if let Ok(pid) =
            self.parameter_table
                .add_feature_param(qualified_name, feature_id, value, unit)
            && let Some(param) = self.parameter_table.get_mut(pid)
        {
            param.read_only = true;
            param.last_known = Some(value);
        }
    }
}

// -- Pattern execution helper -------------------------------------------------

/// Return `true` if the feature is a pattern variant.
pub(super) fn is_pattern_feature(feature: &Feature) -> bool {
    matches!(
        feature,
        Feature::LinearPattern { .. }
            | Feature::CircularPattern { .. }
            | Feature::CurvePattern { .. }
            | Feature::Mirror { .. }
    )
}

/// Extract seed feature IDs from a pattern [`Feature`] variant.
fn extract_seed_ids(feature: &Feature) -> Vec<FeatureId> {
    match feature {
        Feature::LinearPattern { seed_features, .. }
        | Feature::CircularPattern { seed_features, .. }
        | Feature::CurvePattern { seed_features, .. }
        | Feature::Mirror { seed_features, .. } => seed_features.clone(),
        _ => Vec::new(),
    }
}

impl FeatureTree {
    /// Extract data from the tree and delegate to [`crate::pattern_execute`].
    pub(super) fn execute_pattern_feature(
        &mut self,
        id: FeatureId,
        body_breps: &mut HashMap<BodyId, BRep>,
        tolerance: &Tolerance,
    ) -> Result<(), ReplayError> {
        let node = &self.features[id];
        let feat = node.feature.clone();
        let op_id = node.operation_id;
        let body = node.body;
        let seed_ids = extract_seed_ids(&feat);
        let seed_feats: Vec<_> = seed_ids
            .iter()
            .filter_map(|&s| self.features.get(s).map(|n| n.feature.clone()))
            .collect();
        let seed_ops: Vec<_> = seed_ids
            .iter()
            .filter_map(|&s| self.features.get(s).map(|n| n.operation_id))
            .collect();
        let mut inst_ops = self.features[id].instance_op_ids.clone();
        let mut ctx = crate::pattern_execute::PatternContext {
            pattern_feature: &feat,
            seed_features: &seed_feats,
            seed_op_ids: &seed_ops,
            instance_op_ids: &mut inst_ops,
            body_breps,
            param_table: &self.parameter_table,
            tolerance,
            next_op_id: &mut self.next_op_id,
            pattern_op_id: op_id,
            seed_body: body,
            features: &self.features,
            seed_feature_ids: &seed_ids,
        };
        let r = crate::pattern_execute::execute_pattern(&mut ctx);
        self.features[id].instance_op_ids = inst_ops;
        r.map(|_| ())
    }
}

// -- Config-aware parameter accessors -----------------------------------------

impl FeatureTree {
    /// Effective parameter value for the active config.
    pub fn effective_value(&self, param: ParameterId) -> f64 {
        self.parameter_table
            .get(param)
            .map(|p| p.value)
            .unwrap_or(0.0)
    }

    /// Effective suppression state for a feature considering the active config.
    pub fn effective_suppressed(&self, feature: FeatureId) -> bool {
        if let Some(overridden) =
            config_eval::resolve_suppression(feature, self.active_config, &self.configs)
        {
            return overridden;
        }
        self.features
            .get(feature)
            .map(|n| n.suppressed)
            .unwrap_or(false)
    }

    /// Evaluate all parameters for a specific config (expensive -- full evaluation).
    pub fn evaluate_config(&self, config: Option<ConfigId>) -> HashMap<ParameterId, f64> {
        config_eval::evaluate_config(config, &self.configs, &self.parameter_table)
    }
}
