//! Convenience API for creating pattern features with ergonomic spacing.
//!
//! This module adds builder methods to [`FeatureTree`] for the four pattern
//! feature variants (linear, circular, mirror) without requiring callers to
//! manually pre-allocate parameter IDs or construct the inner types.

use std::collections::HashSet;

use slotmap::SlotMap;

use crate::FeatureTree;
use crate::error::HistoryError;
use crate::expression::{BinOp, Expression};
use crate::feature::{AxisDesc, Feature, Vec3Desc};
use crate::node::{FeatureId, FeatureNode};
use crate::parameter::ParameterId;
use crate::pattern::{InstanceOverride, MirrorPlaneDesc, PatternDirection};
use crate::pattern_validate::validate_seed_features;
use crate::unit::Unit;

/// Collect the union of all seed [`FeatureId`]s that are already claimed
/// by an existing pattern feature in the tree.
///
/// Used by the convenience API to satisfy the "no shared seeds" rule checked
/// by [`validate_seed_features`].
pub(crate) fn collect_existing_seeds(
    features: &SlotMap<FeatureId, FeatureNode>,
) -> HashSet<FeatureId> {
    features
        .values()
        .filter_map(|node| match &node.feature {
            Feature::LinearPattern { seed_features, .. }
            | Feature::CircularPattern { seed_features, .. }
            | Feature::CurvePattern { seed_features, .. }
            | Feature::Mirror { seed_features, .. } => Some(seed_features.iter().copied()),
            _ => None,
        })
        .flatten()
        .collect()
}

impl FeatureTree {
    /// Add a linear pattern using pitch-based spacing.
    ///
    /// Creates three feature-scoped parameters:
    /// - `{label}.primary_count` — number of instances (including the seed).
    /// - `{label}.primary_pitch` — distance between consecutive instances.
    /// - `{label}.primary_total_distance` — total span; driven by the
    ///   expression `pitch * (count - 1)` so that changing count or pitch
    ///   automatically re-evaluates the total distance.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::DuplicateParameterName`] if any generated
    /// parameter name conflicts with an existing entry in the parameter table.
    pub fn add_linear_pattern_with_pitch(
        &mut self,
        seed_features: Vec<FeatureId>,
        direction: Vec3Desc,
        count: f64,
        pitch: f64,
        label: &str,
    ) -> Result<FeatureId, HistoryError> {
        let existing_seeds = collect_existing_seeds(&self.features);
        validate_seed_features(
            &seed_features,
            &self.features,
            &self.timeline,
            &existing_seeds,
        )?;

        let prefix = label.to_lowercase().replace(' ', "_");

        let count_id = self.parameter_table.add_global(
            &format!("{prefix}.primary_count"),
            count,
            Unit::Ratio,
        )?;

        // pitch is stored in SI (metres); the caller passes display-unit mm.
        let pitch_si = Unit::Millimeters.to_si(pitch);
        let pitch_id = self.parameter_table.add_global(
            &format!("{prefix}.primary_pitch"),
            pitch_si,
            Unit::Millimeters,
        )?;

        // total_distance = pitch * (count - 1), stored as an expression so that
        // changing count or pitch automatically re-evaluates the total distance.
        let total_si = pitch_si * (count - 1.0).max(0.0);
        let total_id = self.parameter_table.add_global(
            &format!("{prefix}.primary_total_distance"),
            total_si,
            Unit::Millimeters,
        )?;
        // pitch * (count - 1)
        let expr = Expression::BinOp(
            Box::new(Expression::Param(pitch_id)),
            BinOp::Mul,
            Box::new(Expression::BinOp(
                Box::new(Expression::Param(count_id)),
                BinOp::Sub,
                Box::new(Expression::Literal(1.0)),
            )),
        );
        self.parameter_table.set_expression(total_id, expr)?;

        let feature = crate::feature::Feature::LinearPattern {
            seed_features,
            primary: PatternDirection {
                direction,
                count: count_id,
                total_distance: total_id,
            },
            secondary: None,
            geometry_pattern: false,
            instance_overrides: std::collections::HashMap::new(),
        };

        Ok(self.add_feature(feature, label))
    }

    /// Add a circular pattern that fills a full 360° turn (2π radians).
    ///
    /// Creates two feature-scoped parameters:
    /// - `{label}.count` — number of instances (including the seed).
    /// - `{label}.total_angle` — set to τ (2π radians = full circle).
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::DuplicateParameterName`] if any generated
    /// parameter name conflicts with an existing entry in the parameter table.
    pub fn add_circular_pattern_full_circle(
        &mut self,
        seed_features: Vec<FeatureId>,
        axis: AxisDesc,
        count: f64,
        label: &str,
    ) -> Result<FeatureId, HistoryError> {
        let existing_seeds = collect_existing_seeds(&self.features);
        validate_seed_features(
            &seed_features,
            &self.features,
            &self.timeline,
            &existing_seeds,
        )?;

        let prefix = label.to_lowercase().replace(' ', "_");

        let count_id =
            self.parameter_table
                .add_global(&format!("{prefix}.count"), count, Unit::Ratio)?;

        let total_angle_id = self.parameter_table.add_global(
            &format!("{prefix}.total_angle"),
            std::f64::consts::TAU,
            Unit::Radians,
        )?;

        let feature = crate::feature::Feature::CircularPattern {
            seed_features,
            axis,
            count: count_id,
            total_angle: total_angle_id,
            geometry_pattern: false,
            instance_overrides: std::collections::HashMap::new(),
        };

        Ok(self.add_feature(feature, label))
    }

    /// Add a circular pattern using an angular pitch between instances.
    ///
    /// Creates three feature-scoped parameters:
    /// - `{label}.count` — number of instances (including the seed).
    /// - `{label}.pitch_angle` — angular spacing between consecutive instances
    ///   (in radians).
    /// - `{label}.total_angle` — total angular span; driven by the expression
    ///   `pitch_angle * (count - 1)` so that changing count or pitch
    ///   automatically re-evaluates the total angle.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::DuplicateParameterName`] if any generated
    /// parameter name conflicts with an existing entry in the parameter table.
    pub fn add_circular_pattern_with_pitch(
        &mut self,
        seed_features: Vec<FeatureId>,
        axis: AxisDesc,
        count: f64,
        pitch_angle: f64,
        label: &str,
    ) -> Result<FeatureId, HistoryError> {
        let existing_seeds = collect_existing_seeds(&self.features);
        validate_seed_features(
            &seed_features,
            &self.features,
            &self.timeline,
            &existing_seeds,
        )?;

        let prefix = label.to_lowercase().replace(' ', "_");

        let count_id =
            self.parameter_table
                .add_global(&format!("{prefix}.count"), count, Unit::Ratio)?;

        // pitch_angle is in radians (caller passes radians directly).
        let pitch_id = self.parameter_table.add_global(
            &format!("{prefix}.pitch_angle"),
            pitch_angle,
            Unit::Radians,
        )?;

        // total_angle = pitch_angle * (count - 1), stored as an expression so
        // that changing count or pitch_angle automatically re-evaluates.
        let total_angle = pitch_angle * (count - 1.0).max(0.0);
        let total_angle_id = self.parameter_table.add_global(
            &format!("{prefix}.total_angle"),
            total_angle,
            Unit::Radians,
        )?;
        // pitch_angle * (count - 1)
        let expr = Expression::BinOp(
            Box::new(Expression::Param(pitch_id)),
            BinOp::Mul,
            Box::new(Expression::BinOp(
                Box::new(Expression::Param(count_id)),
                BinOp::Sub,
                Box::new(Expression::Literal(1.0)),
            )),
        );
        self.parameter_table.set_expression(total_angle_id, expr)?;

        let feature = crate::feature::Feature::CircularPattern {
            seed_features,
            axis,
            count: count_id,
            total_angle: total_angle_id,
            geometry_pattern: false,
            instance_overrides: std::collections::HashMap::new(),
        };

        Ok(self.add_feature(feature, label))
    }

    /// Add a mirror feature that reflects seed features about a plane.
    ///
    /// Validates that `seed_features` is non-empty, all features exist and
    /// belong to the same body, form a contiguous timeline range, start with a
    /// [`Feature::Sketch`], contain at least one Boss/Cut operation, and are not
    /// already owned by another pattern.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError`] if any seed-feature validation rule is violated.
    pub fn add_mirror(
        &mut self,
        seed_features: Vec<FeatureId>,
        plane: MirrorPlaneDesc,
        merge_result: bool,
        label: &str,
    ) -> Result<FeatureId, HistoryError> {
        let existing_seeds = collect_existing_seeds(&self.features);
        validate_seed_features(
            &seed_features,
            &self.features,
            &self.timeline,
            &existing_seeds,
        )?;

        let feature = crate::feature::Feature::Mirror {
            seed_features,
            plane,
            merge_result,
        };
        Ok(self.add_feature(feature, label))
    }

    /// Suppress a specific instance of a pattern feature.
    ///
    /// If `pattern_id` does not refer to a pattern feature with an
    /// `instance_overrides` map, this is a no-op.  The suppression is stored
    /// directly in the feature's data; callers must call
    /// [`modify_feature`](Self::modify_feature) (or [`replay`](Self::replay))
    /// to propagate the change to the BRep.
    pub fn skip_pattern_instance(&mut self, pattern_id: FeatureId, instance_idx: usize) {
        let Some(node) = self.features.get_mut(pattern_id) else {
            return;
        };
        match &mut node.feature {
            crate::feature::Feature::LinearPattern {
                instance_overrides, ..
            }
            | crate::feature::Feature::CircularPattern {
                instance_overrides, ..
            }
            | crate::feature::Feature::CurvePattern {
                instance_overrides, ..
            } => {
                let entry =
                    instance_overrides
                        .entry(instance_idx)
                        .or_insert_with(|| InstanceOverride {
                            suppressed: false,
                            param_overrides: std::collections::BTreeMap::new(),
                        });
                entry.suppressed = true;
            }
            _ => {}
        }
    }

    /// Override a parameter value for a specific pattern instance.
    ///
    /// The override is stored in the `instance_overrides` map of the pattern
    /// feature.  If `pattern_id` does not refer to a pattern feature with an
    /// `instance_overrides` map, this is a no-op.
    pub fn override_pattern_instance_param(
        &mut self,
        pattern_id: FeatureId,
        instance_idx: usize,
        param: ParameterId,
        value: f64,
    ) {
        let Some(node) = self.features.get_mut(pattern_id) else {
            return;
        };
        match &mut node.feature {
            crate::feature::Feature::LinearPattern {
                instance_overrides, ..
            }
            | crate::feature::Feature::CircularPattern {
                instance_overrides, ..
            }
            | crate::feature::Feature::CurvePattern {
                instance_overrides, ..
            } => {
                let entry =
                    instance_overrides
                        .entry(instance_idx)
                        .or_insert_with(|| InstanceOverride {
                            suppressed: false,
                            param_overrides: std::collections::BTreeMap::new(),
                        });
                entry.param_overrides.insert(param, value);
            }
            _ => {}
        }
    }
}

#[cfg(test)]
#[path = "pattern_api_tests.rs"]
mod tests;
