//! Additional convenience API for creating pattern features.
//!
//! This module extends [`FeatureTree`] with pattern creation methods that
//! accept different parameterisation styles (total-distance, angular span,
//! curve-driven) and grid-index helper utilities.

use crate::FeatureTree;
use crate::error::HistoryError;
use crate::feature::{AxisDesc, Feature, Vec3Desc};
use crate::node::FeatureId;
use crate::pattern::PatternDirection;
use crate::pattern_api::collect_existing_seeds;
use crate::pattern_validate::validate_seed_features;
use crate::unit::Unit;
use popkern_core::naming::EntityName;

impl FeatureTree {
    /// Add a linear pattern using total-distance spacing.
    ///
    /// Creates two feature-scoped parameters as plain values (no expression
    /// linking between them):
    /// - `{label}.primary_count` — number of instances (including the seed).
    /// - `{label}.primary_total_distance` — total span from seed to last
    ///   instance.
    ///
    /// Use this variant when the total distance is the primary design intent
    /// and the pitch is derived (or irrelevant).
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::DuplicateParameterName`] if any generated
    /// parameter name conflicts with an existing entry in the parameter table.
    pub fn add_linear_pattern_with_total(
        &mut self,
        seed_features: Vec<FeatureId>,
        direction: Vec3Desc,
        count: f64,
        total_distance: f64,
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

        // total_distance is stored in SI (metres); the caller passes mm.
        let total_si = Unit::Millimeters.to_si(total_distance);
        let total_id = self.parameter_table.add_global(
            &format!("{prefix}.primary_total_distance"),
            total_si,
            Unit::Millimeters,
        )?;

        let feature = Feature::LinearPattern {
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

    /// Add a circular pattern using an explicit total angular span.
    ///
    /// Creates two feature-scoped parameters as plain values (no expression
    /// linking between them):
    /// - `{label}.count` — number of instances (including the seed).
    /// - `{label}.total_angle` — total angular span in radians.
    ///
    /// Use this variant when the total span is the primary design intent and
    /// the angular pitch is derived (or irrelevant).
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::DuplicateParameterName`] if any generated
    /// parameter name conflicts with an existing entry in the parameter table.
    pub fn add_circular_pattern_with_span(
        &mut self,
        seed_features: Vec<FeatureId>,
        axis: AxisDesc,
        count: f64,
        total_angle: f64,
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

        // total_angle is in radians (caller passes radians directly).
        let total_angle_id = self.parameter_table.add_global(
            &format!("{prefix}.total_angle"),
            total_angle,
            Unit::Radians,
        )?;

        let feature = Feature::CircularPattern {
            seed_features,
            axis,
            count: count_id,
            total_angle: total_angle_id,
            geometry_pattern: false,
            instance_overrides: std::collections::HashMap::new(),
        };

        Ok(self.add_feature(feature, label))
    }

    /// Add a curve-driven pattern that places instances along an edge curve.
    ///
    /// Creates one feature-scoped parameter:
    /// - `{label}.count` — number of instances (including the seed).
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::DuplicateParameterName`] if any generated
    /// parameter name conflicts with an existing entry in the parameter table.
    pub fn add_curve_pattern(
        &mut self,
        seed_features: Vec<FeatureId>,
        curve: EntityName,
        count: f64,
        equal_spacing: bool,
        align_to_curve: bool,
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

        let feature = Feature::CurvePattern {
            seed_features,
            curve,
            count: count_id,
            equal_spacing,
            align_to_curve,
            geometry_pattern: false,
            instance_overrides: std::collections::HashMap::new(),
        };

        Ok(self.add_feature(feature, label))
    }

    /// Suppress a specific instance in a 2D grid pattern using row/column
    /// coordinates.
    ///
    /// The flat index is computed as `row * count_secondary + col` (row-major),
    /// where `count_secondary` is read from the secondary direction's count
    /// parameter in the parameter table.
    ///
    /// If `pattern_id` does not refer to a [`Feature::LinearPattern`] with a
    /// secondary direction, or if the secondary count parameter is missing or
    /// has a non-positive value, this is a no-op.
    pub fn skip_pattern_instance_2d(&mut self, pattern_id: FeatureId, row: usize, col: usize) {
        // Look up the secondary count to compute the flat index.
        let secondary_count_id = match self.features.get(pattern_id).map(|n| &n.feature) {
            Some(Feature::LinearPattern {
                secondary: Some(sec),
                ..
            }) => sec.count,
            _ => return,
        };

        let count2 = match self.parameter_table.get(secondary_count_id) {
            Some(p) if p.value > 0.0 => p.value as usize,
            _ => return,
        };

        let flat = row * count2 + col;
        self.skip_pattern_instance(pattern_id, flat);
    }
}

#[cfg(test)]
#[path = "pattern_api_extras_tests.rs"]
mod tests;
