//! Health check for the feature tree.
//!
//! Validates that all [`EntityName`] references in the tree can be resolved
//! against the current BRep. This is a "dry run" that checks references
//! without actually executing features.

use crate::FeatureTree;
use crate::node::FeatureId;
use crate::timeline::feature_entity_names;
use popkern_core::naming::EntityName;
use popkern_core::naming_resolve::{Resolution, resolve_name};

/// A single reference validation issue.
#[derive(Debug)]
pub struct ReferenceIssue {
    /// The feature that has the broken reference.
    pub feature_id: FeatureId,
    /// The feature's label.
    pub feature_label: String,
    /// The [`EntityName`] that could not be resolved.
    pub unresolved_name: EntityName,
    /// Human-readable description of the problem.
    pub message: String,
}

/// Result of validating all references in the tree.
#[derive(Debug)]
pub struct HealthReport {
    /// Features with broken references.
    pub issues: Vec<ReferenceIssue>,
    /// Total number of references checked.
    pub total_references: usize,
    /// Number of references that resolved successfully.
    pub resolved_references: usize,
}

impl HealthReport {
    /// Whether all references are healthy.
    pub fn is_healthy(&self) -> bool {
        self.issues.is_empty()
    }
}

impl FeatureTree {
    /// Validate that all [`EntityName`] references in the tree can be resolved
    /// against the current BRep.
    ///
    /// This is a "dry run" that checks references without actually executing
    /// features. Call after replay to verify the tree is healthy, or before
    /// replay to predict failures.
    ///
    /// Returns `None` if the tree has no current BRep (call replay first).
    pub fn validate_references(&self) -> Option<HealthReport> {
        let brep = self.current_brep.as_ref()?;

        let mut issues = Vec::new();
        let mut total = 0;
        let mut resolved = 0;

        let end = self.rollback_index.unwrap_or(self.timeline.len());

        for i in 0..end {
            let id = self.timeline[i];
            let node = &self.features[id];
            if node.suppressed {
                continue;
            }

            let names = feature_entity_names(&node.feature);
            for name in names {
                total += 1;
                match resolve_name(name, brep) {
                    Resolution::Exact(_) | Resolution::Split(_) => {
                        resolved += 1;
                    }
                    Resolution::NotFound => {
                        issues.push(ReferenceIssue {
                            feature_id: id,
                            feature_label: node.label.clone(),
                            unresolved_name: name.clone(),
                            message: format!("Reference in '{}' cannot be resolved", node.label),
                        });
                    }
                }
            }
        }

        Some(HealthReport {
            issues,
            total_references: total,
            resolved_references: resolved,
        })
    }
}

#[cfg(test)]
#[path = "health_tests.rs"]
mod tests;
