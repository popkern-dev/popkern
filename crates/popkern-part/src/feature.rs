//! Feature log: every operation applied to a [`Part`](crate::Part)
//! is recorded with the parameters that fully define it, per the
//! UI-forward "undoable by construction" rule. v1 records; replay /
//! parametric re-execution is a future milestone (see the design
//! doc's exploration section).

/// One recorded modeling operation.
#[derive(Debug, Clone)]
pub struct Feature {
    /// Short human-readable name ("extrude", "union", "fillet").
    pub name: &'static str,
    /// Parameter summary, derivable-from-parameters per UI rule 5
    /// (e.g. "radius 0.08 on 9 edges: circular ∧ between(Cylinder, Plane)").
    pub summary: String,
}

/// Append-only log of the features that produced a part.
#[derive(Debug, Clone, Default)]
pub struct FeatureLog {
    entries: Vec<Feature>,
}

impl FeatureLog {
    /// Record a feature.
    pub(crate) fn push(&mut self, name: &'static str, summary: impl Into<String>) {
        self.entries.push(Feature {
            name,
            summary: summary.into(),
        });
    }

    /// The recorded features, oldest first.
    pub fn entries(&self) -> &[Feature] {
        &self.entries
    }

    /// Merge another part's history (used by boolean combinations:
    /// the result's log is `self`'s log, then the other body's log,
    /// then the combination itself).
    pub(crate) fn extend_from(&mut self, other: &FeatureLog) {
        self.entries.extend(other.entries.iter().cloned());
    }
}
