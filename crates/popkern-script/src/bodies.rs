//! Body bookkeeping (spec §9): tool models live in the same tree under
//! their own `BodyId`; bodies are labelled, looked up, activated, and
//! checked for liveness (a boolean consumes its tool).

use popkern_core::topology::BRep;
use popkern_history::SnapshotOrigin;
use popkern_history::feature::Feature;

use crate::body::Body;
use crate::error::ScriptError;
use crate::model::Model;
use crate::record::record_feature;

impl Model {
    /// Record an embedded body snapshot (`raw_body`, spec §5) as a new
    /// body; the snapshot is validated at replay and named afresh.
    pub fn raw_body(&mut self, brep: BRep, origin: SnapshotOrigin) -> Result<Body, ScriptError> {
        let feature = Feature::RawBody {
            snapshot: Box::new(brep),
            origin,
        };
        let id = record_feature(
            &mut self.tree,
            self.tolerance,
            feature,
            "",
            None,
            Vec::new(),
        )?;
        let body = self
            .tree
            .feature(id)
            .map(|n| n.body)
            .ok_or(popkern_history::HistoryError::FeatureNotFound)?;
        self.active = Some(body);
        Ok(Body::new(body))
    }

    /// Label a body (`DuplicateLabel` if the label names another body).
    pub fn label_body(&mut self, name: &str, body: Body) -> Result<(), ScriptError> {
        Ok(self.tree.label_body(name, body.id)?)
    }

    /// Look a body up by label.
    pub fn body(&self, label: &str) -> Result<Body, ScriptError> {
        self.tree
            .body_by_label(label)
            .map(Body::new)
            .ok_or_else(|| popkern_history::HistoryError::LabelNotFound(label.to_string()).into())
    }

    /// A body's label, if any.
    pub fn body_label(&self, body: Body) -> Option<String> {
        self.tree.body_label(body.id).map(str::to_string)
    }

    /// Bodies alive at the end of the timeline (consumed tools excluded).
    pub fn alive_bodies(&mut self) -> Vec<Body> {
        if self.tree.current_brep().is_none() {
            let _ = self.tree.replay(self.tolerance);
        }
        self.tree
            .bodies()
            .into_iter()
            .filter(|&b| self.tree.body_brep(b).is_some())
            .map(Body::new)
            .collect()
    }

    /// Make `body` the one new features attach to (`BodyNotAlive` if
    /// a boolean consumed it).
    pub fn set_active(&mut self, body: Body) -> Result<(), ScriptError> {
        if !body.is_alive(self) {
            return Err(ScriptError::BodyNotAlive(body.id.0 as u64));
        }
        self.active = Some(body.id);
        Ok(())
    }
}

impl Body {
    /// Whether the body still exists at the end of the timeline.
    pub fn is_alive(&self, model: &Model) -> bool {
        model.tree.body_brep(self.id).is_some()
    }
}
