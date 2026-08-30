//! Typed edits for `modify` (spec §7): the same parameter shapes as
//! the op that recorded the feature, applied to the document without
//! re-recording.

use popkern_history::feature::Feature;
use popkern_history::feature_blend::ContinuityKind;

use crate::scalar::Scalar;
use crate::selection::{EdgeSelector, FaceSelector};

/// An edit to apply to a recorded feature.
#[derive(Debug, Clone, Default)]
pub struct Edit {
    pub(crate) scalars: Vec<(String, Scalar)>,
    pub(crate) continuity: Option<ContinuityKind>,
    pub(crate) edges: Option<EdgeSelector>,
    pub(crate) faces: Option<FaceSelector>,
}

impl Edit {
    /// Start an empty edit.
    pub fn new() -> Edit {
        Edit::default()
    }
    /// Set a scalar role (`radius`, `distance`, `thickness`, `width`, …).
    pub fn scalar(mut self, role: &str, value: impl Into<Scalar>) -> Edit {
        self.scalars.push((role.to_string(), value.into()));
        self
    }
    /// Fillet radius.
    pub fn radius(self, value: impl Into<Scalar>) -> Edit {
        self.scalar("radius", value)
    }
    /// Chamfer/offset/extrude distance.
    pub fn distance(self, value: impl Into<Scalar>) -> Edit {
        self.scalar("distance", value)
    }
    /// Shell/thicken thickness.
    pub fn thickness(self, value: impl Into<Scalar>) -> Edit {
        self.scalar("thickness", value)
    }
    /// Fillet continuity.
    pub fn continuity(mut self, kind: ContinuityKind) -> Edit {
        self.continuity = Some(kind);
        self
    }
    /// Re-select the feature's edges (resolved against the body the
    /// feature sees, i.e. the state just before it).
    pub fn edges(mut self, selector: EdgeSelector) -> Edit {
        self.edges = Some(selector);
        self
    }
    /// Re-select the feature's faces.
    pub fn faces(mut self, selector: FaceSelector) -> Edit {
        self.faces = Some(selector);
        self
    }
    /// The roles this edit can apply to `feature`, for argument errors.
    pub(crate) fn applies_to(&self, feature: &Feature) -> Result<(), String> {
        if self.continuity.is_some() && !matches!(feature, Feature::Fillet { .. }) {
            return Err("continuity applies to fillets only".into());
        }
        if self.edges.is_some()
            && !matches!(feature, Feature::Fillet { .. } | Feature::Chamfer { .. })
        {
            return Err("edges applies to fillets and chamfers".into());
        }
        if self.faces.is_some() && !matches!(feature, Feature::Shell { .. } | Feature::Draft { .. })
        {
            return Err("faces applies to shells and drafts".into());
        }
        Ok(())
    }
}
