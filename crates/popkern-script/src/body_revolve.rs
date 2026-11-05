//! Revolve (spec §5): a body's last sketch turned about an axis, by an
//! angle or up to a face (of the revolved body or, brep-bench §4.4, of
//! another: a rib turned until it meets a stop), as a new body or a
//! boss into the face's body.

use popkern_history::feature::{AxisDesc, Feature, OperationMode, RevolveEnd as TreeRevolveEnd};
use popkern_history::unit::Unit;

use crate::body::{Body, RevolveEnd};
use crate::body_recipes::AxisSpec;
use crate::error::ScriptError;
use crate::model::Model;
use popkern_history::node::SelectionRecord;

impl Body {
    /// Revolve the body's last sketch about `axis`.
    pub fn revolve(
        self,
        model: &mut Model,
        axis: AxisDesc,
        end: RevolveEnd,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let face = self.last_sketch_face(model)?;
        let mut records = Vec::new();
        let end_condition = self.revolve_end(model, end, &mut records)?;
        let feature = Feature::Revolve {
            face,
            axis,
            end_condition,
            mode: OperationMode::NewBody,
        };
        self.record(model, feature, label, records)
    }

    /// Revolve the body's last sketch about a queried or literal axis,
    /// with a label.
    pub fn revolve_with(
        self,
        model: &mut Model,
        axis: AxisSpec,
        end: RevolveEnd,
        label: &str,
    ) -> Result<Body, ScriptError> {
        self.revolve_as(model, axis, end, OperationMode::NewBody, label)
    }

    /// Revolve with a mode: a boss into the body owning an up-to face
    /// (or the body the sketch sits on), or a new body.
    pub fn revolve_as(
        self,
        model: &mut Model,
        axis: AxisSpec,
        end: RevolveEnd,
        mode: OperationMode,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let face = self.last_sketch_face(model)?;
        let mut records = Vec::new();
        let axis = self.axis(model, axis, &mut records)?;
        let end_condition = self.revolve_end(model, end, &mut records)?;
        let feature = Feature::Revolve {
            face,
            axis,
            end_condition,
            mode,
        };
        self.record(model, feature, label, records)
    }

    fn revolve_end(
        self,
        model: &mut Model,
        end: RevolveEnd,
        records: &mut Vec<SelectionRecord>,
    ) -> Result<TreeRevolveEnd, ScriptError> {
        Ok(match end {
            RevolveEnd::Blind(a) => TreeRevolveEnd::Blind {
                angle: Body::bind(model, "angle", &a, Unit::Degrees)?,
            },
            RevolveEnd::FullCircle => TreeRevolveEnd::FullCircle,
            RevolveEnd::MidPlane(a) => TreeRevolveEnd::MidPlane {
                half_angle: Body::bind(model, "half_angle", &a, Unit::Degrees)?,
            },
            RevolveEnd::UpToSurface(sel) => TreeRevolveEnd::UpToSurface {
                face: self.one_face_recorded(model, "end_face", &sel, records)?,
            },
            RevolveEnd::UpToSurfaceOn(on, sel) => TreeRevolveEnd::UpToSurface {
                face: on.one_face_recorded(model, "end_face", &sel, records)?,
            },
        })
    }
}
