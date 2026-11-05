//! The other body a feature names, whose state its selection is
//! checked against: an extrude's or revolve's end face, a sweep's or
//! pipe's path, a fill's wire.

use popkern_history::feature::{ExtrudeEnd, Feature};
use popkern_history::node::BodyId;
use popkern_script::Model;

use super::Emitter;

impl Emitter {
    /// The body `feature` (at timeline `index`, on `body`) names
    /// besides its own, when there is one.
    pub(super) fn other_body(
        model: &Model,
        feature: &Feature,
        index: usize,
        body: BodyId,
    ) -> Option<BodyId> {
        // An extrude up to (or offset from) another body's face:
        // that body is named in the end, and the face's selection
        // is checked against it.
        match feature {
            Feature::Extrude {
                end_condition:
                    ExtrudeEnd::UpToSurface { face } | ExtrudeEnd::OffsetFromSurface { face, .. },
                ..
            }
            | Feature::Revolve {
                end_condition: popkern_history::feature::RevolveEnd::UpToSurface { face },
                ..
            } => Self::owner_before(model, face.operation, index).filter(|o| *o != body),
            // A sweep's or pipe's path edge lives on the path body
            // (a wire, or a solid's edge), named in the call.
            Feature::Sweep { path_edge, .. } | Feature::Pipe { path_edge, .. } => {
                Self::owner_before(model, path_edge.operation, index).filter(|o| *o != body)
            }
            Feature::Primitive(p) => match &p.shape {
                popkern_history::feature_primitive::PrimitiveShape::SheetFill { edges } => {
                    edges.first().and_then(|e| {
                        Self::owner_before(model, e.operation, index).filter(|o| *o != body)
                    })
                }
                _ => None,
            },
            _ => None,
        }
    }
}
