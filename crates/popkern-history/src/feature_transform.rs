//! Serialisable rigid motion for [`Feature::Transform`](crate::Feature::Transform)
//! (PK-362, scripted-modeling phase 1b).
//!
//! A transform feature moves the whole current body by a similarity
//! motion — the semantics of `popkern-part`'s `translate` / `rotate_*`
//! / `transformed`. Translation offsets and rotation angles are
//! [`ParameterId`]s so the motion stays parametric; a rotation axis is
//! an [`AxisDesc`] so it may follow a model edge; the general arm is a
//! literal 4×4 matrix and inherits `BRep::transform`'s similarity gate
//! (non-uniform scale and shear are declined transactionally because
//! analytic surfaces do not survive them).

use crate::error::ReplayError;
use crate::feature::AxisDesc;
use crate::parameter::{ParameterId, ParameterTable};
use crate::resolve::resolve_axis;
use popkern_core::geometry::transform::Transform3;
use popkern_core::naming::EntityName;
use popkern_core::topology::BRep;
use serde::{Deserialize, Serialize};

/// The motion a [`Feature::Transform`](crate::Feature::Transform) applies.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum RigidMotionDesc {
    /// Translate by `(dx, dy, dz)` (each resolved via parameter table).
    Translate {
        /// Offset along world X.
        dx: ParameterId,
        /// Offset along world Y.
        dy: ParameterId,
        /// Offset along world Z.
        dz: ParameterId,
    },
    /// Rotate by `angle` (radians, resolved via parameter table) about
    /// `axis`, right-handed about the axis direction.
    Rotate {
        /// The rotation axis (a model edge or a literal line).
        axis: AxisDesc,
        /// Rotation angle in radians.
        angle: ParameterId,
    },
    /// A literal row-major 4×4 homogeneous matrix. Must be a similarity
    /// (rigid motion, optional mirror, uniform scale); replay declines
    /// anything else.
    Matrix([[f64; 4]; 4]),
}

impl RigidMotionDesc {
    /// Every parameter the motion references, with a field label for
    /// dependency reporting.
    pub(crate) fn param_ids(&self) -> Vec<(&'static str, ParameterId)> {
        match self {
            RigidMotionDesc::Translate { dx, dy, dz } => {
                vec![("dx", *dx), ("dy", *dy), ("dz", *dz)]
            }
            RigidMotionDesc::Rotate { angle, .. } => vec![("angle", *angle)],
            RigidMotionDesc::Matrix(_) => vec![],
        }
    }

    /// Every stable entity name the motion references (a rotation
    /// axis edge).
    pub(crate) fn entity_names(&self) -> Vec<&EntityName> {
        match self {
            RigidMotionDesc::Rotate {
                axis: AxisDesc::Edge(name),
                ..
            } => vec![name],
            RigidMotionDesc::Rotate { .. }
            | RigidMotionDesc::Translate { .. }
            | RigidMotionDesc::Matrix(_) => vec![],
        }
    }

    /// Mutable access to the rotation-axis edge name, if any (for
    /// pattern-instance remapping).
    pub(crate) fn entity_name_mut(&mut self) -> Option<&mut EntityName> {
        match self {
            RigidMotionDesc::Rotate {
                axis: AxisDesc::Edge(name),
                ..
            } => Some(name),
            _ => None,
        }
    }

    /// Resolve the motion against the current body and parameter table
    /// into a concrete [`Transform3`].
    ///
    /// # Errors
    ///
    /// Parameter lookup or axis resolution failures, and a
    /// [`RigidMotionDesc::Matrix`] that is not a similarity.
    pub(crate) fn resolve(
        &self,
        brep: &BRep,
        param_table: &ParameterTable,
    ) -> Result<Transform3, ReplayError> {
        let get = |id: ParameterId| -> Result<f64, ReplayError> {
            param_table
                .get(id)
                .map(|p| p.value)
                .ok_or_else(|| ReplayError::FeatureFailed {
                    label: "transform".into(),
                    error: format!("parameter {id:?} not found").into(),
                })
        };
        match self {
            RigidMotionDesc::Translate { dx, dy, dz } => {
                Ok(Transform3::translation(get(*dx)?, get(*dy)?, get(*dz)?))
            }
            RigidMotionDesc::Rotate { axis, angle } => {
                let (origin, direction) = resolve_axis(axis, brep)?;
                Transform3::rotation_around_axis(origin, direction, get(*angle)?).map_err(|e| {
                    ReplayError::FeatureFailed {
                        label: "transform".into(),
                        error: e.to_string().into(),
                    }
                })
            }
            RigidMotionDesc::Matrix(matrix) => {
                let t = Transform3 { matrix: *matrix };
                if t.similarity_scale().is_none() {
                    return Err(ReplayError::FeatureFailed {
                        label: "transform".into(),
                        error: "transform matrix is not a similarity (rigid motion, optional \
                                mirror, uniform scale); non-uniform scale and shear do not \
                                preserve analytic surface forms"
                            .into(),
                    });
                }
                Ok(t)
            }
        }
    }
}

#[cfg(test)]
#[path = "feature_transform_tests.rs"]
mod tests;
