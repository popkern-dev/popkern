//! Neutral-plane descriptor for [`Feature::Draft`](crate::Feature::Draft)
//! (PK-363).
//!
//! A draft's neutral plane — the plane whose intersection with each
//! drafted face stays put — was historically a model face by stable
//! name. `popkern-part`'s [`Part::draft`](popkern_part::Part::draft)
//! instead derives an implicit plane tied to no face: perpendicular to
//! the pull direction through the part's lowest point along it. That
//! recipe had no document spelling, so a Part-style draft could not be
//! recorded in a feature tree. [`NeutralPlaneDesc::Literal`] closes the
//! gap; the pull direction is the plane normal in both arms, exactly
//! as replay always treated the face arm.
//!
//! **Serialisation.** The enum is `#[serde(untagged)]` so that the
//! face arm serialises as the bare [`EntityName`] object version-3
//! documents already carry — every pre-PK-363 `.pkn` loads unchanged
//! with no envelope bump — and a literal plane serialises as an
//! `{origin, normal}` object. The two shapes share no field name, so
//! the untagged dispatch is unambiguous.

use crate::error::ReplayError;
use crate::resolve::resolve_face;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::naming::EntityName;
use popkern_core::topology::BRep;
use serde::{Deserialize, Serialize};

/// The neutral plane of a draft feature.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum NeutralPlaneDesc {
    /// A planar model face, by stable name. The plane's normal is
    /// the face's outward normal (the surface orientation flag is
    /// honoured).
    Face(EntityName),
    /// A literal plane defined by an origin point and a normal. The
    /// normal need not be unit length; it is normalised at replay
    /// and a zero-length normal declines transactionally.
    Literal {
        /// A point on the neutral plane, in world space.
        origin: [f64; 3],
        /// The plane normal in world space (also the pull direction).
        normal: [f64; 3],
    },
}

impl NeutralPlaneDesc {
    /// The model face this descriptor references, if any.
    pub fn entity_name(&self) -> Option<&EntityName> {
        match self {
            NeutralPlaneDesc::Face(name) => Some(name),
            NeutralPlaneDesc::Literal { .. } => None,
        }
    }

    /// Mutable access to the referenced face name, if any (pattern
    /// remapping).
    pub(crate) fn entity_name_mut(&mut self) -> Option<&mut EntityName> {
        match self {
            NeutralPlaneDesc::Face(name) => Some(name),
            NeutralPlaneDesc::Literal { .. } => None,
        }
    }

    /// Resolve to a world-space `(origin, unit normal)` pair against
    /// the body the draft will be applied to.
    ///
    /// # Errors
    ///
    /// The face arm fails when the name does not resolve or the face
    /// is not planar; the literal arm fails on a zero-length normal.
    pub(crate) fn resolve(&self, brep: &BRep) -> Result<(Point3, Vec3), ReplayError> {
        match self {
            NeutralPlaneDesc::Face(name) => {
                let face = resolve_face(name, brep)?;
                crate::replay::replay_fns::extract_plane_from_face(brep, face)
            }
            NeutralPlaneDesc::Literal { origin, normal } => {
                let n = Vec3::new(normal[0], normal[1], normal[2]);
                let len = n.length();
                if !(len.is_finite() && len > 0.0) {
                    return Err(ReplayError::FeatureFailed {
                        label: "draft".into(),
                        error: format!(
                            "literal neutral plane normal {normal:?} is degenerate (length {len})"
                        ),
                    });
                }
                Ok((
                    Point3::new(origin[0], origin[1], origin[2]),
                    n * (1.0 / len),
                ))
            }
        }
    }
}

#[cfg(test)]
#[path = "feature_draft_tests.rs"]
mod tests;
