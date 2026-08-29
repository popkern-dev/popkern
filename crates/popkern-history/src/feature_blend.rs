//! Serialisable blend recipes for [`Feature::Fillet`](crate::Feature::Fillet)
//! and [`Feature::Chamfer`](crate::Feature::Chamfer).
//!
//! These descriptors mirror the `popkern-blend` option types with
//! serde (the same pattern as [`JoinStyleKind`](crate::JoinStyleKind))
//! and are converted to kernel types at the replay boundary. Every
//! field is `#[serde(default)]`-compatible so trees saved before the
//! recipe existed load unchanged: the defaults reproduce exactly what
//! replay used to hardcode (`FilletOptions::default()`, an equal
//! chamfer).
//!
//! Radii and distances are [`ParameterId`]s so a variable-radius law
//! or a two-distance chamfer stays parametric (driven by the parameter
//! table, expressions, configurations) exactly like the constant case.

use crate::parameter::ParameterId;
use popkern_core::naming::EntityName;
use serde::{Deserialize, Serialize};

/// Everything beyond "constant radius" that a fillet feature carries:
/// the radius law, the requested continuity, the end-cap style, the
/// cross-section profile, and the strip-interference policy.
///
/// The default recipe is a constant-radius, `G1`, open-ended,
/// circular, transactionally-declining fillet — the historical replay
/// behaviour.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct FilletRecipeDesc {
    /// How the radius varies along each edge. The feature's `radius`
    /// parameter is the constant value, or the start value of a law.
    #[serde(default)]
    pub radius_law: FilletRadiusLaw,
    /// Requested blend-boundary continuity.
    #[serde(default)]
    pub continuity: ContinuityKind,
    /// End-cap style for open fillet ends.
    #[serde(default)]
    pub cap: CapKindDesc,
    /// Cross-section profile.
    #[serde(default)]
    pub profile: BlendProfileKind,
    /// Resolution policy for non-completable strip interference.
    #[serde(default)]
    pub interference: InterferencePolicyKind,
}

/// Radius law along a filleted edge (mirrors
/// [`popkern_blend::FilletRadius`]).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub enum FilletRadiusLaw {
    /// Constant radius (the feature's `radius` parameter).
    #[default]
    Constant,
    /// Linear from the feature's `radius` at the edge start to `end`
    /// at the edge end.
    Linear {
        /// Radius at the end vertex (resolved via parameter table).
        end: ParameterId,
    },
    /// Radius defined at explicit stations along the edge.
    ByStations(Vec<RadiusStation>),
}

/// One station of a [`FilletRadiusLaw::ByStations`] law.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RadiusStation {
    /// Normalised edge parameter in `[0, 1]`.
    pub t: f64,
    /// Radius at this station (resolved via parameter table).
    pub radius: ParameterId,
}

/// Requested blend-boundary continuity (mirrors
/// [`popkern_blend::ContinuityLevel`]).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum ContinuityKind {
    /// Positional only.
    G0,
    /// Tangent continuity — the rolling-ball arc sweep (default).
    #[default]
    G1,
    /// Curvature continuity — the degree-5 sweep.
    G2,
}

/// End-cap style for open fillet ends (mirrors
/// [`popkern_blend::CapKind`]; the kernel's `UserDefined` placeholder
/// is not representable in a document).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum CapKindDesc {
    /// No cap surface (default).
    #[default]
    Open,
    /// Planar cap perpendicular to the spine.
    Flat,
    /// Hemispherical cap.
    Rounded,
    /// Tangent-matched NURBS cap.
    TangentMatch,
}

/// Cross-section profile (mirrors [`popkern_blend::BlendProfile`]).
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub enum BlendProfileKind {
    /// Circular rolling-ball section (default).
    #[default]
    Circular,
    /// Conic section with shoulder fullness `rho` in `(0, 1)`.
    Conic {
        /// Shoulder fullness; `0.5` is parabolic.
        rho: f64,
    },
    /// Curvature-continuous degree-5 section (equivalent to
    /// requesting `G2`).
    CurvatureContinuous,
}

/// Strip-interference policy (mirrors
/// [`popkern_blend::InterferencePolicy`]).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum InterferencePolicyKind {
    /// Decline the whole operation transactionally (default).
    #[default]
    Decline,
    /// Drop interfering pairs and serve the rest, with warnings.
    ServeNonInterfering,
}

/// Chamfer geometry beyond the equal-distance bevel (mirrors
/// [`popkern_blend::ChamferKind`]). The feature's `distance` parameter
/// is the equal distance, or the first distance of the other kinds.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub enum ChamferKindDesc {
    /// Equal-distance 45° bevel (default).
    #[default]
    Equal,
    /// Two distances: the feature's `distance` measured on `face`,
    /// `distance_2` on the other support.
    TwoDistance {
        /// Reference face for the first distance, by stable name.
        face: EntityName,
        /// Distance on the second face (resolved via parameter table).
        distance_2: ParameterId,
    },
    /// The feature's `distance` measured on `face`, at `angle`
    /// (radians) from it.
    DistanceAngle {
        /// Reference face for the distance measurement, by stable name.
        face: EntityName,
        /// Angle from the reference face (resolved via parameter table).
        angle: ParameterId,
    },
}

impl ContinuityKind {
    /// The kernel continuity level this descriptor requests.
    pub fn to_kernel(self) -> popkern_blend::ContinuityLevel {
        match self {
            ContinuityKind::G0 => popkern_blend::ContinuityLevel::G0,
            ContinuityKind::G1 => popkern_blend::ContinuityLevel::G1,
            ContinuityKind::G2 => popkern_blend::ContinuityLevel::G2,
        }
    }
}

impl CapKindDesc {
    /// The kernel cap kind this descriptor requests.
    pub fn to_kernel(self) -> popkern_blend::CapKind {
        match self {
            CapKindDesc::Open => popkern_blend::CapKind::Open,
            CapKindDesc::Flat => popkern_blend::CapKind::Flat,
            CapKindDesc::Rounded => popkern_blend::CapKind::Rounded,
            CapKindDesc::TangentMatch => popkern_blend::CapKind::TangentMatch,
        }
    }
}

impl BlendProfileKind {
    /// The kernel profile this descriptor requests.
    pub fn to_kernel(self) -> popkern_blend::BlendProfile {
        match self {
            BlendProfileKind::Circular => popkern_blend::BlendProfile::Circular,
            BlendProfileKind::Conic { rho } => popkern_blend::BlendProfile::Conic { rho },
            BlendProfileKind::CurvatureContinuous => {
                popkern_blend::BlendProfile::CurvatureContinuous
            }
        }
    }
}

impl InterferencePolicyKind {
    /// The kernel policy this descriptor requests.
    pub fn to_kernel(self) -> popkern_blend::InterferencePolicy {
        match self {
            InterferencePolicyKind::Decline => popkern_blend::InterferencePolicy::Decline,
            InterferencePolicyKind::ServeNonInterfering => {
                popkern_blend::InterferencePolicy::ServeNonInterfering
            }
        }
    }
}

impl FilletRecipeDesc {
    /// Kernel fillet options for this recipe: the kernel defaults with
    /// the recipe's continuity, cap, profile, and interference policy
    /// applied. The engineering toggles that `FilletOptions` also
    /// carries (Gregory variants, soft supports, apex construction)
    /// are not document state and keep their defaults.
    pub fn to_kernel_options(&self) -> popkern_blend::FilletOptions {
        popkern_blend::FilletOptions {
            continuity: self.continuity.to_kernel(),
            cap_kind: self.cap.to_kernel(),
            profile: self.profile.to_kernel(),
            interference_policy: self.interference.to_kernel(),
            ..popkern_blend::FilletOptions::default()
        }
    }

    /// Every parameter the recipe references, beyond the feature's own
    /// `radius`, with a field label for dependency reporting.
    pub(crate) fn param_ids(&self) -> Vec<(&'static str, ParameterId)> {
        match &self.radius_law {
            FilletRadiusLaw::Constant => vec![],
            FilletRadiusLaw::Linear { end } => vec![("end_radius", *end)],
            FilletRadiusLaw::ByStations(stations) => stations
                .iter()
                .map(|s| ("station_radius", s.radius))
                .collect(),
        }
    }
}

impl ChamferKindDesc {
    /// Every parameter the kind references, beyond the feature's own
    /// `distance`, with a field label for dependency reporting.
    pub(crate) fn param_ids(&self) -> Vec<(&'static str, ParameterId)> {
        match self {
            ChamferKindDesc::Equal => vec![],
            ChamferKindDesc::TwoDistance { distance_2, .. } => vec![("distance_2", *distance_2)],
            ChamferKindDesc::DistanceAngle { angle, .. } => vec![("angle", *angle)],
        }
    }

    /// The reference face this kind names, if any (for entity
    /// reference reporting and health checks).
    pub fn reference_face(&self) -> Option<&EntityName> {
        match self {
            ChamferKindDesc::Equal => None,
            ChamferKindDesc::TwoDistance { face, .. }
            | ChamferKindDesc::DistanceAngle { face, .. } => Some(face),
        }
    }
}

#[cfg(test)]
#[path = "feature_blend_tests.rs"]
mod tests;
