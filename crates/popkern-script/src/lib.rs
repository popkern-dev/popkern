//! `popkern-script` — the scripting surface recorder.
//!
//! Code-first authoring that records into the `popkern-history`
//! feature tree (design D1: feature-tree-as-truth). Every call on a
//! [`Model`] or [`Body`] resolves its selections against the live
//! body, constructs one `Feature`, appends it, replays it, and either
//! serves the new body or removes the feature again and raises a
//! structured [`ScriptError`] (transactional per call).
//!
//! The normative contract is `docs/architecture/scripting-surface.md`;
//! this crate is its Rust skin and the substrate every other skin
//! (Rhai, Python, MCP) wraps. Selections are **data**
//! ([`EdgeSelector`], [`FaceSelector`]) with a canonical spelling, so
//! a document can always be shown back as code.

#![forbid(unsafe_code)]

pub mod bodies;
pub mod body;
pub mod body_recipes;
pub mod body_sweep;
pub mod configs;
pub mod entity_spelling;
pub mod error;
pub mod error_payload;
pub mod feature_edit;
pub mod migration;
pub mod model;
pub mod params;
pub mod plane;
pub mod profile;
mod record;
pub mod scalar;
pub mod selection;
pub mod selection_parse;

pub use body::{Body, End, RevolveEnd};
pub use body_recipes::{AxisSpec, ChamferSpec, RadiusLaw, ShellOptions};
pub use body_sweep::{HelixParams, LoftOptions};
pub use error::{ScriptError, SelectionKind};
pub use error_payload::surface_parameter;
pub use feature_edit::Edit;
pub use migration::ParamInfo;
pub use model::{FeatureHandle, Model};
pub use params::ConfigInfo;
pub use plane::Plane;
pub use popkern_history::{
    BlendDeclineKind, BlendDeclined, BlendFailureReason, BlendFailureRecord, FeatureFailure,
    RemediationRecord,
};
pub use profile::ProfileBuilder;
pub use scalar::Scalar;
pub use selection::{
    EdgePredicate, EdgeSelector, EntityRef, FacePredicate, FaceSelector, SelectionSnapshot,
};
pub use selection_parse::{SelectionParseError, parse_edges, parse_faces};

/// Everything a script needs in scope.
pub mod prelude {
    pub use crate::body::{Body, End, RevolveEnd};
    pub use crate::body_recipes::{AxisSpec, ChamferSpec, RadiusLaw, ShellOptions};
    pub use crate::body_sweep::{HelixParams, LoftOptions};
    pub use crate::error::ScriptError;
    pub use crate::feature_edit::Edit;
    pub use crate::model::{FeatureHandle, Model};
    pub use crate::plane::Plane;
    pub use crate::scalar::Scalar;
    pub use crate::selection::{EdgeSelector, FaceSelector, edges, faces};
    pub use popkern_history::SnapshotOrigin;
    pub use popkern_history::feature::{
        AxisDesc, BasePlane, BoolOpKind, JoinStyleKind, OperationMode, ShellDirectionKind,
        TrimPolicyKind, Vec3Desc,
    };
    pub use popkern_history::feature_blend::{
        BlendProfileKind, CapKindDesc, ChamferKindDesc, ContinuityKind, FilletRadiusLaw,
        InterferencePolicyKind,
    };
    pub use popkern_history::unit::Unit;
    pub use popkern_part::{Dir, SurfaceKind};
}
