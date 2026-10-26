//! Generative design as a part feature (PK-956,
//! docs/architecture/generative-features-plan.md): a part names its
//! design space, the faces to hold, load and keep, and a material; the
//! analysis solver's density optimisation finds the stiffest design at
//! a volume fraction, and its surface becomes a planar-faced body later
//! features build on.

pub mod cache;
pub mod components;
pub mod facet;
pub mod generate;
pub mod history;
mod history_store;
pub mod naming;
pub mod request;
pub mod surface;

pub use cache::{
    Update, defer_with, design_history, export_cache, export_histories, export_touched,
    generate_cached, import_cache, import_histories, observe_with, touched_designs,
};
pub use generate::{
    GenerativeError, GenerativeReport, GenerativeRequest, GenerativeResult, LoadGroup,
    SupportGroup, generate,
};
pub use history::{DesignHistory, DisplayMesh, HistoryFrame};
pub use naming::name_design;

/// The support kinds a request uses (the analysis setup's).
pub use popkern_analysis::setup::SupportKind;
/// The build direction a request prints along (the solver's).
pub use popkern_analysis::structural::generative::print::BuildDirection;
