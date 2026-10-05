//! popkern-project: projects as folders of ordered migration scripts
//! (live-projects-and-history.md).
//!
//! A project is a folder:
//!
//! ```text
//! my-bracket/
//!   popkern.toml            the manifest (name, format, optional base)
//!   migrations/
//!     001_base.rhai         builds the document (model() or assembly())
//!     002_pockets.rhai      let m = previous(); … edits it
//!   imports/vendor.step     an optional base document
//! ```
//!
//! [`ProjectRunner`] builds the document by running every migration on
//! the document the earlier ones left, keeping an in-memory checkpoint
//! after each so an edit reruns only from the first changed migration.
//! The built document carries an applied-migration record
//! ([`ProjectRecord`]) in its `.pkn` envelope, so `popkern build`
//! applies only the pending migrations to a saved document and refuses
//! (or rebuilds) one whose history diverged from the folder. All of
//! this is shared by the CLI, the browser engine and the server: the
//! store is a [`ProjectFs`] — a folder on disk ([`NativeFs`]) or memory
//! ([`MemoryFs`]).

#![forbid(unsafe_code)]

pub mod error;
pub mod files;
pub mod manifest;
pub mod project;
pub mod record;
pub mod runner;
pub mod scaffold;

pub use error::ProjectError;
pub use files::{MemoryFs, NativeFs, ProjectFs, normalize};
pub use manifest::{Gallery, MANIFEST, Manifest, PROJECT_FORMAT};
pub use project::{BaseFile, MIGRATIONS, Migration, Project, normalize_source, sha256_hex};
pub use record::{
    BaseRecord, FeatureRef, MigrationRecord, ProjectRecord, Reconcile, read_record, reconcile,
    to_pkn_bytes,
};
pub use runner::{ProjectBuild, ProjectRunner, VersionReport, VersionStatus};
