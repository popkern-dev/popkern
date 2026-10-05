//! Typed project errors (live-projects-and-history.md §8).

use popkern_rhai::ScriptRunError;

/// Why a project could not be read, built, or reconciled with a saved
/// document.
#[derive(Debug, Clone, thiserror::Error)]
pub enum ProjectError {
    /// `popkern.toml` is malformed.
    #[error("{path}:{line}: {reason}")]
    ManifestInvalid {
        /// The manifest's project-relative path.
        path: String,
        /// 1-based line (0 when the whole file is at fault).
        line: usize,
        /// What is wrong.
        reason: String,
    },
    /// A file in `migrations/` is not named `<number>_<name>.rhai`.
    #[error(
        "migrations/{file}: a migration is named <number>_<name>.rhai (for example 010_pockets.rhai)"
    )]
    MigrationName {
        /// The offending file name.
        file: String,
    },
    /// Two migrations share a number, so their order is ambiguous.
    #[error("migrations {first} and {second} share the number {number}: renumber one")]
    MigrationOrder {
        /// The first file.
        first: String,
        /// The second file.
        second: String,
        /// The shared number.
        number: u64,
    },
    /// The project has no migrations and no base.
    #[error("the project has no migrations: add migrations/001_<name>.rhai")]
    Empty,
    /// A path resolves outside the project folder.
    #[error("{path} is outside the project folder")]
    PathOutsideProject {
        /// The path as written.
        path: String,
    },
    /// A project file does not exist.
    #[error("{path} does not exist in the project")]
    FileNotFound {
        /// The project-relative path.
        path: String,
    },
    /// Reading or writing a project file failed.
    #[error("{path}: {detail}")]
    Io {
        /// The project-relative path.
        path: String,
        /// The underlying error.
        detail: String,
    },
    /// A migration script failed (its error carries the line).
    #[error("{migration}: {error}")]
    MigrationFailed {
        /// The migration's file name.
        migration: String,
        /// The script error.
        error: Box<ScriptRunError>,
        /// Upstream migrations whose content changed since this one
        /// last built — where a vanished label usually comes from.
        upstream_changed: Vec<String>,
    },
    /// The base document could not be loaded.
    #[error("base {path}: {detail}")]
    BaseInvalid {
        /// The base's project-relative path.
        path: String,
        /// Why it could not be loaded.
        detail: String,
    },
    /// A saved document's applied migration differs from the folder's.
    #[error(
        "the document applied {migration} with content {recorded}, but the folder's {migration} \
         is {current}: the history diverged (rebuild, or restore the migration)"
    )]
    ProvenanceDiverged {
        /// The migration (or `base`).
        migration: String,
        /// The content hash the document recorded.
        recorded: String,
        /// The folder's content hash.
        current: String,
    },
    /// A saved document applied a migration the folder no longer has.
    #[error(
        "the document applied {name}, which the folder no longer has: its history cannot be \
         reproduced (rebuild from the folder with --rebuild)"
    )]
    MissingMigration {
        /// The recorded migration name.
        name: String,
    },
    /// The document could not be serialized or read.
    #[error("document: {detail}")]
    Document {
        /// Why.
        detail: String,
    },
}

impl ProjectError {
    /// The class string hosts show (`ManifestInvalid`, `MigrationFailed`, …).
    pub fn kind(&self) -> &'static str {
        match self {
            ProjectError::ManifestInvalid { .. } => "ManifestInvalid",
            ProjectError::MigrationName { .. } => "MigrationName",
            ProjectError::MigrationOrder { .. } => "MigrationOrder",
            ProjectError::Empty => "EmptyProject",
            ProjectError::PathOutsideProject { .. } => "PathOutsideProject",
            ProjectError::FileNotFound { .. } => "FileNotFound",
            ProjectError::Io { .. } => "IoFailed",
            ProjectError::MigrationFailed { .. } => "MigrationFailed",
            ProjectError::BaseInvalid { .. } => "BaseInvalid",
            ProjectError::ProvenanceDiverged { .. } => "ProvenanceDiverged",
            ProjectError::MissingMigration { .. } => "MissingMigration",
            ProjectError::Document { .. } => "DocumentError",
        }
    }
}
