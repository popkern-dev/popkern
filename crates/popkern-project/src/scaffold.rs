//! New projects (`popkern new`): the manifest, a first migration that
//! builds a document, and a `.gitignore` for the reserved build folder.

use crate::error::ProjectError;
use crate::files::ProjectFs;
use crate::manifest::{MANIFEST, Manifest};
use crate::project::MIGRATIONS;

/// The first migration of a new project.
pub const FIRST_MIGRATION: &str = "001_base.rhai";

/// The source of a new project's first migration.
pub fn first_migration_source(name: &str) -> String {
    format!(
        "// {name}: the first migration builds the document. Later\n\
         // migrations (002_<name>.rhai, …) start with `let m = previous();`\n\
         // and edit it.\n\
         let m = model();\n\
         m.profile(xy()).rect(0.08, 0.05).sketch(\"plate\").extrude(blind(0.01));\n\
         m\n"
    )
}

/// The source of a later migration's template.
pub fn next_migration_source() -> String {
    "// Edits the document the earlier migrations built.\nlet m = previous();\n\nm\n".to_string()
}

/// Write a new project called `name` into `fs`.
///
/// # Errors
///
/// [`ProjectError::Io`] when the store already holds a manifest (a
/// project is never overwritten) or a write fails.
pub fn create(fs: &dyn ProjectFs, name: &str) -> Result<(), ProjectError> {
    if fs.read(MANIFEST).is_ok() {
        return Err(ProjectError::Io {
            path: MANIFEST.into(),
            detail: "a project already exists here".into(),
        });
    }
    fs.write(MANIFEST, Manifest::new(name).to_text().as_bytes())?;
    fs.write(
        &format!("{MIGRATIONS}/{FIRST_MIGRATION}"),
        first_migration_source(name).as_bytes(),
    )?;
    fs.write(
        ".gitignore",
        b"# Reserved for popkern build output and caches.\n.popkern/\n",
    )?;
    Ok(())
}
