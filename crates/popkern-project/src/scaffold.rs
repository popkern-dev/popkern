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

/// A Rhai string literal for `text`.
fn rhai_string(text: &str) -> String {
    let mut out = String::from("\"");
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// The anchor feature: by label, or by 1-based timeline position when
/// the label is not unique.
fn anchor(label: &str, position: usize, unique: bool) -> String {
    if unique {
        format!("m.feature({})", rhai_string(label))
    } else {
        format!("m.features()[{}]", position.saturating_sub(1))
    }
}

/// The source of a migration that inserts one feature right after the
/// feature labelled `label` — "insert a feature here" at the history
/// playhead (live-projects-and-history.md §12, PK-666). The anchor is
/// found by label; when the label is not `unique` in the tree, by its
/// 1-based timeline `position` instead (exact for the tree it was made
/// from). The closure must record exactly one feature on `body`;
/// features after the anchor replay after it.
pub fn insert_after_source(label: &str, position: usize, unique: bool) -> String {
    let anchor = anchor(label, position, unique);
    format!(
        "// Inserts one feature after {} (step {position} of the history).\n\
         // Record exactly one feature on `body`, for example\n\
         //   body.chamfer(edges().at_height(1.0), 0.05)\n\
         // The features after it replay after the new one.\n\
         let m = previous();\n\
         m.insert_after({anchor}, |m, body| {{\n    body\n}});\n\
         m\n",
        rhai_string(label)
    )
}

/// As [`insert_after_source`] for part `part` of an assembly: the
/// insert runs inside `a.edit_part(part, |m| …)`, so every instance of
/// the part follows and the assembly's mates re-solve.
pub fn insert_after_part_source(part: &str, label: &str, position: usize, unique: bool) -> String {
    let anchor = anchor(label, position, unique);
    format!(
        "// Inserts one feature into part {} after {} (step {position} of its history).\n\
         // Record exactly one feature on `body`, for example\n\
         //   body.chamfer(edges().at_height(1.0), 0.05)\n\
         // Every instance of the part follows; its later features replay after it.\n\
         let a = previous();\n\
         a.edit_part({}, |m| {{\n    m.insert_after({anchor}, |m, body| {{\n        body\n    }});\n}});\n\
         a\n",
        rhai_string(part),
        rhai_string(label),
        rhai_string(part)
    )
}

/// The file stem for an insert-after migration (`insert_after_<label>`,
/// letters, digits, `_` and `-` only).
pub fn insert_after_name(label: &str) -> String {
    let cleaned: String = label
        .to_lowercase()
        .chars()
        .map(|c| if c.is_whitespace() { '_' } else { c })
        .filter(|c| c.is_ascii_alphanumeric() || *c == '_' || *c == '-')
        .collect();
    if cleaned.is_empty() {
        "insert_feature".to_string()
    } else {
        format!("insert_after_{cleaned}")
    }
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
