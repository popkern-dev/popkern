//! Project hosting (live-projects-and-history.md §5.2): the
//! `previous()` verb a migration script starts from, the migration
//! contract (a later migration ends with the very document `previous()`
//! returned), and sandboxed file verbs — `open`, `save`, `to_step`,
//! `assert_step_roundtrip` — that resolve inside the project through a
//! [`ScriptFiles`] store instead of the host's filesystem.

use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::time::Duration;

use popkern_assembly::Document;
use popkern_history::deferred_certification::ContinuityCertification;
use popkern_script::{Model, ProjectReason, ScriptError};
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::assembly_api::{IntoAssemblyResult, ScriptAssembly};
use crate::error::IntoScriptResult;
use crate::host::{ScriptHost, ScriptRunError, authored};
use crate::model_api::ScriptModel;

/// The files a project script may read and write, by project-relative
/// path. Implementations sandbox: a path outside the project fails
/// with [`ProjectReason::PathOutsideProject`].
pub trait ScriptFiles {
    /// The bytes of `path`.
    fn read(&self, path: &str) -> Result<Vec<u8>, ScriptError>;
    /// Replace `path` with `bytes`.
    fn write(&self, path: &str, bytes: &[u8]) -> Result<(), ScriptError>;
    /// Where `path` lives on a filesystem, when the store has one (the
    /// STEP round-trip gate re-imports through the filesystem).
    fn native_path(&self, path: &str) -> Result<Option<PathBuf>, ScriptError>;
}

/// What a migration run produced.
pub struct MigrationOutcome {
    /// The document after the migration.
    pub document: Document,
    /// Lines the script printed, in order.
    pub printed: Vec<String>,
    /// Wall time of the run when a clock was available.
    pub elapsed: Option<Duration>,
}

impl std::fmt::Debug for MigrationOutcome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MigrationOutcome")
            .field("printed", &self.printed)
            .field("elapsed", &self.elapsed)
            .finish()
    }
}

fn project_error(reason: ProjectReason, path: Option<&str>, detail: String) -> ScriptError {
    ScriptError::Project {
        reason,
        path: path.map(str::to_string),
        detail,
    }
}

impl ScriptHost {
    /// A host whose file verbs read and write `files` (a project).
    pub fn with_files(options: crate::HostOptions, files: Rc<dyn ScriptFiles>) -> ScriptHost {
        let mut host = ScriptHost::with_options(options);
        register_project_files(&mut host.engine, files, options.continuity_certification);
        host
    }

    /// Run one project migration named `script`. `previous` is the
    /// document the earlier migrations left (`None` for a first
    /// migration with no base): `previous()` returns it, and the run
    /// must end with that same document — the script may return it,
    /// return nothing, or leave it in scope; other models it creates
    /// are tools. A first migration authors its document as any
    /// script does.
    ///
    /// # Errors
    ///
    /// The script's own error, or a `DocumentError` whose reason is
    /// `no_previous_document`, `document_replaced` or `no_document`.
    pub fn run_migration(
        &mut self,
        source: &str,
        script: &str,
        previous: Option<Document>,
    ) -> Result<MigrationOutcome, ScriptRunError> {
        let certification = self.options.continuity_certification;
        let handle = previous.map(|document| wrap(document, certification));
        let given = handle.clone();
        let name = script.to_string();
        self.engine.register_fn(
            "previous",
            move || -> Result<Dynamic, Box<EvalAltResult>> {
                match &given {
                    Some(h) => Ok(h.clone()),
                    None => Err::<Dynamic, _>(project_error(
                        ProjectReason::NoPreviousDocument,
                        None,
                        format!(
                            "{name} has no previous document: it is the project's first \
                             migration (and the project has no base), so it starts the \
                             document with model() or assembly()"
                        ),
                    ))
                    .script(),
                }
            },
        );
        let ast = self.compile(source, script)?;
        let run = self.eval(&ast, script)?;
        let document = match handle {
            Some(previous) => {
                if is_document(&run.returned) && !same_document(&run.returned, &previous) {
                    return Err(run_error(
                        script,
                        project_error(
                            ProjectReason::DocumentReplaced,
                            None,
                            format!(
                                "{script} returned a different document than previous(): a \
                                 later migration edits the project's document — start from \
                                 previous() and return it (or nothing)"
                            ),
                        ),
                    ));
                }
                unwrap(previous)
            }
            None => match authored(run.returned, run.in_scope) {
                (Some(model), _) => Some(Document::Part(Box::new(model))),
                (None, Some(assembly)) => Some(Document::Assembly(Box::new(assembly))),
                (None, None) => None,
            },
        };
        let document = document.ok_or_else(|| {
            run_error(
                script,
                project_error(
                    ProjectReason::NoDocument,
                    None,
                    format!(
                        "{script} authored no document: the first migration creates one with \
                         model() or assembly() and keeps it in a variable (or returns it)"
                    ),
                ),
            )
        })?;
        Ok(MigrationOutcome {
            document,
            printed: self.printed.borrow().clone(),
            elapsed: run.elapsed,
        })
    }
}

/// A script handle on `document`.
fn wrap(document: Document, certification: ContinuityCertification) -> Dynamic {
    match document {
        Document::Part(mut model) => {
            model.set_continuity_certification(certification);
            Dynamic::from(ScriptModel::new(*model))
        }
        Document::Assembly(assembly) => {
            Dynamic::from(ScriptAssembly::from_assembly(*assembly, certification))
        }
    }
}

/// The document a handle holds.
fn unwrap(handle: Dynamic) -> Option<Document> {
    if let Some(model) = handle.clone().try_cast::<ScriptModel>() {
        return Some(Document::Part(Box::new(model.into_inner())));
    }
    handle
        .try_cast::<ScriptAssembly>()
        .map(|a| Document::Assembly(Box::new(a.into_inner())))
}

fn is_document(value: &Dynamic) -> bool {
    value.is::<ScriptModel>() || value.is::<ScriptAssembly>()
}

fn same_document(a: &Dynamic, b: &Dynamic) -> bool {
    match (
        a.clone().try_cast::<ScriptModel>(),
        b.clone().try_cast::<ScriptModel>(),
    ) {
        (Some(x), Some(y)) => return x.same_document(&y),
        (Some(_), None) | (None, Some(_)) => return false,
        (None, None) => {}
    }
    match (
        a.clone().try_cast::<ScriptAssembly>(),
        b.clone().try_cast::<ScriptAssembly>(),
    ) {
        (Some(x), Some(y)) => Rc::ptr_eq(&x.state, &y.state),
        _ => false,
    }
}

/// A host-level failure as a run error (no script position).
fn run_error(script: &str, error: ScriptError) -> ScriptRunError {
    let payload = error.payload();
    ScriptRunError {
        kind: payload
            .get("kind")
            .and_then(serde_json::Value::as_str)
            .unwrap_or("DocumentError")
            .to_string(),
        message: error.to_string(),
        script: script.to_string(),
        line: None,
        column: None,
        payload: Box::new(payload),
    }
}

/// Register `previous()` for a plain (non-project) run: there is no
/// previous document outside a project migration.
pub(crate) fn register_outside_project(engine: &mut Engine) {
    engine.register_fn("previous", || -> Result<Dynamic, Box<EvalAltResult>> {
        Err::<Dynamic, _>(project_error(
            ProjectReason::NoPreviousDocument,
            None,
            "previous() returns the document earlier migrations built; it is available \
             only in a project migration (popkern build, or a project in popkern live)"
                .into(),
        ))
        .script()
    });
}

/// The file stem of `path` (the STEP product name).
fn stem(path: &str) -> String {
    Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("model")
        .to_string()
}

/// Re-register the file verbs to go through `files`.
fn register_project_files(
    engine: &mut Engine,
    files: Rc<dyn ScriptFiles>,
    certification: ContinuityCertification,
) {
    let f = Rc::clone(&files);
    engine.register_fn("open", move |path: &str| {
        let bytes = f.read(path).script()?;
        let mut model = Model::from_pkn_bytes(&bytes).script()?;
        model.set_continuity_certification(certification);
        Ok::<_, Box<EvalAltResult>>(ScriptModel::new(model))
    });
    let f = Rc::clone(&files);
    engine.register_fn("save", move |m: &mut ScriptModel, path: &str| {
        let bytes = m.borrow_mut().to_pkn_bytes().script()?;
        f.write(path, &bytes).script()
    });
    let f = Rc::clone(&files);
    engine.register_fn("to_step", move |m: &mut ScriptModel, path: &str| {
        let bytes = m.borrow_mut().to_step_bytes(&stem(path)).script()?;
        f.write(path, &bytes).script()
    });
    let f = Rc::clone(&files);
    engine.register_fn(
        "assert_step_roundtrip",
        move |m: &mut ScriptModel, path: &str| match f.native_path(path).script()? {
            Some(native) => m.borrow_mut().assert_step_roundtrip(native).script(),
            None => Err::<(), _>(ScriptError::Io {
                path: path.to_string(),
                detail: "the STEP round-trip gate re-imports through a filesystem; this \
                         project's files are held in memory"
                    .into(),
            })
            .script(),
        },
    );
    let f = Rc::clone(&files);
    engine.register_fn("save", move |a: &mut ScriptAssembly, path: &str| {
        let bytes = a.borrow_mut().assembly.to_pkn_bytes().assembly()?;
        f.write(path, &bytes).script()
    });
    let f = files;
    engine.register_fn("to_step", move |a: &mut ScriptAssembly, path: &str| {
        let bytes = a
            .borrow_mut()
            .assembly
            .to_step_bytes(&stem(path))
            .assembly()?;
        f.write(path, &bytes).script()
    });
}

#[cfg(test)]
#[path = "project_host_tests.rs"]
mod tests;
