//! Projects in a session (live-projects-and-history.md §5.4): the page
//! hands over the project's files, the shared runner builds them with
//! checkpoints kept across calls (an edit to the last migration reruns
//! only it), and the envelope carries one row per migration beside the
//! usual model description.

use std::rc::Rc;

use serde::Serialize;
use serde_json::{Value, json};

use popkern_project::{MemoryFs, Project, ProjectBuild, ProjectError, ProjectFs, ProjectRunner};
use popkern_rhai::Outcome;
use popkern_rhai::outcome::from_run_error;

use crate::envelope::Envelope;
use crate::session::Session;

/// The runner and in-memory store a session keeps between project
/// builds.
pub(crate) struct ProjectState {
    fs: Rc<MemoryFs>,
    runner: ProjectRunner,
}

impl ProjectState {
    /// Replace the runner's script budgets.
    pub(crate) fn set_options(&mut self, options: popkern_rhai::HostOptions) {
        self.runner.set_options(options);
    }
}

/// One migration's row as the page shows it.
#[derive(Debug, Clone, Serialize)]
pub struct VersionRow {
    /// The migration's name (file stem).
    pub name: String,
    /// `built`, `cached`, `failed` or `stale`.
    pub status: &'static str,
    /// Wall time when it ran.
    pub elapsed_ms: Option<f64>,
    /// Labels of the features it recorded.
    pub features: Vec<String>,
    /// Labels of the earlier features it changed.
    pub edited: Vec<String>,
    /// What it printed.
    pub printed: Vec<String>,
}

/// The project half of an envelope.
#[derive(Debug, Clone, Serialize)]
pub struct ProjectView {
    /// The project's name.
    pub name: String,
    /// One row per migration, in order.
    pub versions: Vec<VersionRow>,
    /// Why the build stopped (`kind`, `message`; for a failed
    /// migration also `migration`, `line`, `column`, `error`,
    /// `upstream_changed`).
    pub error: Option<Value>,
}

/// A project error as the envelope's `project.error`.
pub fn project_error_value(error: &ProjectError) -> Value {
    let mut out = json!({ "kind": error.kind(), "message": error.to_string() });
    if let ProjectError::MigrationFailed {
        migration,
        error,
        upstream_changed,
    } = error
    {
        out["migration"] = json!(migration);
        out["line"] = json!(error.line);
        out["column"] = json!(error.column);
        out["error"] = (*error.payload).clone();
        out["upstream_changed"] = json!(upstream_changed);
    }
    out
}

fn view(name: &str, build: &ProjectBuild) -> ProjectView {
    ProjectView {
        name: name.to_string(),
        versions: build
            .versions
            .iter()
            .map(|v| VersionRow {
                name: v.name.clone(),
                status: v.status.as_str(),
                elapsed_ms: v.elapsed.map(|d| d.as_secs_f64() * 1000.0),
                features: v.features.iter().map(|f| f.label.clone()).collect(),
                edited: v.edited.iter().map(|f| f.label.clone()).collect(),
                printed: v.printed.clone(),
            })
            .collect(),
        error: build.error.as_ref().map(project_error_value),
    }
}

/// The session outcome a project error maps to.
fn outcome_of(error: &ProjectError) -> Outcome {
    match error {
        ProjectError::MigrationFailed { error, .. } => from_run_error(error),
        other => Outcome::GateFailure {
            step: "project".into(),
            detail: other.to_string(),
        },
    }
}

/// Decode the page's project transfer: `{"files": {path: text},
/// "base64": {path: data}}` (binary files such as a STEP base travel
/// as base64).
pub fn decode_project_files(json: &str) -> Result<Vec<(String, Vec<u8>)>, String> {
    let value: Value = serde_json::from_str(json).map_err(|e| format!("project files: {e}"))?;
    let mut files = Vec::new();
    if let Some(map) = value.get("files").and_then(Value::as_object) {
        for (path, text) in map {
            let text = text
                .as_str()
                .ok_or_else(|| format!("project file {path}: text expected"))?;
            files.push((path.clone(), text.as_bytes().to_vec()));
        }
    }
    if let Some(map) = value.get("base64").and_then(Value::as_object) {
        for (path, data) in map {
            let data = data
                .as_str()
                .ok_or_else(|| format!("project file {path}: base64 text expected"))?;
            files.push((
                path.clone(),
                base64_decode(data).map_err(|e| format!("{path}: {e}"))?,
            ));
        }
    }
    Ok(files)
}

/// Standard base64 (RFC 4648, `=` padding optional, whitespace ignored).
pub fn base64_decode(text: &str) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(text.len() * 3 / 4);
    let mut buffer = 0u32;
    let mut bits = 0;
    for c in text.bytes() {
        let v = match c {
            b'A'..=b'Z' => c - b'A',
            b'a'..=b'z' => c - b'a' + 26,
            b'0'..=b'9' => c - b'0' + 52,
            b'+' => 62,
            b'/' => 63,
            b'=' | b' ' | b'\n' | b'\r' | b'\t' => continue,
            _ => return Err(format!("invalid base64 character '{}'", c as char)),
        };
        buffer = (buffer << 6) | u32::from(v);
        bits += 6;
        if bits >= 8 {
            bits -= 8;
            out.push((buffer >> bits) as u8);
            buffer &= (1 << bits) - 1;
        }
    }
    Ok(out)
}

impl Session {
    /// Build the project made of `files` (project-relative path →
    /// bytes). Checkpoints persist across calls, so only migrations
    /// from the first changed one rerun. On success the session adopts
    /// the built document; when a migration fails it adopts the last
    /// good document (earlier migrations' result) and the envelope's
    /// outcome classifies the failure.
    pub fn run_project(&mut self, files: Vec<(String, Vec<u8>)>) -> Envelope {
        let fs = match MemoryFs::new(files) {
            Ok(fs) => fs,
            Err(error) => return self.project_failure(&error, None),
        };
        let state = self.project.get_or_insert_with(|| {
            let fs = Rc::new(MemoryFs::default());
            let store: Rc<dyn ProjectFs> = fs.clone();
            ProjectState {
                runner: ProjectRunner::new(store, self.options),
                fs,
            }
        });
        state.fs.replace_with(&fs);
        let project = match Project::load(state.fs.as_ref()) {
            Ok(project) => project,
            Err(error) => return self.project_failure(&error, None),
        };
        let build = match state.runner.build(&project) {
            Ok(build) => build,
            Err(error) => return self.project_failure(&error, Some(&project.manifest.name)),
        };
        let view = view(&project.manifest.name, &build);
        let printed: Vec<String> = build
            .versions
            .iter()
            .flat_map(|v| v.printed.iter().cloned())
            .collect();
        let elapsed: Option<f64> = build
            .versions
            .iter()
            .filter_map(|v| v.elapsed)
            .map(|d| d.as_secs_f64() * 1000.0)
            .reduce(|a, b| a + b);
        let outcome = build.error.as_ref().map_or(Outcome::Served, outcome_of);
        self.migrations = build
            .record
            .migrations
            .iter()
            .flat_map(|m| {
                m.features
                    .iter()
                    .map(move |f| ((f.part.clone(), f.id), m.name.clone()))
            })
            .collect();
        self.assembly_migrations = build
            .versions
            .iter()
            .filter_map(|v| Some((v.name.clone(), v.assembly_steps?)))
            .collect();
        if build.document.is_some() {
            self.document = build.document;
            self.display.clear();
        }
        let mut envelope = self.describe(outcome, printed, elapsed);
        envelope.project = Some(view);
        envelope
    }

    fn project_failure(&mut self, error: &ProjectError, name: Option<&str>) -> Envelope {
        let mut envelope = self.bare(outcome_of(error));
        envelope.project = Some(ProjectView {
            name: name.unwrap_or_default().to_string(),
            versions: Vec::new(),
            error: Some(project_error_value(error)),
        });
        envelope
    }
}
