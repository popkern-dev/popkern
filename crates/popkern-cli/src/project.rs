//! `popkern build <dir>` and `popkern new <dir>` (live-projects-and-history.md
//! §5.4, §6.2): build a project folder's migrations into one document,
//! applying only the pending migrations to a saved `--pkn` document,
//! and scaffold new projects.

use std::path::Path;
use std::rc::Rc;

use popkern_assembly::Document;
use popkern_project::{
    NativeFs, Project, ProjectBuild, ProjectError, ProjectFs, ProjectRunner, Reconcile,
    read_record, reconcile, to_pkn_bytes,
};
use popkern_rhai::HostOptions;
use serde_json::{Value, json};

use crate::args::Cli;

/// Exit code: the project folder itself is malformed (manifest,
/// migration names, paths, a missing applied migration).
pub const EXIT_PROJECT: u8 = 5;
/// Exit code: `--strict` and the saved document's history diverged.
pub const EXIT_DIVERGED: u8 = 6;

/// What a build command did.
struct Summary {
    lines: Vec<String>,
    json: Value,
    code: u8,
}

fn error_json(error: &ProjectError) -> Value {
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

/// The exit code a build error maps to.
fn error_code(error: &ProjectError) -> u8 {
    match error {
        ProjectError::MigrationFailed { error, .. } => {
            popkern_rhai::outcome::from_run_error(error).exit_code()
        }
        ProjectError::ProvenanceDiverged { .. } => EXIT_DIVERGED,
        ProjectError::Document { .. } | ProjectError::Io { .. } => 4,
        _ => EXIT_PROJECT,
    }
}

fn render_error(error: &ProjectError) -> String {
    match error {
        ProjectError::MigrationFailed {
            error,
            upstream_changed,
            ..
        } => {
            let mut text = format!("failed: {}", error.render());
            if !upstream_changed.is_empty() {
                text.push_str(&format!(
                    "\n  changed upstream since it last built: {}",
                    upstream_changed.join(", ")
                ));
            }
            text
        }
        other => format!("failed: {other}"),
    }
}

fn describe(build: &ProjectBuild) -> (Vec<String>, Value) {
    let mut lines = Vec::new();
    let mut versions = Vec::new();
    for v in &build.versions {
        let ms = v.elapsed.map(|d| d.as_secs_f64() * 1000.0);
        let timing = ms.map(|ms| format!("  {:.0} ms", ms)).unwrap_or_default();
        lines.push(format!(
            "  {:<28} {:<7} +{} feature(s), {} edited{timing}",
            v.name,
            v.status.as_str(),
            v.features.len(),
            v.edited.len()
        ));
        for printed in &v.printed {
            lines.push(format!("      {printed}"));
        }
        versions.push(json!({
            "name": v.name, "sha256": v.sha256, "status": v.status.as_str(),
            "elapsed_ms": ms, "printed": v.printed,
            "features": v.features.iter().map(|f| &f.label).collect::<Vec<_>>(),
            "edited": v.edited.iter().map(|f| &f.label).collect::<Vec<_>>(),
        }));
    }
    (lines, Value::Array(versions))
}

/// Load, replay, and return a saved document's bytes as a document.
fn load_saved(bytes: &[u8]) -> Result<Document, ProjectError> {
    let mut document = Document::from_pkn_bytes(bytes).map_err(|e| ProjectError::Document {
        detail: e.to_string(),
    })?;
    for (_, model) in document.models_mut() {
        model.replay().map_err(|e| ProjectError::Document {
            detail: e.to_string(),
        })?;
    }
    Ok(document)
}

fn run_build(cli: &Cli) -> Result<Summary, ProjectError> {
    let fs = Rc::new(NativeFs::new(&cli.script)?);
    let project = Project::load(fs.as_ref())?;
    let store: Rc<dyn ProjectFs> = fs.clone();
    let mut runner = ProjectRunner::new(
        store,
        HostOptions {
            max_operations: cli.max_ops,
            timeout: cli.timeout,
            ..HostOptions::default()
        },
    );
    let mut lines = vec![format!(
        "{} ({} migration(s))",
        project.manifest.name,
        project.migrations.len()
    )];
    let saved = match &cli.pkn {
        Some(path) if path.exists() => Some(std::fs::read(path).map_err(|e| ProjectError::Io {
            path: path.display().to_string(),
            detail: e.to_string(),
        })?),
        _ => None,
    };
    let mut build = match &saved {
        None => runner.build(&project)?,
        Some(bytes) => {
            let record = read_record(bytes)?;
            match reconcile(record.as_ref(), &project) {
                _ if cli.rebuild => {
                    lines.push("  rebuilding from the folder (--rebuild)".into());
                    runner.build(&project)?
                }
                Reconcile::UpToDate | Reconcile::Pending { .. } => {
                    let applied = record.map(|r| r.migrations).unwrap_or_default();
                    let from = applied.len();
                    lines.push(format!(
                        "  the saved document applied {from} migration(s); {} pending",
                        project.migrations.len() - from
                    ));
                    runner.build_from(&project, from, Some((load_saved(bytes)?, applied)))?
                }
                Reconcile::Unrecorded => {
                    lines.push(
                        "  the saved document has no applied-migration record: rebuilding".into(),
                    );
                    runner.build(&project)?
                }
                Reconcile::Diverged(error) if cli.strict => return Err(error),
                Reconcile::Diverged(error) => {
                    lines.push(format!("  {error}\n  rebuilding from the folder"));
                    runner.build(&project)?
                }
                Reconcile::Missing(error) => return Err(error),
            }
        }
    };
    let (version_lines, versions) = describe(&build);
    lines.extend(version_lines);
    let mut outputs = Vec::new();
    let mut code = 0;
    let mut error_value = Value::Null;
    let mut stl = None;
    if let Some(error) = &build.error {
        lines.push(render_error(error));
        error_value = error_json(error);
        code = error_code(error);
    } else if let Some(document) = build.document.as_mut() {
        if let Some(path) = &cli.pkn {
            let bytes = to_pkn_bytes(document, &build.record)?;
            write(path, &bytes)?;
            outputs.push(path.display().to_string());
        }
        if let Some(path) = &cli.step {
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("model")
                .to_string();
            let bytes = document
                .to_step_bytes(&stem)
                .map_err(|e| ProjectError::Document {
                    detail: e.to_string(),
                })?;
            write(path, &bytes)?;
            outputs.push(path.display().to_string());
        }
        if let Some(path) = &cli.stl {
            let options = crate::stl::options(cli.stl_chord, cli.stl_angle);
            let (bytes, s) = crate::stl::document_stl(document, &options)
                .map_err(|detail| ProjectError::Document { detail })?;
            write(path, &bytes)?;
            outputs.push(path.display().to_string());
            lines.push(crate::report::stl_line(&s).trim_end().to_string());
            stl = Some(s);
        }
        lines.push(format!(
            "built: {} migration(s) applied",
            build.record.migrations.len()
        ));
    }
    for output in &outputs {
        lines.push(format!("wrote {output}"));
    }
    let json = json!({
        "project": project.manifest.name,
        "outcome": if build.error.is_none() { "built" } else { "failed" },
        "versions": versions,
        "error": error_value,
        "outputs": outputs,
        "stl": stl,
    });
    Ok(Summary { lines, json, code })
}

fn write(path: &Path, bytes: &[u8]) -> Result<(), ProjectError> {
    std::fs::write(path, bytes).map_err(|e| ProjectError::Io {
        path: path.display().to_string(),
        detail: e.to_string(),
    })
}

/// `popkern build`: print the report and return the exit code.
pub fn build_command(cli: &Cli) -> u8 {
    let summary = run_build(cli).unwrap_or_else(|error| Summary {
        lines: vec![render_error(&error)],
        json: json!({ "outcome": "failed", "error": error_json(&error) }),
        code: error_code(&error),
    });
    if cli.json {
        println!("{}", summary.json);
    } else {
        for line in &summary.lines {
            println!("{line}");
        }
    }
    summary.code
}

/// `popkern new`: scaffold a project and return the exit code.
pub fn new_command(cli: &Cli) -> u8 {
    let dir = &cli.script;
    if let Err(e) = std::fs::create_dir_all(dir) {
        eprintln!("popkern: cannot create '{}': {e}", dir.display());
        return 4;
    }
    let name = cli.name.clone().unwrap_or_else(|| {
        dir.file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| "Project".into())
    });
    let created = NativeFs::new(dir).and_then(|fs| popkern_project::scaffold::create(&fs, &name));
    match created {
        Ok(()) => {
            println!(
                "created {} — build it with: popkern build {}",
                dir.display(),
                dir.display()
            );
            0
        }
        Err(error) => {
            eprintln!("popkern: {error}");
            error_code(&error)
        }
    }
}
