//! The MCP tools: the projects folder (list, read, create, write and
//! delete migration files) and the native kernel (build a project or a
//! loose script, apply parameter values, export STEP).
//!
//! Every tool reads the folder afresh through [`ProjectsRoot`], the
//! same store the page's `./api/projects` uses, so writes go through
//! the same file rules (only `popkern.toml` and `migrations/*.rhai`)
//! and the same atomic replace. A build runs in a fresh native
//! [`Session`] on the request's thread under a wall-clock budget.

use std::collections::BTreeMap;
use std::path::Path;

use popkern_live::envelope::Envelope;
use popkern_live::session::{RunOptions, Session};
use popkern_project::{MANIFEST, MIGRATIONS, Manifest};
use serde_json::{Map, Value, json};

use crate::Context;
use crate::store::{ProjectsRoot, StoreError, StoredProject};

/// A build stops after this long unless the call asks otherwise.
const DEFAULT_TIMEOUT_S: f64 = 300.0;

/// A tool's answer: text for the agent, the same data structured, and
/// whether the tool failed.
pub(crate) struct ToolOutcome {
    structured: Value,
    is_error: bool,
}

impl ToolOutcome {
    fn ok(structured: Value) -> ToolOutcome {
        ToolOutcome {
            structured,
            is_error: false,
        }
    }

    fn failed(message: impl Into<String>) -> ToolOutcome {
        ToolOutcome {
            structured: json!({ "error": message.into() }),
            is_error: true,
        }
    }

    /// The `tools/call` result.
    pub(crate) fn to_json(&self) -> Value {
        let text = serde_json::to_string_pretty(&self.structured).unwrap_or_default();
        json!({
            "content": [{ "type": "text", "text": text }],
            "structuredContent": self.structured,
            "isError": self.is_error,
        })
    }
}

fn string_schema(description: &str) -> Value {
    json!({ "type": "string", "description": description })
}

fn params_schema() -> Value {
    json!({
        "type": "object",
        "description": "Parameter values to apply after the build, by qualified name (`wall`, `fillet-3.radius`), in the parameter's display unit.",
        "additionalProperties": { "type": "number" },
    })
}

fn timeout_schema() -> Value {
    json!({ "type": "number", "description": "Wall-clock budget in seconds (default 300)." })
}

fn tool(name: &str, description: &str, properties: Value, required: &[&str]) -> Value {
    json!({
        "name": name,
        "description": description,
        "inputSchema": { "type": "object", "properties": properties, "required": required },
    })
}

/// The `tools/list` entries.
pub(crate) fn definitions() -> Vec<Value> {
    let project = || string_schema("The project id (its folder name, from list_projects).");
    vec![
        tool(
            "list_projects",
            "List the projects in the served folder: id, name, folder path, migrations, last change.",
            json!({}),
            &[],
        ),
        tool(
            "read_project",
            "Read a project's files (popkern.toml and every migrations/NNN_name.rhai) and its content revision.",
            json!({ "project": project() }),
            &["project"],
        ),
        tool(
            "create_project",
            "Create a project folder with a manifest and a first migration (migrations/001_base.rhai). Without `source` the first migration is the starter plate.",
            json!({
                "name": string_schema("The project name (the folder name is derived from it)."),
                "source": string_schema("The first migration's Rhai source; it must build and end with the model (`let m = model(); ...; m`)."),
            }),
            &["name"],
        ),
        tool(
            "write_file",
            "Create or replace one file of a project: a migration (`migrations/002_holes.rhai`; a bare `002_holes.rhai` means the same) or `popkern.toml`. Migrations run in file-name order.",
            json!({
                "project": project(),
                "path": string_schema("The file path inside the project."),
                "text": string_schema("The file's full new text."),
            }),
            &["project", "path", "text"],
        ),
        tool(
            "delete_file",
            "Delete one migration file from a project.",
            json!({ "project": project(), "path": string_schema("The migration path.") }),
            &["project", "path"],
        ),
        tool(
            "build_project",
            "Build a project's migrations on the native kernel and report how it went: outcome (served, script_error with line, kernel_decline, gate_failure), each migration's status, feature failures, the parameters, and per body its face count, bounds and whether its mesh is watertight.",
            json!({
                "project": project(),
                "params": params_schema(),
                "timeout_s": timeout_schema(),
                "face_names": { "type": "boolean", "description": "Also list every face's stable name per body." },
            }),
            &["project"],
        ),
        tool(
            "run_script",
            "Run a loose Rhai model script (nothing is saved) and report it the way build_project does. Useful to try an idea before writing it into a project.",
            json!({
                "source": string_schema("The script; it ends with the model."),
                "params": params_schema(),
                "timeout_s": timeout_schema(),
                "face_names": { "type": "boolean", "description": "Also list every face's stable name per body." },
            }),
            &["source"],
        ),
        tool(
            "export_step",
            "Build a project and write its bodies as a STEP file into the project's reserved build folder (`.popkern/exports/<id>.step`); answers the file's path.",
            json!({ "project": project(), "params": params_schema(), "timeout_s": timeout_schema() }),
            &["project"],
        ),
    ]
}

fn text_arg<'a>(arguments: &'a Value, key: &str) -> Result<&'a str, ToolOutcome> {
    arguments
        .get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| ToolOutcome::failed(format!("missing string argument `{key}`")))
}

fn store_failure(error: StoreError) -> ToolOutcome {
    match error {
        StoreError::Conflict(_) => ToolOutcome::failed(
            "the project changed while this write was prepared (the page or an editor saved it): read it again and retry",
        ),
        other => ToolOutcome::failed(other.to_string()),
    }
}

fn summary(project: &StoredProject) -> Value {
    let migrations: Vec<&str> = project
        .files
        .keys()
        .filter(|p| p.starts_with("migrations/"))
        .map(String::as_str)
        .collect();
    json!({
        "id": project.id,
        "name": project.name,
        "path": project.path,
        "migrations": migrations,
        "updated": project.updated,
        "revision": project.revision,
    })
}

/// A migration path as the store names it: `NNN_name.rhai` gains the
/// `migrations/` folder.
fn project_path(path: &str) -> String {
    let path = path.trim_start_matches("./");
    if path == MANIFEST || path.starts_with("migrations/") {
        path.to_string()
    } else {
        format!("{MIGRATIONS}/{path}")
    }
}

/// Run `name` with `arguments`; `None` for an unknown tool.
pub(crate) fn call(name: &str, arguments: &Value, context: &Context) -> Option<ToolOutcome> {
    let store = ProjectsRoot::new(&context.options.root);
    let outcome = match name {
        "list_projects" => match store.list() {
            Ok(list) => Ok(ToolOutcome::ok(json!({
                "root": context.options.root.display().to_string(),
                "projects": list.iter().map(summary).collect::<Vec<_>>(),
            }))),
            Err(e) => Err(store_failure(e)),
        },
        "read_project" => read_project(&store, arguments),
        "create_project" => create_project(&store, arguments),
        "write_file" => edit_file(&store, arguments, true),
        "delete_file" => edit_file(&store, arguments, false),
        "build_project" => build_project(&store, arguments).map(|(_, s)| s),
        "run_script" => run_script(arguments),
        "export_step" => export_step(&store, arguments),
        _ => return None,
    };
    Some(outcome.unwrap_or_else(|failure| failure))
}

fn read_project(store: &ProjectsRoot, arguments: &Value) -> Result<ToolOutcome, ToolOutcome> {
    let project = store
        .get(text_arg(arguments, "project")?)
        .map_err(store_failure)?;
    let mut out = summary(&project);
    out["files"] = json!(project.files);
    Ok(ToolOutcome::ok(out))
}

fn create_project(store: &ProjectsRoot, arguments: &Value) -> Result<ToolOutcome, ToolOutcome> {
    let name = text_arg(arguments, "name")?.trim();
    if name.is_empty() {
        return Err(ToolOutcome::failed("the project needs a name"));
    }
    let source = arguments
        .get("source")
        .and_then(Value::as_str)
        .map(str::to_string)
        .unwrap_or_else(|| popkern_project::scaffold::first_migration_source(name));
    let files = BTreeMap::from([
        (MANIFEST.to_string(), Manifest::new(name).to_text()),
        (
            format!(
                "{MIGRATIONS}/{}",
                popkern_project::scaffold::FIRST_MIGRATION
            ),
            source,
        ),
    ]);
    let project = store.create(name, &files).map_err(store_failure)?;
    Ok(ToolOutcome::ok(summary(&project)))
}

/// Write (`write`) or delete one file, revision-checked against the
/// project as just read.
fn edit_file(
    store: &ProjectsRoot,
    arguments: &Value,
    write: bool,
) -> Result<ToolOutcome, ToolOutcome> {
    let id = text_arg(arguments, "project")?;
    let path = project_path(text_arg(arguments, "path")?);
    let project = store.get(id).map_err(store_failure)?;
    let mut files = project.files.clone();
    if write {
        files.insert(path.clone(), text_arg(arguments, "text")?.to_string());
    } else {
        if path == MANIFEST {
            return Err(ToolOutcome::failed("popkern.toml cannot be deleted"));
        }
        if files.remove(&path).is_none() {
            return Err(ToolOutcome::failed(format!("{id} has no file {path}")));
        }
    }
    let saved = store
        .save(id, &files, &project.revision)
        .map_err(store_failure)?;
    let mut out = summary(&saved);
    out[if write { "written" } else { "deleted" }] = json!(path);
    Ok(ToolOutcome::ok(out))
}

fn session(arguments: &Value) -> Session {
    let seconds = arguments
        .get("timeout_s")
        .and_then(Value::as_f64)
        .filter(|s| s.is_finite() && *s > 0.0)
        .unwrap_or(DEFAULT_TIMEOUT_S);
    let options = RunOptions {
        timeout_ms: Some(seconds * 1000.0),
        ..RunOptions::default()
    };
    Session::new(options.host())
}

/// The project's files as the kernel loads them: the exchanged text
/// files, plus the manifest's base file (a STEP to build on).
fn build_files(project: &StoredProject) -> Result<Vec<(String, Vec<u8>)>, ToolOutcome> {
    let mut files: Vec<(String, Vec<u8>)> = project
        .files
        .iter()
        .map(|(p, t)| (p.clone(), t.as_bytes().to_vec()))
        .collect();
    let base = project
        .files
        .get(MANIFEST)
        .and_then(|text| Manifest::parse(text).ok())
        .and_then(|m| m.base);
    if let Some(base) = base {
        let base = popkern_project::normalize(&base)
            .map_err(|e| ToolOutcome::failed(format!("base file: {e}")))?;
        let bytes = std::fs::read(Path::new(&project.path).join(&base))
            .map_err(|e| ToolOutcome::failed(format!("base file {base}: {e}")))?;
        files.push((base, bytes));
    }
    Ok(files)
}

/// Apply `params` to the session's model; the last envelope, or the
/// first that did not serve.
fn apply_params(
    session: &mut Session,
    envelope: Envelope,
    arguments: &Value,
) -> Result<Envelope, ToolOutcome> {
    let Some(params) = arguments.get("params").and_then(Value::as_object) else {
        return Ok(envelope);
    };
    let mut last = envelope;
    for (name, value) in params {
        let value = value
            .as_f64()
            .ok_or_else(|| ToolOutcome::failed(format!("parameter {name}: a number expected")))?;
        last = session.set_param(name, value);
        if last.to_json().get("outcome").and_then(Value::as_str) != Some("served") {
            break;
        }
    }
    Ok(last)
}

fn build_project(
    store: &ProjectsRoot,
    arguments: &Value,
) -> Result<(Session, ToolOutcome), ToolOutcome> {
    let project = store
        .get(text_arg(arguments, "project")?)
        .map_err(store_failure)?;
    let mut session = session(arguments);
    let envelope = session.run_project(build_files(&project)?);
    let envelope = apply_params(&mut session, envelope, arguments)?;
    let mut out = report(&envelope, arguments);
    out["project"]["id"] = json!(project.id);
    out["project"]["path"] = json!(project.path);
    let served = out["outcome"] == "served";
    Ok((
        session,
        ToolOutcome {
            structured: out,
            is_error: !served,
        },
    ))
}

fn run_script(arguments: &Value) -> Result<ToolOutcome, ToolOutcome> {
    let source = text_arg(arguments, "source")?;
    let mut session = session(arguments);
    let envelope = session.run(source, "script.rhai");
    let envelope = apply_params(&mut session, envelope, arguments)?;
    let out = report(&envelope, arguments);
    let served = out["outcome"] == "served";
    Ok(ToolOutcome {
        structured: out,
        is_error: !served,
    })
}

fn export_step(store: &ProjectsRoot, arguments: &Value) -> Result<ToolOutcome, ToolOutcome> {
    let (mut session, built) = build_project(store, arguments)?;
    if built.is_error {
        return Ok(built);
    }
    let id = text_arg(arguments, "project")?;
    let bytes = session
        .export_step(id)
        .map_err(|o| ToolOutcome::failed(format!("STEP export failed: {}", json!(o))))?;
    let folder = Path::new(built.structured["project"]["path"].as_str().unwrap_or(""))
        .join(".popkern")
        .join("exports");
    let file = folder.join(format!("{id}.step"));
    std::fs::create_dir_all(&folder)
        .and_then(|()| std::fs::write(&file, &bytes))
        .map_err(|e| ToolOutcome::failed(format!("cannot write {}: {e}", file.display())))?;
    Ok(ToolOutcome::ok(json!({
        "path": file.display().to_string(),
        "bytes": bytes.len(),
        "bodies": built.structured["bodies"],
    })))
}

/// The envelope as an agent reads it: the outcome and where it
/// failed, migration rows, failed features, parameters, and bodies
/// without their buffers.
fn report(envelope: &Envelope, arguments: &Value) -> Value {
    let all = envelope.to_json();
    let mut out = Map::new();
    for key in [
        "outcome",
        "kind",
        "rendered",
        "payload",
        "step",
        "detail",
        "printed",
        "elapsed_ms",
        "features",
        "params",
        "project",
    ] {
        if let Some(v) = all.get(key).filter(|v| !v.is_null()) {
            out.insert(key.to_string(), v.clone());
        }
    }
    if let Some(report) = all.get("report").filter(|r| !r.is_null()) {
        out.insert(
            "report".into(),
            json!({
                "is_clean": report.get("is_clean"),
                "failures": report.get("failures"),
                "downgraded_blends": report.get("downgraded_blends"),
            }),
        );
    }
    let face_names = arguments
        .get("face_names")
        .and_then(Value::as_bool)
        .unwrap_or(false);
    let bodies: Vec<Value> = envelope
        .bodies
        .iter()
        .map(|b| {
            let mut body = json!({
                "id": b.id,
                "label": b.label,
                "faces": b.face_count,
                "triangles": b.triangle_count,
                "bounds": b.bounds,
                "watertight": b.open_edges == 0,
            });
            if face_names {
                body["face_names"] = json!(b.face_labels);
            }
            body
        })
        .collect();
    out.insert("bodies".into(), json!(bodies));
    if let Some(assembly) = all.get("assembly").filter(|a| !a.is_null()) {
        let parts: Vec<&Value> = assembly
            .get("instances")
            .and_then(Value::as_array)
            .map(|i| i.iter().filter_map(|x| x.get("path")).collect())
            .unwrap_or_default();
        out.insert(
            "assembly".into(),
            json!({ "name": assembly.get("name"), "instances": parts }),
        );
    }
    Value::Object(out)
}
