//! The project runner (live-projects-and-history.md §5.4, §6.1): build
//! a project's document by running its migrations in order, each on
//! the document the earlier ones left, keeping an in-memory checkpoint
//! after every migration so an edit reruns only from the first changed
//! migration.

use std::collections::HashMap;
use std::rc::Rc;
use std::time::Duration;

use popkern_assembly::{Assembly, Document};
use popkern_rhai::outcome::Outcome;
use popkern_rhai::{HostOptions, ScriptHost, ScriptRunError, document_first_failure};

use crate::error::ProjectError;
use crate::files::{ProjectFs, ScriptStore};
use crate::project::{BaseFile, Project, sha256_hex};
use crate::record::{FeatureRef, MigrationRecord, ProjectRecord, attribute};

/// How one migration fared in a build.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VersionStatus {
    /// It ran in this build.
    Built,
    /// Its checkpoint was reused (it and everything before it unchanged).
    Reused,
    /// It failed; the build stopped here.
    Failed,
    /// It did not run because an earlier migration failed.
    Stale,
}

impl VersionStatus {
    /// The status as hosts spell it.
    pub fn as_str(self) -> &'static str {
        match self {
            VersionStatus::Built => "built",
            VersionStatus::Reused => "cached",
            VersionStatus::Failed => "failed",
            VersionStatus::Stale => "stale",
        }
    }
}

/// One migration's line in a build report.
#[derive(Debug, Clone)]
pub struct VersionReport {
    /// The migration's name (file stem).
    pub name: String,
    /// Its content hash.
    pub sha256: String,
    /// How it fared.
    pub status: VersionStatus,
    /// Wall time when it ran and a clock was available.
    pub elapsed: Option<Duration>,
    /// What it printed.
    pub printed: Vec<String>,
    /// Features it recorded.
    pub features: Vec<FeatureRef>,
    /// Earlier features it changed.
    pub edited: Vec<FeatureRef>,
}

/// The result of a build.
pub struct ProjectBuild {
    /// One line per migration, in order.
    pub versions: Vec<VersionReport>,
    /// The document after the last migration that succeeded (the base
    /// alone when the first migration failed; `None` when there is
    /// nothing yet).
    pub document: Option<Document>,
    /// The applied-migration record of `document`.
    pub record: ProjectRecord,
    /// Why the build stopped early, when it did.
    pub error: Option<ProjectError>,
}

impl ProjectBuild {
    /// Whether every migration was applied.
    pub fn complete(&self) -> bool {
        self.error.is_none()
    }
}

struct Checkpoint {
    key: String,
    document: Document,
    report: VersionReport,
    record: MigrationRecord,
}

/// Runs a project's migrations with checkpoints (see the module docs).
pub struct ProjectRunner {
    fs: Rc<dyn ProjectFs>,
    options: HostOptions,
    base: Option<(String, Document)>,
    checkpoints: Vec<Checkpoint>,
    /// For each migration name, the content keys of everything before it
    /// when it last built successfully (failure attribution).
    built_after: HashMap<String, Vec<(String, String)>>,
}

/// The checkpoint key after `previous` and a migration with `sha256`.
fn chain(previous: &str, sha256: &str) -> String {
    sha256_hex(format!("{previous}:{sha256}").as_bytes())
}

/// Replay every part of `document` (a loaded `.pkn` holds recipes, not
/// geometry).
fn replay_all(document: &mut Document) -> Result<(), String> {
    for (part, model) in document.models_mut() {
        model.replay().map_err(|e| match part {
            Some(part) => format!("part '{part}': {e}"),
            None => e.to_string(),
        })?;
    }
    Ok(())
}

/// Load the base document.
fn load_base(base: &BaseFile, name: &str) -> Result<Document, ProjectError> {
    let invalid = |detail: String| ProjectError::BaseInvalid {
        path: base.path.clone(),
        detail,
    };
    let lower = base.path.to_ascii_lowercase();
    if lower.ends_with(".pkn") {
        let mut document =
            Document::from_pkn_bytes(&base.bytes).map_err(|e| invalid(e.to_string()))?;
        replay_all(&mut document).map_err(invalid)?;
        Ok(document)
    } else if lower.ends_with(".step") || lower.ends_with(".stp") {
        let (assembly, _) =
            Assembly::from_step_bytes(&base.bytes, name).map_err(|e| invalid(e.to_string()))?;
        Ok(Document::Assembly(Box::new(assembly)))
    } else {
        Err(invalid("a base is a .pkn or .step file".into()))
    }
}

/// A failed feature left in `document` as a script error.
fn feature_failure(document: &mut Document, script: &str) -> Option<ScriptRunError> {
    match document_first_failure(document)? {
        Outcome::KernelDecline {
            kind,
            rendered,
            payload,
        } => Some(ScriptRunError {
            kind,
            message: rendered,
            script: script.to_string(),
            line: None,
            column: None,
            payload: Box::new(payload),
        }),
        other => Some(ScriptRunError {
            kind: "FeatureFailed".into(),
            message: format!("{other:?}"),
            script: script.to_string(),
            line: None,
            column: None,
            payload: Box::new(serde_json::Value::Null),
        }),
    }
}

impl ProjectRunner {
    /// A runner over the project in `fs`, running scripts under `options`.
    pub fn new(fs: Rc<dyn ProjectFs>, options: HostOptions) -> ProjectRunner {
        ProjectRunner {
            fs,
            options,
            base: None,
            checkpoints: Vec::new(),
            built_after: HashMap::new(),
        }
    }

    /// Replace the script budgets for later builds (checkpoints stay).
    pub fn set_options(&mut self, options: HostOptions) {
        self.options = options;
    }

    /// The project store.
    pub fn fs(&self) -> &Rc<dyn ProjectFs> {
        &self.fs
    }

    /// Build `project`: reuse every checkpoint whose key still matches,
    /// run the rest in order, stop at the first failure.
    ///
    /// # Errors
    ///
    /// A base that cannot be loaded. A failing migration is reported in
    /// [`ProjectBuild::error`], not as an `Err`.
    pub fn build(&mut self, project: &Project) -> Result<ProjectBuild, ProjectError> {
        self.build_from(project, 0, None)
    }

    /// Build `project` starting at migration `from` on `start` — a saved
    /// document that already applied `start.1` (migrations `0..from`).
    /// `from = 0` with no start builds from the base.
    ///
    /// # Errors
    ///
    /// As [`ProjectRunner::build`].
    pub fn build_from(
        &mut self,
        project: &Project,
        from: usize,
        start: Option<(Document, Vec<MigrationRecord>)>,
    ) -> Result<ProjectBuild, ProjectError> {
        let base_key = project
            .base
            .as_ref()
            .map_or_else(|| "none".to_string(), |b| b.sha256.clone());
        let mut record = ProjectRecord::new(project);
        let mut versions = Vec::new();
        let mut keys = vec![("base".to_string(), base_key.clone())];
        let mut key = base_key.clone();
        let mut current: Option<Document> = match &start {
            Some((document, applied)) => {
                record.migrations = applied.clone();
                Some(document.clone())
            }
            None => match &project.base {
                None => None,
                Some(base) => {
                    if self
                        .base
                        .as_ref()
                        .is_none_or(|(sha, _)| *sha != base.sha256)
                    {
                        let document = load_base(base, &project.manifest.name)?;
                        self.base = Some((base.sha256.clone(), document));
                    }
                    self.base.as_ref().map(|(_, d)| d.clone())
                }
            },
        };
        for (i, migration) in project.migrations.iter().enumerate() {
            key = chain(&key, &migration.sha256);
            if i < from {
                keys.push((migration.name.clone(), migration.sha256.clone()));
                continue;
            }
            if start.is_none()
                && let Some(checkpoint) = self.checkpoints.get(i).filter(|c| c.key == key)
            {
                let mut report = checkpoint.report.clone();
                report.status = VersionStatus::Reused;
                report.name.clone_from(&migration.name);
                versions.push(report);
                let mut applied = checkpoint.record.clone();
                applied.name.clone_from(&migration.name);
                record.migrations.push(applied);
                current = Some(checkpoint.document.clone());
                keys.push((migration.name.clone(), migration.sha256.clone()));
                continue;
            }
            self.checkpoints.truncate(i);
            let mut host =
                ScriptHost::with_files(self.options, Rc::new(ScriptStore(self.fs.clone())));
            let script = migration.file_name();
            let outcome = host
                .run_migration(&migration.source, &script, current.clone())
                .and_then(
                    |mut out| match feature_failure(&mut out.document, &script) {
                        Some(error) => Err(error),
                        None => Ok(out),
                    },
                );
            match outcome {
                Ok(out) => {
                    let (features, edited) = attribute(current.as_ref(), &out.document);
                    let report = VersionReport {
                        name: migration.name.clone(),
                        sha256: migration.sha256.clone(),
                        status: VersionStatus::Built,
                        elapsed: out.elapsed,
                        printed: out.printed,
                        features: features.clone(),
                        edited: edited.clone(),
                    };
                    let applied = MigrationRecord {
                        name: migration.name.clone(),
                        sha256: migration.sha256.clone(),
                        features,
                        edited,
                    };
                    // Remembered by name and by content, so a renamed
                    // migration keeps its history.
                    self.built_after
                        .insert(migration.name.clone(), keys.clone());
                    self.built_after
                        .insert(migration.sha256.clone(), keys.clone());
                    if start.is_none() {
                        self.checkpoints.push(Checkpoint {
                            key: key.clone(),
                            document: out.document.clone(),
                            report: report.clone(),
                            record: applied.clone(),
                        });
                    }
                    versions.push(report);
                    record.migrations.push(applied);
                    current = Some(out.document);
                    keys.push((migration.name.clone(), migration.sha256.clone()));
                }
                Err(error) => {
                    let upstream_changed =
                        self.upstream_changed(&migration.name, &migration.sha256, &keys);
                    versions.push(VersionReport {
                        name: migration.name.clone(),
                        sha256: migration.sha256.clone(),
                        status: VersionStatus::Failed,
                        elapsed: None,
                        printed: Vec::new(),
                        features: Vec::new(),
                        edited: Vec::new(),
                    });
                    for later in &project.migrations[i + 1..] {
                        versions.push(VersionReport {
                            name: later.name.clone(),
                            sha256: later.sha256.clone(),
                            status: VersionStatus::Stale,
                            elapsed: None,
                            printed: Vec::new(),
                            features: Vec::new(),
                            edited: Vec::new(),
                        });
                    }
                    return Ok(ProjectBuild {
                        versions,
                        document: current,
                        record,
                        error: Some(ProjectError::MigrationFailed {
                            migration: script,
                            error: Box::new(error),
                            upstream_changed,
                        }),
                    });
                }
            }
        }
        Ok(ProjectBuild {
            versions,
            document: current,
            record,
            error: None,
        })
    }

    /// The upstream migrations (and base) whose content differs from
    /// when this migration (by content, else by name) last built
    /// successfully. Upstream steps compare by position and content,
    /// so renaming an upstream file is not a change.
    fn upstream_changed(&self, name: &str, sha256: &str, now: &[(String, String)]) -> Vec<String> {
        let Some(then) = self
            .built_after
            .get(sha256)
            .or_else(|| self.built_after.get(name))
        else {
            return Vec::new();
        };
        now.iter()
            .enumerate()
            .filter(|(k, (_, sha))| then.get(*k).is_none_or(|(_, old)| old != sha))
            .map(|(_, (upstream, _))| upstream.clone())
            .collect()
    }
}
