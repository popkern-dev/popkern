//! Projects on disk for the page's ServerStore: each subfolder of the
//! root holding a `popkern.toml` is a project, its id the folder name.
//! The page exchanges a project's text files (`popkern.toml` and
//! `migrations/*.rhai`); other files in the folder (imports, a
//! `.popkern/` cache, git) are never touched. A project's revision is
//! the SHA-256 of its exchanged files, so a save names the content it
//! edited and a stale save — another tab, `git pull`, an editor —
//! is refused with the stored project.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use popkern_project::{MANIFEST, MIGRATIONS, sha256_hex};
use serde::Serialize;

/// A project as the page sees it.
#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct StoredProject {
    /// The folder name.
    pub id: String,
    /// The manifest's name.
    pub name: String,
    /// The exchanged files, path -> text.
    pub files: BTreeMap<String, String>,
    /// The content revision.
    pub revision: String,
    /// Last modification time (RFC 3339 UTC).
    pub updated: String,
    /// The project folder on the serving machine.
    pub path: String,
}

/// Why a store operation failed.
#[derive(Debug, thiserror::Error)]
pub enum StoreError {
    /// No such project.
    #[error("no project {0}")]
    NotFound(String),
    /// The project changed since the revision the save named.
    #[error("the project changed since it was opened")]
    Conflict(Box<StoredProject>),
    /// A path is not a project file the page may write.
    #[error(
        "{0} is not a project file (popkern.toml, migrations/NNN_name.rhai or films/name.rhai)"
    )]
    BadPath(String),
    /// Reading or writing failed.
    #[error("{0}")]
    Io(String),
}

fn io(e: std::io::Error) -> StoreError {
    StoreError::Io(e.to_string())
}

/// Folder of a project's films (docs/architecture/presentation-films.md §2).
const FILMS: &str = "films";

/// Whether `path` is a file the page exchanges: the manifest, a
/// migration or a film.
fn exchanged(path: &str) -> bool {
    let script = |f: &str| !f.contains('/') && f.ends_with(".rhai") && !f.starts_with('.');
    path == MANIFEST
        || path.strip_prefix("migrations/").is_some_and(script)
        || path.strip_prefix("films/").is_some_and(script)
}

/// A safe folder name for a project called `name`.
fn slug(name: &str) -> String {
    let s: String = name
        .trim()
        .to_lowercase()
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '-' })
        .collect();
    let s = s.trim_matches('-').to_string();
    if s.is_empty() { "project".into() } else { s }
}

fn valid_id(id: &str) -> bool {
    !id.is_empty()
        && !id.starts_with('.')
        && id
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_')
}

fn rfc3339(time: std::time::SystemTime) -> String {
    let secs = time
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let (days, rem) = (secs / 86_400, secs % 86_400);
    // Civil date from days since 1970-01-01 (Howard Hinnant's algorithm).
    let z = days as i64 + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!(
        "{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z",
        rem / 3600,
        (rem % 3600) / 60,
        rem % 60
    )
}

/// Projects under a root folder.
pub struct ProjectsRoot {
    root: PathBuf,
}

impl ProjectsRoot {
    /// The projects under `root`.
    pub fn new(root: impl Into<PathBuf>) -> ProjectsRoot {
        ProjectsRoot { root: root.into() }
    }

    fn folder(&self, id: &str) -> Result<PathBuf, StoreError> {
        if !valid_id(id) {
            return Err(StoreError::NotFound(id.to_string()));
        }
        let folder = self.root.join(id);
        if !folder.join(MANIFEST).is_file() {
            return Err(StoreError::NotFound(id.to_string()));
        }
        Ok(folder)
    }

    /// Every project, newest first.
    pub fn list(&self) -> Result<Vec<StoredProject>, StoreError> {
        let mut out = Vec::new();
        for entry in std::fs::read_dir(&self.root).map_err(io)? {
            let entry = entry.map_err(io)?;
            let id = entry.file_name().to_string_lossy().into_owned();
            if valid_id(&id) && entry.path().join(MANIFEST).is_file() {
                out.push(self.get(&id)?);
            }
        }
        out.sort_by(|a, b| b.updated.cmp(&a.updated));
        Ok(out)
    }

    /// Project `id`.
    pub fn get(&self, id: &str) -> Result<StoredProject, StoreError> {
        let folder = self.folder(id)?;
        let mut files = BTreeMap::new();
        files.insert(
            MANIFEST.to_string(),
            std::fs::read_to_string(folder.join(MANIFEST)).map_err(io)?,
        );
        let mut newest = std::fs::metadata(folder.join(MANIFEST))
            .and_then(|m| m.modified())
            .map_err(io)?;
        for dir in [MIGRATIONS, FILMS] {
            let Ok(entries) = std::fs::read_dir(folder.join(dir)) else {
                continue;
            };
            for entry in entries.flatten() {
                let name = entry.file_name().to_string_lossy().into_owned();
                let path = format!("{dir}/{name}");
                if exchanged(&path) {
                    files.insert(path, std::fs::read_to_string(entry.path()).map_err(io)?);
                    if let Ok(t) = entry.metadata().and_then(|m| m.modified()) {
                        newest = newest.max(t);
                    }
                }
            }
        }
        let name = popkern_project::Manifest::parse(&files[MANIFEST])
            .map(|m| m.name)
            .unwrap_or_else(|_| id.to_string());
        Ok(StoredProject {
            id: id.to_string(),
            name,
            revision: revision(&files),
            files,
            updated: rfc3339(newest),
            path: folder.display().to_string(),
        })
    }

    /// A new project called `name` holding `files`.
    pub fn create(
        &self,
        name: &str,
        files: &BTreeMap<String, String>,
    ) -> Result<StoredProject, StoreError> {
        let base = slug(name);
        let mut id = base.clone();
        let mut n = 2;
        while self.root.join(&id).exists() {
            id = format!("{base}-{n}");
            n += 1;
        }
        let folder = self.root.join(&id);
        std::fs::create_dir_all(folder.join(MIGRATIONS)).map_err(io)?;
        std::fs::write(folder.join(".gitignore"), ".popkern/\n").map_err(io)?;
        write_files(&folder, files)?;
        self.get(&id)
    }

    /// Replace `id`'s exchanged files, provided its revision is still
    /// `revision`; the updated project, or [`StoreError::Conflict`]
    /// carrying the stored one.
    pub fn save(
        &self,
        id: &str,
        files: &BTreeMap<String, String>,
        revision: &str,
    ) -> Result<StoredProject, StoreError> {
        let current = self.get(id)?;
        if current.revision != revision {
            return Err(StoreError::Conflict(Box::new(current)));
        }
        let folder = self.folder(id)?;
        for path in current.files.keys() {
            if !files.contains_key(path) && path != MANIFEST {
                std::fs::remove_file(folder.join(path)).map_err(io)?;
            }
        }
        write_files(&folder, files)?;
        self.get(id)
    }

    /// Move project `id` aside to `<root>/.trash/` (never deleted
    /// outright from the page).
    pub fn remove(&self, id: &str) -> Result<(), StoreError> {
        let folder = self.folder(id)?;
        let trash = self.root.join(".trash");
        std::fs::create_dir_all(&trash).map_err(io)?;
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0);
        std::fs::rename(folder, trash.join(format!("{id}-{stamp}"))).map_err(io)
    }
}

fn revision(files: &BTreeMap<String, String>) -> String {
    let mut all = Vec::new();
    for (path, text) in files {
        all.extend_from_slice(path.as_bytes());
        all.push(0);
        all.extend_from_slice(text.as_bytes());
        all.push(0);
    }
    sha256_hex(&all)
}

fn write_files(folder: &Path, files: &BTreeMap<String, String>) -> Result<(), StoreError> {
    for path in files.keys() {
        if !exchanged(path) {
            return Err(StoreError::BadPath(path.clone()));
        }
    }
    for (path, text) in files {
        let target = folder.join(path);
        if let Some(parent) = target.parent() {
            std::fs::create_dir_all(parent).map_err(io)?;
        }
        let temp = target.with_extension("tmp-popkern");
        std::fs::write(&temp, text).map_err(io)?;
        std::fs::rename(&temp, &target).map_err(io)?;
    }
    Ok(())
}

#[cfg(test)]
#[path = "store_tests.rs"]
mod tests;
