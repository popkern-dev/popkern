//! Project file stores: where a project's files live (a folder on
//! disk, or memory in the browser), always addressed by
//! project-relative path and sandboxed to the project.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::path::{Component, Path, PathBuf};
use std::rc::Rc;

use popkern_rhai::ScriptFiles;
use popkern_script::{ProjectReason, ScriptError};

use crate::error::ProjectError;

/// The files of one project. Every path is project-relative with `/`
/// separators; implementations refuse paths that leave the project.
pub trait ProjectFs {
    /// The bytes of `path`.
    fn read(&self, path: &str) -> Result<Vec<u8>, ProjectError>;
    /// Create or replace `path` (creating its folders).
    fn write(&self, path: &str, bytes: &[u8]) -> Result<(), ProjectError>;
    /// The file names directly in folder `dir` (empty when it does not
    /// exist), sorted.
    fn list(&self, dir: &str) -> Result<Vec<String>, ProjectError>;
    /// Where `path` lives on a filesystem, when the store has one.
    fn native_path(&self, path: &str) -> Result<Option<PathBuf>, ProjectError>;
}

/// `path` normalised to project-relative `/`-separated components:
/// `.` components dropped; absolute paths, drive prefixes and any `..`
/// refused (a sandbox does not reason about where `..` lands).
pub fn normalize(path: &str) -> Result<String, ProjectError> {
    let outside = || ProjectError::PathOutsideProject {
        path: path.to_string(),
    };
    if path.contains('\\') {
        return Err(outside());
    }
    let mut parts = Vec::new();
    for component in Path::new(path).components() {
        match component {
            Component::Normal(part) => parts.push(part.to_str().ok_or_else(outside)?.to_string()),
            Component::CurDir => {}
            Component::ParentDir | Component::RootDir | Component::Prefix(_) => {
                return Err(outside());
            }
        }
    }
    if parts.is_empty() {
        return Err(outside());
    }
    Ok(parts.join("/"))
}

/// A project folder on disk.
#[derive(Debug, Clone)]
pub struct NativeFs {
    root: PathBuf,
}

impl NativeFs {
    /// The project rooted at `root` (which must exist).
    pub fn new(root: impl AsRef<Path>) -> Result<NativeFs, ProjectError> {
        let root = root.as_ref();
        let root = root.canonicalize().map_err(|e| ProjectError::Io {
            path: root.display().to_string(),
            detail: e.to_string(),
        })?;
        Ok(NativeFs { root })
    }

    /// The project folder.
    pub fn root(&self) -> &Path {
        &self.root
    }

    /// The native path of project path `path`, checked to stay inside
    /// the root once symbolic links resolve.
    fn resolve(&self, path: &str) -> Result<PathBuf, ProjectError> {
        let relative = normalize(path)?;
        let full = self.root.join(&relative);
        // The deepest existing ancestor must canonicalise inside the
        // root (a symlinked folder could point anywhere).
        let mut probe = full.as_path();
        while !probe.exists() {
            match probe.parent() {
                Some(parent) => probe = parent,
                None => break,
            }
        }
        let real = probe.canonicalize().map_err(|e| ProjectError::Io {
            path: relative.clone(),
            detail: e.to_string(),
        })?;
        if !real.starts_with(&self.root) {
            return Err(ProjectError::PathOutsideProject {
                path: path.to_string(),
            });
        }
        Ok(full)
    }
}

impl ProjectFs for NativeFs {
    fn read(&self, path: &str) -> Result<Vec<u8>, ProjectError> {
        let full = self.resolve(path)?;
        std::fs::read(&full).map_err(|e| match e.kind() {
            std::io::ErrorKind::NotFound => ProjectError::FileNotFound {
                path: path.to_string(),
            },
            _ => ProjectError::Io {
                path: path.to_string(),
                detail: e.to_string(),
            },
        })
    }

    fn write(&self, path: &str, bytes: &[u8]) -> Result<(), ProjectError> {
        let full = self.resolve(path)?;
        let io = |e: std::io::Error| ProjectError::Io {
            path: path.to_string(),
            detail: e.to_string(),
        };
        if let Some(parent) = full.parent() {
            std::fs::create_dir_all(parent).map_err(io)?;
        }
        std::fs::write(&full, bytes).map_err(io)
    }

    fn list(&self, dir: &str) -> Result<Vec<String>, ProjectError> {
        let full = self.resolve(dir)?;
        if !full.is_dir() {
            return Ok(Vec::new());
        }
        let entries = std::fs::read_dir(&full).map_err(|e| ProjectError::Io {
            path: dir.to_string(),
            detail: e.to_string(),
        })?;
        let mut names: Vec<String> = entries
            .filter_map(Result::ok)
            .filter(|e| e.file_type().is_ok_and(|t| t.is_file()))
            .filter_map(|e| e.file_name().to_str().map(str::to_string))
            .collect();
        names.sort();
        Ok(names)
    }

    fn native_path(&self, path: &str) -> Result<Option<PathBuf>, ProjectError> {
        self.resolve(path).map(Some)
    }
}

/// A project held in memory (the browser's store hands its files over
/// as a map).
#[derive(Debug, Clone, Default)]
pub struct MemoryFs {
    files: Rc<RefCell<BTreeMap<String, Vec<u8>>>>,
}

impl MemoryFs {
    /// A store holding `files` (project-relative path → bytes).
    pub fn new(
        files: impl IntoIterator<Item = (String, Vec<u8>)>,
    ) -> Result<MemoryFs, ProjectError> {
        let mut map = BTreeMap::new();
        for (path, bytes) in files {
            map.insert(normalize(&path)?, bytes);
        }
        Ok(MemoryFs {
            files: Rc::new(RefCell::new(map)),
        })
    }

    /// Replace every file with `other`'s (a host handing over the
    /// project's current state; runners holding this store see it).
    pub fn replace_with(&self, other: &MemoryFs) {
        let files = other.files.borrow().clone();
        *self.files.borrow_mut() = files;
    }

    /// Every file, path → bytes (what a host writes back).
    pub fn files(&self) -> BTreeMap<String, Vec<u8>> {
        self.files.borrow().clone()
    }
}

impl ProjectFs for MemoryFs {
    fn read(&self, path: &str) -> Result<Vec<u8>, ProjectError> {
        let key = normalize(path)?;
        self.files
            .borrow()
            .get(&key)
            .cloned()
            .ok_or_else(|| ProjectError::FileNotFound {
                path: path.to_string(),
            })
    }

    fn write(&self, path: &str, bytes: &[u8]) -> Result<(), ProjectError> {
        let key = normalize(path)?;
        self.files.borrow_mut().insert(key, bytes.to_vec());
        Ok(())
    }

    fn list(&self, dir: &str) -> Result<Vec<String>, ProjectError> {
        let prefix = format!("{}/", normalize(dir)?);
        Ok(self
            .files
            .borrow()
            .keys()
            .filter_map(|k| k.strip_prefix(&prefix))
            .filter(|rest| !rest.contains('/'))
            .map(str::to_string)
            .collect())
    }

    fn native_path(&self, _: &str) -> Result<Option<PathBuf>, ProjectError> {
        Ok(None)
    }
}

/// A project store as the script host's file verbs see it.
pub(crate) struct ScriptStore(pub(crate) Rc<dyn ProjectFs>);

fn script_error(error: ProjectError) -> ScriptError {
    match error {
        ProjectError::PathOutsideProject { path } => ScriptError::Project {
            reason: ProjectReason::PathOutsideProject,
            detail: format!("{path} is outside the project folder"),
            path: Some(path),
        },
        ProjectError::FileNotFound { path } => ScriptError::Project {
            reason: ProjectReason::FileNotFound,
            detail: format!("{path} does not exist in the project"),
            path: Some(path),
        },
        other => ScriptError::Io {
            path: String::new(),
            detail: other.to_string(),
        },
    }
}

impl ScriptFiles for ScriptStore {
    fn read(&self, path: &str) -> Result<Vec<u8>, ScriptError> {
        self.0.read(path).map_err(script_error)
    }

    fn write(&self, path: &str, bytes: &[u8]) -> Result<(), ScriptError> {
        self.0.write(path, bytes).map_err(script_error)
    }

    fn native_path(&self, path: &str) -> Result<Option<PathBuf>, ScriptError> {
        self.0.native_path(path).map_err(script_error)
    }
}

#[cfg(test)]
#[path = "files_tests.rs"]
mod tests;
