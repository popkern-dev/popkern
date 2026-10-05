//! A project as read from its store: the manifest, the optional base
//! document, and the migrations in order, each with the content hash
//! that keys its checkpoint and its applied-migration record.

use sha2::{Digest, Sha256};

use crate::error::ProjectError;
use crate::files::ProjectFs;
use crate::manifest::{MANIFEST, Manifest};

/// The folder holding migration scripts.
pub const MIGRATIONS: &str = "migrations";

/// One migration script.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Migration {
    /// The ordering number (`010` in `010_pockets.rhai`).
    pub number: u64,
    /// The file stem (`010_pockets`), the migration's name.
    pub name: String,
    /// The project-relative path.
    pub path: String,
    /// The source, line endings normalised to `\n`, BOM removed.
    pub source: String,
    /// SHA-256 of the normalised source, hex.
    pub sha256: String,
}

impl Migration {
    /// The script name diagnostics carry (`010_pockets.rhai`).
    pub fn file_name(&self) -> String {
        format!("{}.rhai", self.name)
    }
}

/// The project's base document file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseFile {
    /// The project-relative path.
    pub path: String,
    /// The file's bytes.
    pub bytes: Vec<u8>,
    /// SHA-256 of the bytes, hex.
    pub sha256: String,
}

/// A loaded project.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Project {
    /// The manifest.
    pub manifest: Manifest,
    /// The base document file, when the manifest names one.
    pub base: Option<BaseFile>,
    /// The migrations in order.
    pub migrations: Vec<Migration>,
}

/// SHA-256 of `bytes` as lowercase hex.
pub fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// `source` with a leading BOM removed and `\r\n` / `\r` line endings
/// as `\n`: the form a migration's hash is taken over, so a checkout
/// with different line endings is the same migration.
pub fn normalize_source(source: &str) -> String {
    source
        .strip_prefix('\u{feff}')
        .unwrap_or(source)
        .replace("\r\n", "\n")
        .replace('\r', "\n")
}

/// The number and name of a migration file name, when it is one.
fn parse_migration_name(file: &str) -> Option<(u64, String)> {
    let stem = file.strip_suffix(".rhai")?;
    let (digits, rest) = stem.split_once('_')?;
    if digits.is_empty() || !digits.chars().all(|c| c.is_ascii_digit()) || rest.is_empty() {
        return None;
    }
    if !rest
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
    {
        return None;
    }
    Some((digits.parse().ok()?, stem.to_string()))
}

impl Project {
    /// Read the project in `fs`.
    ///
    /// # Errors
    ///
    /// A missing or malformed manifest, a misnamed or duplicate-numbered
    /// migration, an unreadable base, or a project with nothing to build.
    pub fn load(fs: &dyn ProjectFs) -> Result<Project, ProjectError> {
        let text = fs.read(MANIFEST)?;
        let text = String::from_utf8(text).map_err(|_| ProjectError::ManifestInvalid {
            path: MANIFEST.into(),
            line: 0,
            reason: "not UTF-8 text".into(),
        })?;
        let manifest = Manifest::parse(&text)?;
        let base = match &manifest.base {
            None => None,
            Some(path) => {
                let bytes = fs.read(path)?;
                Some(BaseFile {
                    path: crate::files::normalize(path)?,
                    sha256: sha256_hex(&bytes),
                    bytes,
                })
            }
        };
        let mut migrations: Vec<Migration> = Vec::new();
        for file in fs.list(MIGRATIONS)? {
            if file.starts_with('.') {
                continue;
            }
            let (number, name) = parse_migration_name(&file)
                .ok_or_else(|| ProjectError::MigrationName { file: file.clone() })?;
            let path = format!("{MIGRATIONS}/{file}");
            let bytes = fs.read(&path)?;
            let text = String::from_utf8(bytes).map_err(|_| ProjectError::Io {
                path: path.clone(),
                detail: "not UTF-8 text".into(),
            })?;
            let source = normalize_source(&text);
            migrations.push(Migration {
                number,
                name,
                path,
                sha256: sha256_hex(source.as_bytes()),
                source,
            });
        }
        migrations.sort_by_key(|m| m.number);
        for pair in migrations.windows(2) {
            if pair[0].number == pair[1].number {
                return Err(ProjectError::MigrationOrder {
                    first: pair[0].file_name(),
                    second: pair[1].file_name(),
                    number: pair[0].number,
                });
            }
        }
        if migrations.is_empty() && base.is_none() {
            return Err(ProjectError::Empty);
        }
        Ok(Project {
            manifest,
            base,
            migrations,
        })
    }

    /// The migration called `name`.
    pub fn migration(&self, name: &str) -> Option<&Migration> {
        self.migrations.iter().find(|m| m.name == name)
    }
}
