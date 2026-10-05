//! `popkern.toml`: the project manifest.
//!
//! The manifest is a small, strict subset of TOML — `[section]`
//! headers and `key = value` lines whose value is a basic string, an
//! integer, or a boolean; `#` comments. Every key is known: a typo is
//! an error with its line, not a silently ignored setting. (The subset
//! keeps the kernel bundle free of a TOML dependency; a full TOML file
//! using only these forms reads the same.)
//!
//! ```toml
//! [project]
//! name = "Bracket"
//! format = 1
//! base = "imports/vendor.step"   # optional: a STEP or .pkn file
//!
//! [gallery]                       # examples only
//! title = "Bracket"
//! blurb = "…"
//! demonstrates = "…"
//! thumbnail = "thumbnail.png"
//! ```

use std::collections::BTreeMap;

use crate::error::ProjectError;

/// The project-format version this kernel reads and writes.
pub const PROJECT_FORMAT: u64 = 1;

/// The manifest's file name.
pub const MANIFEST: &str = "popkern.toml";

/// A parsed manifest.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Manifest {
    /// The project's display name.
    pub name: String,
    /// The project-format version.
    pub format: u64,
    /// The base document (a STEP or `.pkn` file), project-relative.
    pub base: Option<String>,
    /// Gallery metadata (example projects).
    pub gallery: Option<Gallery>,
}

/// How an example project presents itself in the gallery.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Gallery {
    /// Title (defaults to the project name).
    pub title: Option<String>,
    /// One-sentence description.
    pub blurb: Option<String>,
    /// What the example demonstrates.
    pub demonstrates: Option<String>,
    /// Thumbnail image, project-relative.
    pub thumbnail: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
enum Value {
    Text(String),
    Integer(i64),
    Bool(bool),
}

impl Manifest {
    /// A manifest for a new project called `name`.
    pub fn new(name: &str) -> Manifest {
        Manifest {
            name: name.to_string(),
            format: PROJECT_FORMAT,
            base: None,
            gallery: None,
        }
    }

    /// Parse manifest `text`.
    ///
    /// # Errors
    ///
    /// [`ProjectError::ManifestInvalid`] with the offending line.
    pub fn parse(text: &str) -> Result<Manifest, ProjectError> {
        let invalid = |line: usize, reason: String| ProjectError::ManifestInvalid {
            path: MANIFEST.into(),
            line,
            reason,
        };
        let mut sections: BTreeMap<String, BTreeMap<String, (usize, Value)>> = BTreeMap::new();
        let mut current: Option<String> = None;
        for (index, raw) in text.lines().enumerate() {
            let number = index + 1;
            let line = strip_comment(raw).trim();
            if line.is_empty() {
                continue;
            }
            if let Some(name) = line.strip_prefix('[').and_then(|l| l.strip_suffix(']')) {
                let name = name.trim();
                if !matches!(name, "project" | "gallery") {
                    return Err(invalid(
                        number,
                        format!("unknown section [{name}] (expected [project] or [gallery])"),
                    ));
                }
                if sections.contains_key(name) {
                    return Err(invalid(number, format!("section [{name}] appears twice")));
                }
                sections.insert(name.to_string(), BTreeMap::new());
                current = Some(name.to_string());
                continue;
            }
            let Some(section) = &current else {
                return Err(invalid(number, "a setting before any [section]".into()));
            };
            let (key, value) = line
                .split_once('=')
                .ok_or_else(|| invalid(number, "expected key = value".into()))?;
            let key = key.trim();
            let value = parse_value(value.trim()).map_err(|reason| invalid(number, reason))?;
            let table = sections.entry(section.clone()).or_default();
            if table.insert(key.to_string(), (number, value)).is_some() {
                return Err(invalid(number, format!("{key} is set twice")));
            }
        }
        let mut project = sections
            .remove("project")
            .ok_or_else(|| invalid(0, "missing the [project] section".into()))?;
        let name = take_text(&mut project, "name")?
            .ok_or_else(|| invalid(0, "[project] needs name = \"…\"".into()))?;
        let format = match project.remove("format") {
            None => PROJECT_FORMAT,
            Some((_, Value::Integer(v))) if v >= 1 && v as u64 <= PROJECT_FORMAT => v as u64,
            Some((line, other)) => {
                return Err(invalid(
                    line,
                    format!(
                        "format = {other:?} is not a project format this kernel reads (1..={PROJECT_FORMAT})"
                    ),
                ));
            }
        };
        let base = take_text(&mut project, "base")?;
        reject_unknown(&project, "project")?;
        let gallery = match sections.remove("gallery") {
            None => None,
            Some(mut table) => {
                let gallery = Gallery {
                    title: take_text(&mut table, "title")?,
                    blurb: take_text(&mut table, "blurb")?,
                    demonstrates: take_text(&mut table, "demonstrates")?,
                    thumbnail: take_text(&mut table, "thumbnail")?,
                };
                reject_unknown(&table, "gallery")?;
                Some(gallery)
            }
        };
        Ok(Manifest {
            name,
            format,
            base,
            gallery,
        })
    }

    /// The manifest as text (what `popkern new` writes).
    pub fn to_text(&self) -> String {
        let mut out = format!(
            "[project]\nname = {}\nformat = {}\n",
            quote(&self.name),
            self.format
        );
        if let Some(base) = &self.base {
            out.push_str(&format!("base = {}\n", quote(base)));
        }
        if let Some(g) = &self.gallery {
            out.push_str("\n[gallery]\n");
            for (key, value) in [
                ("title", &g.title),
                ("blurb", &g.blurb),
                ("demonstrates", &g.demonstrates),
                ("thumbnail", &g.thumbnail),
            ] {
                if let Some(v) = value {
                    out.push_str(&format!("{key} = {}\n", quote(v)));
                }
            }
        }
        out
    }
}

fn take_text(
    table: &mut BTreeMap<String, (usize, Value)>,
    key: &str,
) -> Result<Option<String>, ProjectError> {
    match table.remove(key) {
        None => Ok(None),
        Some((_, Value::Text(s))) => Ok(Some(s)),
        Some((line, _)) => Err(ProjectError::ManifestInvalid {
            path: MANIFEST.into(),
            line,
            reason: format!("{key} must be a \"string\""),
        }),
    }
}

fn reject_unknown(
    table: &BTreeMap<String, (usize, Value)>,
    section: &str,
) -> Result<(), ProjectError> {
    match table.iter().next() {
        None => Ok(()),
        Some((key, (line, _))) => Err(ProjectError::ManifestInvalid {
            path: MANIFEST.into(),
            line: *line,
            reason: format!("unknown setting {key} in [{section}]"),
        }),
    }
}

/// `line` without a trailing `#` comment (a `#` inside a string stays).
fn strip_comment(line: &str) -> &str {
    let mut in_string = false;
    let mut escaped = false;
    for (i, c) in line.char_indices() {
        match c {
            '\\' if in_string && !escaped => {
                escaped = true;
                continue;
            }
            '"' if !escaped => in_string = !in_string,
            '#' if !in_string => return &line[..i],
            _ => {}
        }
        escaped = false;
    }
    line
}

fn parse_value(text: &str) -> Result<Value, String> {
    if let Some(body) = text.strip_prefix('"') {
        let body = body
            .strip_suffix('"')
            .ok_or_else(|| "unterminated string".to_string())?;
        let mut out = String::new();
        let mut chars = body.chars();
        while let Some(c) = chars.next() {
            match c {
                '\\' => match chars.next() {
                    Some('"') => out.push('"'),
                    Some('\\') => out.push('\\'),
                    Some('n') => out.push('\n'),
                    Some('t') => out.push('\t'),
                    other => return Err(format!("unsupported escape \\{}", other.unwrap_or(' '))),
                },
                '"' => return Err("a string may not contain an unescaped \"".into()),
                c => out.push(c),
            }
        }
        return Ok(Value::Text(out));
    }
    match text {
        "true" => return Ok(Value::Bool(true)),
        "false" => return Ok(Value::Bool(false)),
        _ => {}
    }
    text.replace('_', "")
        .parse::<i64>()
        .map(Value::Integer)
        .map_err(|_| format!("{text} is not a \"string\", an integer, or true/false"))
}

fn quote(text: &str) -> String {
    let escaped = text
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\t', "\\t");
    format!("\"{escaped}\"")
}

#[cfg(test)]
#[path = "manifest_tests.rs"]
mod tests;
