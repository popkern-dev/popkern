//! `popkern film` (docs/architecture/presentation-films.md §8): build a
//! project or a script natively, run its films against the built
//! document (validation and `--json` need no browser), and render one to
//! a file with headless Chrome driving the page's own export.

use std::path::Path;

use popkern_live::{HostOptions, Session};
use popkern_rhai::film::film_name;
use serde_json::{Value, json};

use crate::args::Cli;

/// The presets the page's export menu offers (film-export.js).
pub const PRESETS: &[&str] = &[
    "mp4-1080p",
    "mp4-vertical",
    "mp4-720p",
    "webm",
    "gif",
    "png",
];

/// What `popkern film` works on.
#[derive(Debug, Clone)]
pub struct FilmSource {
    /// Project files (path → bytes), a single script as a one-migration
    /// project.
    pub files: Vec<(String, Vec<u8>)>,
    /// The project's id for the page: its folder name.
    pub id: String,
    /// Films: `(name, path in files)`.
    pub films: Vec<(String, String)>,
}

fn read_dir_files(
    dir: &Path,
    prefix: &str,
    out: &mut Vec<(String, Vec<u8>)>,
) -> std::io::Result<()> {
    let Ok(entries) = std::fs::read_dir(dir.join(prefix)) else {
        return Ok(());
    };
    let mut names: Vec<String> = entries
        .flatten()
        .filter_map(|e| e.file_name().to_str().map(str::to_string))
        .collect();
    names.sort();
    for name in names {
        let path = dir.join(prefix).join(&name);
        if path.is_file() {
            out.push((format!("{prefix}/{name}"), std::fs::read(path)?));
        }
    }
    Ok(())
}

/// Read a project folder or a script and its films.
pub fn load(path: &Path) -> Result<FilmSource, String> {
    if path.is_dir() {
        let mut files = vec![(
            "popkern.toml".to_string(),
            std::fs::read(path.join("popkern.toml"))
                .map_err(|e| format!("{}: {e}", path.join("popkern.toml").display()))?,
        )];
        for dir in ["migrations", "films", "imports"] {
            read_dir_files(path, dir, &mut files).map_err(|e| e.to_string())?;
        }
        let films = files
            .iter()
            .filter(|(p, _)| p.starts_with("films/") && p.ends_with(".rhai"))
            .map(|(p, _)| (film_name(p), p.clone()))
            .collect();
        let canonical = std::fs::canonicalize(path).map_err(|e| e.to_string())?;
        let id = canonical
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("project")
            .to_string();
        return Ok(FilmSource { files, id, films });
    }
    let source = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let file = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("script.rhai");
    let stem = file.strip_suffix(".rhai").unwrap_or(file).to_string();
    let dir = path.parent().unwrap_or(Path::new("."));
    let mut films = Vec::new();
    let mut files = vec![
        (
            "popkern.toml".to_string(),
            popkern_project::Manifest::new(&stem).to_text().into_bytes(),
        ),
        (
            format!(
                "migrations/001_{}.rhai",
                stem.replace(
                    |c: char| !c.is_ascii_alphanumeric() && c != '_' && c != '-',
                    "_"
                )
            ),
            source,
        ),
    ];
    let mut siblings: Vec<String> = std::fs::read_dir(if dir.as_os_str().is_empty() {
        Path::new(".")
    } else {
        dir
    })
    .map_err(|e| e.to_string())?
    .flatten()
    .filter_map(|e| e.file_name().to_str().map(str::to_string))
    .filter(|f| {
        f.ends_with(".film.rhai")
            && (f == &format!("{stem}.film.rhai") || f.starts_with(&format!("{stem}.")))
    })
    .collect();
    siblings.sort();
    for sibling in siblings {
        let name = film_name(&sibling);
        let bytes = std::fs::read(dir.join(&sibling)).map_err(|e| e.to_string())?;
        let path = format!("films/{name}.rhai");
        files.push((path.clone(), bytes));
        films.push((name, path));
    }
    Ok(FilmSource {
        files,
        id: stem.replace(
            |c: char| !c.is_ascii_alphanumeric() && c != '_' && c != '-',
            "-",
        ),
        films,
    })
}

/// Build the document and run every selected film natively: the
/// storyboards, or the exit code and message of the first failure.
pub fn run_films(cli: &Cli, source: &FilmSource) -> Result<Vec<Value>, (u8, String)> {
    let mut session = Session::new(HostOptions {
        max_operations: cli.max_ops,
        timeout: cli.timeout,
        ..HostOptions::default()
    });
    let envelope = session.run_project(source.files.clone());
    if !envelope.outcome.is_served() {
        let code = match &envelope.outcome {
            popkern_rhai::Outcome::ScriptError { .. } => 2,
            popkern_rhai::Outcome::GateFailure { .. } => 4,
            _ => 3,
        };
        return Err((
            code,
            format!(
                "the document did not build: {}",
                serde_json::to_string(&envelope.outcome).unwrap_or_default()
            ),
        ));
    }
    let selected: Vec<&(String, String)> = source
        .films
        .iter()
        .filter(|(name, _)| cli.film.as_deref().is_none_or(|f| f == name))
        .collect();
    if selected.is_empty() {
        let known: Vec<&str> = source.films.iter().map(|(n, _)| n.as_str()).collect();
        return Err((
            2,
            match &cli.film {
                Some(f) => format!(
                    "no film '{f}'; the films are {}",
                    if known.is_empty() {
                        "none".to_string()
                    } else {
                        known.join(", ")
                    }
                ),
                None => {
                    "no films: add films/<name>.rhai (or stem.<name>.film.rhai beside a script)"
                        .to_string()
                }
            },
        ));
    }
    let mut out = Vec::new();
    for (name, path) in selected {
        let text = source
            .files
            .iter()
            .find(|(p, _)| p == path)
            .map(|(_, b)| String::from_utf8_lossy(b).into_owned())
            .unwrap_or_default();
        match session.film(&text, path, None) {
            Ok(reply) => out.push(
                json!({ "name": name, "path": path, "film": reply.film, "printed": reply.printed }),
            ),
            Err(outcome) => {
                return Err((
                    2,
                    format!(
                        "film {name}: {}",
                        serde_json::to_string(&outcome).unwrap_or_default()
                    ),
                ));
            }
        }
    }
    Ok(out)
}

/// `popkern film`: exit code.
pub fn film_command(cli: &Cli) -> u8 {
    let source = match load(&cli.script) {
        Ok(source) => source,
        Err(message) => {
            eprintln!("popkern film: {message}");
            return 5;
        }
    };
    let films = match run_films(cli, &source) {
        Ok(films) => films,
        Err((code, message)) => {
            eprintln!("popkern film: {message}");
            return code;
        }
    };
    if cli.json {
        println!(
            "{}",
            serde_json::to_string_pretty(&films).unwrap_or_default()
        );
    } else {
        for f in &films {
            let film = &f["film"];
            println!(
                "film {}: {} step(s), {:.1} s at {} fps, {}x{}",
                f["name"].as_str().unwrap_or(""),
                film["steps"].as_array().map_or(0, Vec::len),
                film["duration"].as_f64().unwrap_or(0.0),
                film["fps"],
                film["size"][0],
                film["size"][1],
            );
        }
    }
    let Some(out) = &cli.out else {
        return 0;
    };
    if films.len() != 1 {
        eprintln!(
            "popkern film: --out renders one film; pick it with --film (films: {})",
            source
                .films
                .iter()
                .map(|(n, _)| n.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        );
        return 2;
    }
    if !PRESETS.contains(&cli.preset.as_str()) {
        eprintln!("popkern film: --preset is one of {}", PRESETS.join(", "));
        return 2;
    }
    let path = films[0]["path"].as_str().unwrap_or_default().to_string();
    match crate::film_browser::render(cli, &source, &path) {
        Ok(bytes) => match std::fs::write(out, &bytes) {
            Ok(()) => {
                println!(
                    "wrote {} ({:.1} MB)",
                    out.display(),
                    bytes.len() as f64 / 1e6
                );
                0
            }
            Err(e) => {
                eprintln!("popkern film: {}: {e}", out.display());
                4
            }
        },
        Err(message) => {
            eprintln!("popkern film: {message}");
            4
        }
    }
}
