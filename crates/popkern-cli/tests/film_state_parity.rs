//! PK-1011: the Rust film state (`popkern_rhai::film::state`) plays
//! every gallery film exactly as `web/live/film-state.js` does: the
//! camera, every item's offset / turns / opacity / visibility /
//! highlight (the items a film touches; untouched ones are left out of
//! the fixtures), joints, explode, the section and its plane, and the
//! caption, at many times including every step boundary, to 1e-9.
//!
//! The JavaScript's answers are golden fixtures under
//! `tests/fixtures/film_state/`. To regenerate after a film or the
//! player changes: `POPKERN_FILM_FIXTURES=1 cargo test -p popkern-cli
//! --test film_state_parity` (writes the inputs and runs
//! `node tools/film-state-fixtures.mjs` on them).

use std::path::{Path, PathBuf};

use popkern_live::{HostOptions, Session};
use popkern_rhai::film::state::{camera_at, section_at, state_at};
use serde_json::{Value, json};

fn examples() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../popkern-rhai/examples")
}

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/film_state")
}

/// The scripts that carry films.
fn filmed_scripts() -> Vec<PathBuf> {
    let mut stems: Vec<String> = std::fs::read_dir(examples())
        .expect("examples")
        .flatten()
        .filter_map(|e| e.file_name().to_str().map(str::to_string))
        .filter_map(|f| f.strip_suffix(".film.rhai").map(str::to_string))
        .map(|f| f.split('.').next().unwrap_or("").to_string())
        .collect();
    stems.sort();
    stems.dedup();
    stems
        .into_iter()
        .map(|s| examples().join(format!("{s}.rhai")))
        .filter(|p| p.exists())
        .collect()
}

/// Selectors and the visible bounds the player captures.
fn scene_inputs(envelope: &popkern_live::Envelope) -> (Vec<String>, [[f64; 3]; 2]) {
    let mut paths = Vec::new();
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    let mut grow = |b: &[[f64; 3]; 2]| {
        for i in 0..3 {
            lo[i] = lo[i].min(b[0][i]);
            hi[i] = hi[i].max(b[1][i]);
        }
    };
    match &envelope.assembly {
        None => {
            for (i, body) in envelope.bodies.iter().enumerate() {
                paths.push(format!("body {i}"));
                grow(&body.bounds);
            }
        }
        Some(assembly) => {
            for instance in assembly.instances.iter().filter(|i| i.is_part) {
                paths.push(instance.path.clone());
            }
            for body in &envelope.bodies {
                grow(&body.bounds);
            }
        }
    }
    (paths, [lo, hi])
}

fn times(film: &Value) -> Vec<f64> {
    let duration = film["duration"].as_f64().unwrap_or(0.0);
    let mut out: Vec<f64> = (0..=((duration * 7.0) as usize))
        .map(|k| k as f64 / 7.0)
        .collect();
    for step in film["steps"].as_array().into_iter().flatten() {
        let start = step["start"].as_f64().unwrap_or(0.0);
        let end = start + step["duration"].as_f64().unwrap_or(0.0);
        for t in [start - 1e-4, start, start + 0.1, end - 0.1, end - 1e-4] {
            out.push(t.max(0.0));
        }
    }
    for track in film["tracks"].as_array().into_iter().flatten() {
        let start = track["start"].as_f64().unwrap_or(0.0);
        let end = start + track["duration"].as_f64().unwrap_or(0.0);
        for t in [start, (start + end) / 2.0, end] {
            out.push(t);
        }
    }
    out.push(duration);
    out.push(duration + 0.5);
    out.sort_by(f64::total_cmp);
    out.dedup();
    out
}

/// `a` and `b` agree: numbers to `tol`, everything else exactly.
fn agree(a: &Value, b: &Value, tol: f64, at: &str) -> Result<(), String> {
    match (a, b) {
        (Value::Number(x), Value::Number(y)) => {
            let (x, y) = (
                x.as_f64().unwrap_or(f64::NAN),
                y.as_f64().unwrap_or(f64::NAN),
            );
            if (x - y).abs() <= tol * (1.0 + x.abs().max(y.abs())) {
                Ok(())
            } else {
                Err(format!("{at}: {x} vs {y}"))
            }
        }
        (Value::Array(x), Value::Array(y)) if x.len() == y.len() => x
            .iter()
            .zip(y)
            .enumerate()
            .try_for_each(|(i, (p, q))| agree(p, q, tol, &format!("{at}[{i}]"))),
        (Value::Object(x), Value::Object(y)) => {
            for key in x.keys().chain(y.keys()) {
                let (p, q) = (
                    x.get(key).unwrap_or(&Value::Null),
                    y.get(key).unwrap_or(&Value::Null),
                );
                agree(p, q, tol, &format!("{at}.{key}"))?;
            }
            Ok(())
        }
        _ if a == b => Ok(()),
        _ => Err(format!("{at}: {a} vs {b}")),
    }
}

#[test]
fn every_gallery_film_plays_as_the_page_plays_it() {
    let regenerate = std::env::var_os("POPKERN_FILM_FIXTURES").is_some();
    std::fs::create_dir_all(fixtures()).expect("fixture dir");
    let mut checked = 0;
    for script in filmed_scripts() {
        let source = popkern_cli::film::load(&script).expect("load");
        let mut session = Session::new(HostOptions::default());
        let envelope = session.run_project(source.files.clone());
        assert!(envelope.outcome.is_served(), "{}: built", script.display());
        let (paths, bounds) = scene_inputs(&envelope);
        for (name, path) in &source.films {
            let text = source
                .files
                .iter()
                .find(|(p, _)| p == path)
                .map(|(_, b)| String::from_utf8_lossy(b).into_owned())
                .expect("film text");
            let reply = session
                .film(&text, path, None)
                .unwrap_or_else(|e| panic!("{name}: {e:?}"));
            let film = reply.film;
            let film_json = serde_json::to_value(&film).expect("film json");
            let stem = script.file_stem().and_then(|s| s.to_str()).unwrap_or("x");
            let file = fixtures().join(format!("{stem}.{name}.json"));
            let times = times(&film_json);
            if regenerate {
                let input =
                    json!({ "film": film_json, "paths": paths, "bounds": bounds, "times": times });
                std::fs::write(&file, serde_json::to_string(&input).expect("json")).expect("write");
                let status = std::process::Command::new("node")
                    .arg(
                        Path::new(env!("CARGO_MANIFEST_DIR"))
                            .join("../../tools/film-state-fixtures.mjs"),
                    )
                    .arg(&file)
                    .status()
                    .expect("node runs");
                assert!(status.success(), "node fixture writer");
            }
            let fixture: Value = serde_json::from_str(
                &std::fs::read_to_string(&file)
                    .unwrap_or_else(|_| panic!("{}: missing; regenerate", file.display())),
            )
            .expect("fixture json");
            agree(&fixture["film"], &film_json, 0.0, "film").unwrap_or_else(|e| {
                panic!("{}: the film changed ({e}); regenerate", file.display())
            });
            assert_eq!(
                fixture["paths"],
                json!(paths),
                "{}: scene paths",
                file.display()
            );
            let states = fixture["states"].as_array().expect("states");
            assert_eq!(states.len(), times.len());
            for (expected, &t) in states.iter().zip(&times) {
                let camera = camera_at(&film, t);
                let state = state_at(&film, t, &paths);
                let plane = state.section.as_ref().and_then(|s| section_at(s, &bounds));
                let ours = json!({
                    "t": t,
                    "camera": camera,
                    "items": state
                        .items
                        .iter()
                        .filter(|(_, i)| {
                            i.offset != [0.0; 3]
                                || !i.turns.is_empty()
                                || i.opacity != 1.0
                                || !i.visible
                                || i.highlight != 0.0
                        })
                        .collect::<Vec<_>>(),
                    "joints": state.joints,
                    "explode": state.explode,
                    "section": state.section,
                    "plane": plane,
                    "caption": state.caption,
                });
                agree(expected, &ours, 1e-9, &format!("{stem}.{name} t={t}"))
                    .unwrap_or_else(|e| panic!("{e}"));
                checked += 1;
            }
        }
    }
    assert!(checked > 500, "only {checked} states compared");
}
