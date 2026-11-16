//! PK-412: the examples gallery cannot rot. One table drives three
//! things — the native gate on every script (served, no failures, no
//! downgraded blends, watertight display mesh within the angular facet
//! bound, STEP through the OCC oracle when present), the manifest the popkern live page reads
//! (`web/live/examples.json` + the script copies beside it),
//! and the thumbnails. `POPKERN_UPDATE_GALLERY=1` regenerates the
//! manifest, copies and thumbnails; without it the test fails on any
//! drift, so a script edit that changes a feature count is a diff in
//! review.

use std::path::{Path, PathBuf};

use popkern_cli::{Outcome, export_document};
use popkern_live::{HostOptions, Session};
use serde_json::{Value, json};

#[path = "gallery/table.rs"]
mod table;
use table::GALLERY;
#[path = "gallery/film_parity.rs"]
mod film_parity;
#[path = "gallery/run_command.rs"]
mod run_command;

/// Project examples (live-projects-and-history.md §12): folder name in
/// crates/popkern-project/gallery, title, blurb, expected feature count,
/// and what it demonstrates. The page lists them beside the scripts and
/// builds them with `run_project`.
const GALLERY_PROJECTS: &[(&str, &str, &str, usize, &str)] = &[(
    "boss_plate_in_steps",
    "Boss plate in three steps",
    "The boss plate as a project: plate, boss, rounds — each migration starts from the document the earlier ones built.",
    6,
    "projects, previous(), migrations",
)];

fn project_gallery_dir() -> PathBuf {
    root().join("crates/popkern-project/gallery")
}

/// The files of project folder `dir`, project-relative, sorted.
fn project_files(dir: &Path) -> Vec<String> {
    let mut out = vec!["popkern.toml".to_string()];
    let mut migrations: Vec<String> = std::fs::read_dir(dir.join("migrations"))
        .expect("migrations")
        .filter_map(|e| {
            e.ok()?
                .file_name()
                .to_str()
                .map(|n| format!("migrations/{n}"))
        })
        .collect();
    migrations.sort();
    out.extend(migrations);
    // Films (presentation-films.md §2) travel with the project.
    if let Ok(entries) = std::fs::read_dir(dir.join("films")) {
        let mut films: Vec<String> = entries
            .filter_map(|e| e.ok()?.file_name().to_str().map(|n| format!("films/{n}")))
            .filter(|p| p.ends_with(".rhai"))
            .collect();
        films.sort();
        out.extend(films);
    }
    out
}

/// A gallery script's films: `(name, file)` for every
/// `stem.<name>.film.rhai` (or `stem.film.rhai`) beside it.
fn script_films(stem: &str) -> Vec<(String, String)> {
    let mut out: Vec<(String, String)> = std::fs::read_dir(examples_dir())
        .expect("examples")
        .filter_map(|e| e.ok()?.file_name().to_str().map(str::to_string))
        .filter(|f| {
            f.ends_with(".film.rhai")
                && (f == &format!("{stem}.film.rhai") || f.starts_with(&format!("{stem}.")))
        })
        .map(|f| (popkern_rhai::film::film_name(&f), f))
        .collect();
    out.sort();
    out
}

/// Every film plays against the session's built document.
fn assert_films_run(label: &str, session: &mut Session, films: &[(String, String)]) {
    for (name, source) in films {
        if let Err(outcome) = session.film(source, name, None) {
            panic!("[{label}] film {name}: {outcome:?}");
        }
    }
}

/// Build every project example; its manifest entry.
fn project_entries() -> Vec<Value> {
    let mut entries = Vec::new();
    for (dir, title, blurb, features, demonstrates) in GALLERY_PROJECTS {
        let folder = project_gallery_dir().join(dir);
        let files: Vec<(String, Vec<u8>)> = project_files(&folder)
            .into_iter()
            .map(|p| {
                let bytes = std::fs::read(folder.join(&p)).expect("project file");
                (p, bytes)
            })
            .collect();
        let mut session = Session::new(HostOptions::default());
        let envelope = session.run_project(files.clone());
        assert!(
            envelope.outcome.is_served(),
            "[{dir}] {:?}",
            envelope.outcome
        );
        let films: Vec<(String, String)> = files
            .iter()
            .filter(|(p, _)| p.starts_with("films/"))
            .map(|(p, b)| (p.clone(), String::from_utf8_lossy(b).into_owned()))
            .collect();
        assert_films_run(dir, &mut session, &films);
        assert_eq!(envelope.features, *features, "[{dir}] feature count");
        assert!(
            envelope.bodies.iter().all(|b| b.open_edges == 0),
            "[{dir}] watertight"
        );
        entries.push(json!({
            "project": dir,
            "files": files.iter().map(|(p, _)| p).collect::<Vec<_>>(),
            "title": title,
            "blurb": blurb,
            "features": features,
            "demonstrates": demonstrates,
        }));
    }
    entries
}

/// f32 buffer rounding on unit normals: ~1e-7 rad per component.
const FACET_SPREAD_SLACK: f64 = 1.0e-4;

/// Worst per-triangle vertex-normal spread (radians) over the page's
/// display buffers: the angular facet bound the viewer relies on for
/// smooth shading and straight zebra stripes. Pole vertices count —
/// a flipped pole normal reads as a ~180° spread.
fn worst_facet_spread(normals: &[f32], indices: &[u32]) -> f64 {
    let normal = |i: u32| {
        let k = 3 * i as usize;
        [
            f64::from(normals[k]),
            f64::from(normals[k + 1]),
            f64::from(normals[k + 2]),
        ]
    };
    let angle = |a: [f64; 3], b: [f64; 3]| {
        let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let la = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
        let lb = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
        (dot / (la * lb)).clamp(-1.0, 1.0).acos()
    };
    indices
        .chunks(3)
        .map(|t| {
            let n = [normal(t[0]), normal(t[1]), normal(t[2])];
            angle(n[0], n[1])
                .max(angle(n[1], n[2]))
                .max(angle(n[2], n[0]))
        })
        .fold(0.0, f64::max)
}

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

fn examples_dir() -> PathBuf {
    root().join("crates/popkern-rhai/examples")
}

fn page_dir() -> PathBuf {
    root().join("web/live")
}

fn updating() -> bool {
    std::env::var("POPKERN_UPDATE_GALLERY").is_ok_and(|v| v == "1")
}

fn tmp(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("gallery");
    std::fs::create_dir_all(&dir).expect("dir");
    dir.join(name)
}

/// One gallery script through the native gate: served through the
/// live session, clean, every body watertight within the display facet
/// bound, the `run` command's exports (STEP and thumbnail) with the OCC
/// oracle on the STEP when present, then its films played.
fn check_script(stem: &str) {
    let Some(&(_, _, _, features, _)) = GALLERY.iter().find(|row| row.0 == stem) else {
        panic!("[{stem}] not in the gallery table");
    };
    let features = &features;
    let script = examples_dir().join(format!("{stem}.rhai"));
    let source = std::fs::read_to_string(&script).unwrap_or_else(|e| panic!("{stem}: {e}"));
    // Native engine: served, clean, every body watertight.
    let mut session = Session::new(HostOptions::default());
    let envelope = session.run(&source, &format!("{stem}.rhai"));
    assert!(
        envelope.outcome.is_served(),
        "[{stem}] {:?}",
        envelope.outcome
    );
    assert_eq!(envelope.features, *features, "[{stem}] feature count");
    let report = envelope.report.as_ref().expect("report");
    assert_eq!(report["failures"], 0, "[{stem}]");
    assert_eq!(
        report["downgraded_blends"], 0,
        "[{stem}] no script leans on a downgrade"
    );
    assert!(!envelope.bodies.is_empty(), "[{stem}] a body on screen");
    for body in &envelope.bodies {
        assert_eq!(body.open_edges, 0, "[{stem}] watertight mesh");
        let worst = worst_facet_spread(&body.normals, &body.indices);
        assert!(
            worst <= popkern_mesh::DISPLAY_MAX_FACET_ANGLE + FACET_SPREAD_SLACK,
            "[{stem}] display mesh facet spread {:.2}° over the {:.1}° bound",
            worst.to_degrees(),
            popkern_mesh::DISPLAY_MAX_FACET_ANGLE.to_degrees()
        );
    }
    // The run command's exports, STEP + thumbnail, of the document the
    // session built (before a film poses it), oracle on the STEP. The
    // command's own build is the same host (PK-1223: it built every
    // script a second time); `the_run_command_writes_step_and_thumbnail`
    // keeps the command whole.
    let step = tmp(&format!("{stem}.step"));
    let png = if updating() {
        page_dir().join("examples").join(format!("{stem}.png"))
    } else {
        tmp(&format!("{stem}.png"))
    };
    let cli = run_command::cli(&script, &step, &png);
    let document = session
        .document_mut()
        .expect("the session keeps its document");
    let summary = export_document(&cli, document);
    assert!(
        matches!(summary.outcome, Outcome::Served),
        "[{stem}] {:?}",
        summary.outcome
    );
    match run_command::oracle_passes(&step) {
        Some(passed) => assert!(passed, "[{stem}] OCC oracle rejected the STEP"),
        None => eprintln!("[pk412] {stem}: OCC oracle absent"),
    }
    // Its films play as the page plays them (the JavaScript's fixtures).
    let films: Vec<(String, String)> = script_films(stem)
        .into_iter()
        .map(|(name, file)| {
            let text = std::fs::read_to_string(examples_dir().join(file)).expect("film");
            (name, text)
        })
        .collect();
    film_parity::films_play_as_the_page_plays_them(stem, &mut session, &envelope, &films);
}

/// One test per gallery script, so the suite runs them side by side
/// (in one test they ran serially, ~37 min, the suite's longest tail)
/// and a failure names its script.
macro_rules! gallery_scripts {
    ($($stem:ident),* $(,)?) => {
        $(
            #[test]
            fn $stem() {
                check_script(stringify!($stem));
            }
        )*
        /// The stems with a test of their own.
        const TESTED: &[&str] = &[$(stringify!($stem)),*];
    };
}

mod script {
    use super::check_script;

    gallery_scripts!(
        boss_plate,
        revolved_knob,
        shelled_enclosure,
        chamfered_bracket,
        bored_block,
        tapered_fillet,
        pyramid_apex,
        mixed_radius_corner,
        parametric_plate,
        project_box,
        ball_bearing,
        involute_gear,
        centrifugal_impeller,
        polyhedral_dice,
        coil_spring,
        propeller,
        keyed_shaft,
        planetary_gearbox,
        threaded_stud,
        hex_bolt,
        countersunk_joint,
        fastener_kit,
        angled_fitting,
        door_hinge,
        piano_hinge,
        living_hinges,
        hinged_box,
        hinge_assembly,
        hinged_case,
        cycloidal_joint,
        six_axis_arm,
        generative_bracket,
    );

    /// Every filmed example is a gallery script, so its films are
    /// compared, every fixture belongs to one of its films, and the
    /// fixtures hold many states.
    #[test]
    fn every_film_is_a_gallery_scripts_and_has_its_fixture() {
        let mut fixtures = Vec::new();
        for entry in std::fs::read_dir(super::examples_dir())
            .expect("examples")
            .flatten()
        {
            let file = entry.file_name().to_string_lossy().into_owned();
            let Some(rest) = file.strip_suffix(".film.rhai") else {
                continue;
            };
            let stem = rest.split('.').next().unwrap_or(rest);
            assert!(
                TESTED.contains(&stem),
                "{file}: its script is not in the gallery"
            );
            fixtures.push(format!(
                "{stem}.{}.json",
                popkern_rhai::film::film_name(&file)
            ));
        }
        fixtures.sort();
        let mut found: Vec<String> = std::fs::read_dir(super::film_parity::fixtures())
            .expect("fixtures")
            .flatten()
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .collect();
        found.sort();
        assert_eq!(found, fixtures, "film fixtures and gallery films differ");
        let states = super::film_parity::fixture_states();
        assert!(states > 500, "only {states} states in the fixtures");
    }

    #[test]
    fn every_gallery_script_has_its_own_test() {
        let table: Vec<&str> = super::GALLERY.iter().map(|row| row.0).collect();
        assert_eq!(
            TESTED,
            table.as_slice(),
            "gallery_scripts! must list the table's stems in order"
        );
    }
}

/// The manifest the page reads and the copies beside it match the table
/// (entries are the table's own rows: no script runs here).
#[test]
fn the_page_manifest_and_copies_match_the_gallery_table() {
    assert!(GALLERY.len() >= 8, "at least eight gallery scripts");
    let mut entries = Vec::new();
    for (stem, title, blurb, features, demonstrates) in GALLERY {
        let films = script_films(stem);
        let mut entry = json!({
            "file": format!("{stem}.rhai"),
            "title": title,
            "blurb": blurb,
            "features": features,
            "demonstrates": demonstrates,
            "thumbnail": format!("examples/{stem}.png"),
        });
        if !films.is_empty() {
            entry["films"] = films
                .iter()
                .map(|(name, file)| json!({ "name": name, "file": file }))
                .collect();
        }
        entries.push(entry);
    }
    entries.extend(project_entries());
    // The manifest and the script copies the page serves.
    let manifest = json!({ "examples": entries });
    let manifest_text = format!(
        "{}\n",
        serde_json::to_string_pretty(&manifest).expect("json")
    );
    let manifest_path = page_dir().join("examples.json");
    let copies = page_dir().join("examples");
    if updating() {
        std::fs::create_dir_all(&copies).expect("dir");
        std::fs::write(&manifest_path, &manifest_text).expect("write manifest");
        for (stem, ..) in GALLERY {
            std::fs::copy(
                examples_dir().join(format!("{stem}.rhai")),
                copies.join(format!("{stem}.rhai")),
            )
            .expect("copy");
            for (_, file) in script_films(stem) {
                std::fs::copy(examples_dir().join(&file), copies.join(&file)).expect("copy film");
            }
        }
        for (dir, ..) in GALLERY_PROJECTS {
            for path in project_files(&project_gallery_dir().join(dir)) {
                let target = copies.join(dir).join(&path);
                std::fs::create_dir_all(target.parent().expect("parent")).expect("dir");
                std::fs::copy(project_gallery_dir().join(dir).join(&path), target).expect("copy");
            }
        }
        return;
    }
    let committed = std::fs::read_to_string(&manifest_path).unwrap_or_default();
    assert_eq!(
        committed, manifest_text,
        "web/live/examples.json drifted from the gallery table: run with \
         POPKERN_UPDATE_GALLERY=1"
    );
    let committed: Value = serde_json::from_str(&committed).expect("manifest json");
    for (stem, ..) in GALLERY {
        let source = std::fs::read(examples_dir().join(format!("{stem}.rhai"))).expect("read");
        let copy = std::fs::read(copies.join(format!("{stem}.rhai"))).unwrap_or_default();
        assert_eq!(
            source, copy,
            "[{stem}] web/live/examples copy drifted: run with POPKERN_UPDATE_GALLERY=1"
        );
        assert!(
            copies.join(format!("{stem}.png")).exists(),
            "[{stem}] thumbnail missing: run with POPKERN_UPDATE_GALLERY=1"
        );
        for (_, file) in script_films(stem) {
            let source = std::fs::read(examples_dir().join(&file)).expect("read");
            let copy = std::fs::read(copies.join(&file)).unwrap_or_default();
            assert_eq!(
                source, copy,
                "[{file}] web/live/examples copy drifted: run with POPKERN_UPDATE_GALLERY=1"
            );
        }
    }
    for (dir, ..) in GALLERY_PROJECTS {
        for path in project_files(&project_gallery_dir().join(dir)) {
            let source = std::fs::read(project_gallery_dir().join(dir).join(&path)).expect("read");
            let copy = std::fs::read(copies.join(dir).join(&path)).unwrap_or_default();
            assert_eq!(
                source, copy,
                "[{dir}/{path}] web/live/examples copy drifted: run with POPKERN_UPDATE_GALLERY=1"
            );
        }
    }
    assert_eq!(
        committed["examples"].as_array().map(Vec::len),
        Some(GALLERY.len() + GALLERY_PROJECTS.len())
    );
}

/// The page's default model ships with the gallery boss plate's tour
/// (a fresh Scratch's `films/tour.rhai`): the two copies never drift.
#[test]
fn the_default_film_is_the_boss_plate_tour() {
    let gallery = std::fs::read(examples_dir().join("boss_plate.tour.film.rhai")).expect("tour");
    let page = std::fs::read(page_dir().join("default.film.rhai")).expect("default film");
    assert_eq!(
        gallery, page,
        "web/live/default.film.rhai drifted from the boss plate's tour"
    );
}
