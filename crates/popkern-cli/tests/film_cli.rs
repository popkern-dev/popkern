//! `popkern film` (docs/architecture/presentation-films.md §8): a script
//! and its sibling films, a project folder's `films/`, `--json`, `--film`
//! selection, and the exit codes of a bad film and an unknown name.

use std::path::{Path, PathBuf};
use std::process::Command;

fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("film_cli")
        .join(name);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("dir");
    dir
}

fn popkern(args: &[&str]) -> (i32, String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_popkern"))
        .args(args)
        .output()
        .expect("popkern");
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

const BLOCKS: &str = r#"let a = assembly("blocks");
let block = a.part("block", |m| { m.profile(xy()).rect(1.0, 1.0).sketch("s").extrude(blind(1.0)); });
a.place(block, "low");
a.place(block, "high").translate(0.0, 0.0, 2.0);
a
"#;

const SPIN: &str = r#"let film = storyboard(#{ size: [640, 360] });
film.step("Spin", "", 2.0, |s| { s.orbit(90.0); s.move("high", #{ from: [0.0, 0.0, 3.0] }); });
film
"#;

#[test]
fn a_script_and_its_films_print_their_storyboards() {
    let dir = scratch("script");
    std::fs::write(dir.join("blocks.rhai"), BLOCKS).expect("write");
    std::fs::write(dir.join("blocks.spin.film.rhai"), SPIN).expect("write");
    std::fs::write(
        dir.join("blocks.hold.film.rhai"),
        "let film = storyboard(); film.hold(1.5); film",
    )
    .expect("write");
    let script = dir.join("blocks.rhai");
    let (code, out, err) = popkern(&["film", script.to_str().expect("utf8"), "--json"]);
    assert_eq!(code, 0, "{err}");
    let films: serde_json::Value = serde_json::from_str(&out).expect("json");
    let names: Vec<&str> = films
        .as_array()
        .expect("array")
        .iter()
        .map(|f| f["name"].as_str().unwrap_or(""))
        .collect();
    assert_eq!(names, ["hold", "spin"]);
    assert_eq!(films[1]["film"]["size"], serde_json::json!([640, 360]));
    assert_eq!(films[1]["film"]["tracks"][1]["channel"]["target"], "high");
    let (code, out, _) = popkern(&["film", script.to_str().expect("utf8"), "--film", "spin"]);
    assert_eq!(code, 0);
    assert!(
        out.contains("film spin: 1 step(s), 2.0 s at 30 fps, 640x360"),
        "{out}"
    );
    let (code, _, err) = popkern(&["film", script.to_str().expect("utf8"), "--film", "nope"]);
    assert_eq!(code, 2);
    assert!(
        err.contains("no film 'nope'; the films are hold, spin"),
        "{err}"
    );
    // --out renders one film: two need --film.
    let (code, _, err) = popkern(&["film", script.to_str().expect("utf8"), "--out", "x.mp4"]);
    assert_eq!(code, 2);
    assert!(err.contains("--out renders one film"), "{err}");
}

#[test]
fn a_project_films_folder_and_a_bad_film() {
    let dir = scratch("project").join("blocks");
    std::fs::create_dir_all(dir.join("migrations")).expect("dir");
    std::fs::create_dir_all(dir.join("films")).expect("dir");
    std::fs::write(dir.join("popkern.toml"), "[project]\nname = \"blocks\"\n").expect("write");
    std::fs::write(dir.join("migrations/001_blocks.rhai"), BLOCKS).expect("write");
    std::fs::write(dir.join("films/spin.rhai"), SPIN).expect("write");
    let (code, out, err) = popkern(&["film", dir.to_str().expect("utf8")]);
    assert_eq!(code, 0, "{err}");
    assert!(out.contains("film spin"), "{out}");
    std::fs::write(
        dir.join("films/broken.rhai"),
        "let film = storyboard();\nfilm.step(\"a\", 1.0, |s| { s.fade(\"lid\", 1.0, 0.0); });\nfilm\n",
    )
    .expect("write");
    let (code, _, err) = popkern(&["film", dir.to_str().expect("utf8"), "--film", "broken"]);
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("no 'lid'") && err.contains("low"), "{err}");
}
