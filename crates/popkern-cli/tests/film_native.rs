//! PK-1011 / PK-1012: native film output and still renders
//! (docs/architecture/cli-render-export.md): `popkern film --out` with
//! no browser writes an MP4 with one sample per frame and the film's
//! duration, a GIF with the preset's frame count, a PNG frame at the
//! film's size; a script without a film renders the built-in turntable,
//! the same bytes twice; `popkern render` writes the requested image;
//! what only the browser can draw is refused with exit 4.

use std::io::Cursor;
use std::path::{Path, PathBuf};
use std::process::Command;

fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("film_native")
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

/// 1.2 s at 10 fps: a caption, an orbit, a move, a fade and a section.
const SHORT: &str = r#"let film = storyboard(#{ size: [320, 240], fps: 10, look: "toon" });
film.step("Short", "a test film", 1.2, |s| {
    s.orbit(60.0);
    s.move("high", #{ from: [0.0, 0.0, 2.0] });
    s.fade("low", 1.0, 0.5);
    s.section(#{ from: "z", from_at: 1.0, to_at: 0.6 });
});
film
"#;

fn script_with_film(name: &str) -> PathBuf {
    let dir = scratch(name);
    std::fs::write(dir.join("blocks.rhai"), BLOCKS).expect("script");
    std::fs::write(dir.join("blocks.short.film.rhai"), SHORT).expect("film");
    dir.join("blocks.rhai")
}

fn run(args: &[&str]) -> String {
    let (code, out, err) = popkern(args);
    assert_eq!(code, 0, "{args:?}\nstdout {out}\nstderr {err}");
    out
}

#[test]
fn an_mp4_holds_a_sample_per_frame_for_the_film_duration() {
    let script = script_with_film("mp4");
    let out = script.with_file_name("short.mp4");
    run(&[
        "film",
        script.to_str().expect("utf8"),
        "--out",
        out.to_str().expect("utf8"),
        "--preset",
        "mp4-720p",
    ]);
    let bytes = std::fs::read(&out).expect("mp4");
    assert_eq!(&bytes[4..8], b"ftyp");
    let size = bytes.len() as u64;
    let reader = mp4::Mp4Reader::read_header(Cursor::new(bytes), size).expect("parses");
    let track = reader.tracks().get(&1).expect("video");
    assert_eq!(track.sample_count(), 12, "1.2 s at 10 fps");
    assert_eq!((track.width(), track.height()), (1280, 720));
    assert!((track.duration().as_secs_f64() - 1.2).abs() < 1e-6);
}

#[test]
fn a_gif_and_a_png_frame_follow_their_presets() {
    let script = script_with_film("gif");
    let gif = script.with_file_name("short.gif");
    run(&[
        "film",
        script.to_str().expect("utf8"),
        "--out",
        gif.to_str().expect("utf8"),
        "--preset",
        "gif",
    ]);
    let bytes = std::fs::read(&gif).expect("gif");
    assert!(bytes.starts_with(b"GIF89a"));
    let mut decoder = gif::DecodeOptions::new()
        .read_info(Cursor::new(&bytes))
        .expect("reads");
    assert_eq!(
        (decoder.width(), decoder.height()),
        (640, 480),
        "640 wide, the film's aspect"
    );
    let mut frames = 0;
    while decoder.read_next_frame().expect("frame").is_some() {
        frames += 1;
    }
    assert_eq!(frames, 18, "1.2 s at the preset's 15 fps");
    let png = script.with_file_name("frame.png");
    run(&[
        "film",
        script.to_str().expect("utf8"),
        "--out",
        png.to_str().expect("utf8"),
        "--preset",
        "png",
        "--at",
        "0.6",
    ]);
    let image = image::open(&png).expect("png");
    assert_eq!(
        (image.width(), image.height()),
        (320, 240),
        "the film's size"
    );
}

#[test]
fn a_script_without_a_film_renders_the_turntable_the_same_twice() {
    let dir = scratch("turntable");
    let script = dir.join("blocks.rhai");
    std::fs::write(&script, BLOCKS).expect("script");
    let mut renders = Vec::new();
    for k in 0..2 {
        let out = dir.join(format!("turn{k}.png"));
        let stdout = run(&[
            "film",
            script.to_str().expect("utf8"),
            "--out",
            out.to_str().expect("utf8"),
            "--preset",
            "png",
            "--at",
            "1.0",
        ]);
        assert!(stdout.contains("film turntable"), "{stdout}");
        renders.push(std::fs::read(&out).expect("png"));
    }
    assert_eq!(renders[0], renders[1], "deterministic");
    // Without --out the film list stays as the scripts say.
    let (code, _, err) = popkern(&["film", script.to_str().expect("utf8"), "--json"]);
    assert_eq!(code, 2, "{err}");
}

#[test]
fn render_writes_a_still_at_the_requested_size() {
    let script = script_with_film("render");
    for (name, extra) in [
        ("studio.png", vec![]),
        (
            "toon.png",
            vec![
                "--look",
                "toon",
                "--azimuth",
                "60",
                "--elevation",
                "10",
                "--zoom",
                "1.5",
                "--background",
                "#20242c",
            ],
        ),
        ("frame.png", vec!["--film", "short", "--at", "0.9"]),
    ] {
        let out = script.with_file_name(name);
        let mut args = vec![
            "render",
            script.to_str().expect("utf8"),
            "--out",
            out.to_str().expect("utf8"),
            "--width",
            "300",
            "--height",
            "200",
        ];
        args.extend(extra);
        run(&args);
        let image = image::open(&out).expect("png").to_rgba8();
        assert_eq!(image.dimensions(), (300, 200), "{name}");
        let background = image.get_pixel(0, 0).0;
        let drawn = image.pixels().filter(|p| p.0 != background).count();
        assert!(
            drawn * 33 > 300 * 200,
            "{name}: the model is drawn ({drawn} px)"
        );
    }
    let (code, _, err) = popkern(&["render", script.to_str().expect("utf8")]);
    assert_eq!(code, 2, "--out is required: {err}");
    let (code, _, err) = popkern(&[
        "render",
        script.to_str().expect("utf8"),
        "--out",
        "x.png",
        "--look",
        "neon",
    ]);
    assert_eq!(code, 2, "{err}");
}

#[test]
fn what_only_the_browser_draws_is_refused() {
    let script = script_with_film("webm");
    let out = script.with_file_name("short.webm");
    let (code, _, err) = popkern(&[
        "film",
        script.to_str().expect("utf8"),
        "--out",
        out.to_str().expect("utf8"),
        "--preset",
        "webm",
    ]);
    assert_eq!(code, 4);
    assert!(err.contains("--browser"), "{err}");
}
