//! PK-503: `popkern run --zebra[=mesh|surface]` renders the PNG as a
//! reflection-line zebra view and reports the mesh-vs-surface normal
//! disagreement of the same frame, in text and JSON.

use std::path::{Path, PathBuf};

use popkern_cli::{Cli, Outcome, execute, parse_args, report};
use popkern_viz::zebra::ZebraNormals;

fn boss_plate() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../popkern-rhai/examples/boss_plate.rhai")
        .canonicalize()
        .expect("boss plate example")
}

fn tmp(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("popkern_cli_zebra");
    std::fs::create_dir_all(&dir).expect("dir");
    dir.join(name)
}

fn parsed(args: &[&str]) -> Cli {
    let owned: Vec<String> = args.iter().map(|s| s.to_string()).collect();
    parse_args(&owned).unwrap_or_else(|e| panic!("{e}"))
}

#[test]
fn zebra_flags_parse() {
    let cli = parsed(&["run", "a.rhai"]);
    assert_eq!(cli.zebra, None);
    assert_eq!(cli.stripes, 16);
    assert!((cli.facet_angle_deg - 12.0).abs() < 1e-9);
    assert_eq!(
        parsed(&["run", "a.rhai", "--zebra"]).zebra,
        Some(ZebraNormals::Mesh)
    );
    assert_eq!(
        parsed(&["run", "a.rhai", "--zebra=mesh"]).zebra,
        Some(ZebraNormals::Mesh)
    );
    let cli = parsed(&[
        "run",
        "a.rhai",
        "--zebra=surface",
        "--stripes",
        "10",
        "--facet-angle",
        "6",
    ]);
    assert_eq!(cli.zebra, Some(ZebraNormals::Surface));
    assert_eq!(cli.stripes, 10);
    assert!((cli.facet_angle_deg - 6.0).abs() < 1e-12);
    assert!(parse_args(&["run".into(), "a.rhai".into(), "--zebra=sideways".into()]).is_err());
    assert!(parse_args(&["run".into(), "a.rhai".into(), "--stripes".into()]).is_err());
}

#[test]
fn zebra_render_reports_the_mode_disagreement() {
    let png = tmp("boss_plate_zebra.png");
    let script = boss_plate();
    let cli = parsed(&[
        "run",
        script.to_str().expect("utf8"),
        "--png",
        png.to_str().expect("utf8"),
        "--zebra=surface",
        "--width",
        "320",
        "--height",
        "240",
    ]);
    let summary = execute(&cli);
    assert!(matches!(summary.outcome, Outcome::Served), "{summary:?}");
    assert!(png.exists());
    let z = summary.zebra.as_ref().expect("zebra summary");
    assert_eq!(z.normals, "surface");
    assert!(z.covered_pixels > 1000);
    assert_eq!(z.evaluated_pixels + z.reassigned_pixels, z.covered_pixels);
    assert!(z.reassigned_pixels < z.covered_pixels / 20, "{z:?}");
    assert!(z.ray_hit_pixels as f64 >= 0.99 * z.evaluated_pixels as f64);
    // The boss plate's rounds under the 12° display bound: the facet
    // error is bounded by about half the bound.
    assert!(z.max_angle_deg > 0.0 && z.max_angle_deg < 12.0, "{z:?}");
    let text = report::text(&summary);
    assert!(
        text.contains("zebra (surface normals, 16 stripes, 12° facets)"),
        "{text}"
    );
    let json = serde_json::to_value(&summary).expect("json");
    assert!(json.pointer("/zebra/max_angle_deg").is_some());
    // A plain render carries no zebra section.
    let plain = parsed(&[
        "run",
        script.to_str().expect("utf8"),
        "--png",
        tmp("boss_plate_plain.png").to_str().expect("utf8"),
    ]);
    let summary = execute(&plain);
    assert!(summary.zebra.is_none());
    let json = serde_json::to_value(&summary).expect("json");
    assert!(json.get("zebra").is_none());
}
