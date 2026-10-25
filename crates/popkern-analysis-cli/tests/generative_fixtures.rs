//! Generative design through the CLI (PH-116): a `[generative]` setup
//! optimises, writes the design VTU and reports its history; bad
//! combinations decline.

use std::path::PathBuf;

use popkern_analysis_cli::run::run;

#[path = "support/report.rs"]
mod report;
use report::{assert_valid, repo};

fn run_text(name: &str, text: &str) -> popkern_analysis_cli::run::RunOutput {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("generative");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let setup = dir.join(name);
    std::fs::write(&setup, text).unwrap_or_else(|e| panic!("{e}"));
    run(
        &repo("tests/fixtures/analysis/block_two_bores.step"),
        &setup,
    )
}

fn fixture() -> String {
    std::fs::read_to_string(repo(
        "tests/fixtures/analysis/block_two_bores.generative.toml",
    ))
    .unwrap_or_else(|e| panic!("{e}"))
    .replace("max_iterations = 300", "max_iterations = 30")
        + "\n[options]\ntarget_cell_size = 1.25\n"
}

/// PH-127: without a cell size the cells are a third of the filter
/// radius (a 3.75 mm filter: 1.25 mm cells, and no `FilterUnderResolved`);
/// a filter so fine that the block needs more cells than the budget
/// declines naming the count.
#[test]
fn generative_cells_default_to_a_third_of_the_filter() {
    let base = std::fs::read_to_string(repo(
        "tests/fixtures/analysis/block_two_bores.generative.toml",
    ))
    .unwrap_or_else(|e| panic!("{e}"))
    .replace("max_iterations = 300", "max_iterations = 5");
    let out = run_text(
        "default.toml",
        &base.replace("filter_mm = 1.5", "filter_mm = 3.75"),
    );
    assert_valid(&out.report);
    assert_eq!(out.report["status"], "solved", "{}", out.report);
    let cell = out.report["result"]["cell_size"].as_f64().unwrap_or(0.0);
    assert!((cell - 1.25).abs() < 1e-12, "{cell}");
    assert!(!out.report.to_string().contains("FilterUnderResolved"));
    let out = run_text(
        "fine.toml",
        &base.replace("filter_mm = 1.5", "filter_mm = 0.3"),
    );
    assert_valid(&out.report);
    assert_eq!(out.report["status"], "declined", "{}", out.report);
    assert_eq!(
        out.report["error"]["kind"], "CellBudgetExceeded",
        "{}",
        out.report
    );
}

/// The bored block at 1.25 mm cells: the design meets its volume, stiffens
/// well past the uniform start, keeps its top layer, and the VTU holds
/// one hexahedron per entry with its density.
#[test]
fn a_generative_setup_optimises_and_writes_its_design() {
    let out = run_text("block.toml", &fixture());
    assert_valid(&out.report);
    assert_eq!(out.report["status"], "solved", "{}", out.report);
    let r = &out.report["result"];
    let fraction = r["volume_fraction"].as_f64().unwrap_or(0.0);
    assert!((fraction - 0.3).abs() < 1e-9, "{fraction}");
    let history = r["history"].as_array().cloned().unwrap_or_default();
    assert!(!history.is_empty());
    let first = history[0]["compliance"].as_f64().unwrap_or(0.0);
    let last = r["compliance_nmm"].as_f64().unwrap_or(f64::INFINITY);
    assert!(last < 0.5 * first, "{first} -> {last}");
    assert!(r["kept_fraction"].as_f64().unwrap_or(0.0) > 0.0);
    // PH-125: the block is thicker than the filter's 3 mm member over
    // most of it.
    let text = out.report.to_string();
    assert!(!text.contains("ThinDesignSpace"), "{text}");
    // PH-124: its 1.5 mm filter spans 1.2 cells of 1.25 mm.
    assert!(text.contains("FilterUnderResolved"), "{text}");
    // PH-118, PH-119: crisp by default, and the build direction's
    // unsupported volume reported.
    assert!(r["grey_fraction"].is_f64(), "{r}");
    assert!(r["unsupported_fraction"].is_f64(), "{r}");
    // PH-128: the thinnest members, at least a cell.
    assert!(
        r["member_thickness_mm"].as_f64().unwrap_or(0.0) >= 1.25,
        "{r}"
    );
    assert!(
        history
            .iter()
            .any(|h| h["beta"].as_f64().unwrap_or(0.0) > 0.0)
    );
    let entries = r["entries"].as_u64().unwrap_or(0);
    let vtu = out.vtu.unwrap_or_default();
    assert!(vtu.contains(&format!("NumberOfCells=\"{entries}\"")));
    assert!(vtu.contains("Name=\"density\""));
    assert!(vtu.contains("Name=\"floor\""));
    // PH-117: the design's surface is closed, and the strength preview
    // runs on it as a part with the setup's supports and loads.
    let mesh = &r["design_mesh"];
    assert_eq!(mesh["closed"], true, "{mesh}");
    assert!(mesh["triangles"].as_u64().unwrap_or(0) > 0);
    let check = &r["design_check"];
    assert_eq!(check["status"], "solved", "{check}");
    assert!(check["safety"]["governing"].is_u64(), "{check}");
}

/// Out-of-range settings, a combination with `[modal]`, a kept face
/// that is not defined, and a build direction that is not a signed axis
/// all decline as setup errors.
#[test]
fn bad_generative_setups_decline() {
    let base = fixture();
    for (name, text) in [
        (
            "fraction.toml",
            base.replace("volume_fraction = 0.3", "volume_fraction = 1.5"),
        ),
        (
            "filter.toml",
            base.replace("filter_mm = 1.5", "filter_mm = 0"),
        ),
        ("modal.toml", base.clone() + "\n[modal]\nmodes = 3\n"),
        (
            "keep.toml",
            base.replace("keep = [\"top\"]", "keep = [\"lid\"]"),
        ),
        (
            "build.toml",
            base.replace("build_direction = \"+z\"", "build_direction = \"up\""),
        ),
    ] {
        let out = run_text(name, &text);
        assert_valid(&out.report);
        assert_eq!(out.report["status"], "declined", "{name}");
        assert_eq!(out.report["error"]["kind"], "SetupFile", "{name}");
    }
    // The robust formulation works on a crisp design only.
    let grey_robust = base.replace(
        "build_direction = \"+z\"",
        "build_direction = \"+z\"\ncrisp = false\nrobust = true",
    );
    let out = run_text("robust.toml", &grey_robust);
    assert_valid(&out.report);
    assert_eq!(out.report["status"], "declined");
    assert_eq!(out.report["error"]["kind"], "InvalidGenerative");
}

/// `robust = true` solves the eroded design: the report carries its
/// compliance, above the design's own.
#[test]
fn a_robust_setup_reports_its_eroded_design() {
    let text = fixture().replace(
        "build_direction = \"+z\"",
        "build_direction = \"+z\"\nrobust = true",
    );
    let out = run_text("robust_ok.toml", &text);
    assert_valid(&out.report);
    assert_eq!(out.report["status"], "solved", "{}", out.report);
    let r = &out.report["result"];
    let eroded = r["eroded_compliance_nmm"].as_f64().unwrap_or(0.0);
    let own = r["compliance_nmm"].as_f64().unwrap_or(f64::INFINITY);
    assert!(eroded > own, "{eroded} against {own}");
}
