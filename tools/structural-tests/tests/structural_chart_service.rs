//! Structural test: chart artifacts are handled only by the chart service
//! (PK-1228, design docs/architecture/chart-service.md §7).
//!
//! A RATCHET over every non-test source file outside
//! `crates/popkern-core/src/geometry/chart/`: occurrences of the patterns
//! that do a chart's period or singular-row handling by hand (periodicity
//! queries, period arithmetic, singular-row tests, the moved wrap) are
//! counted per file against `chart_service_baseline.txt`. A file fails when
//! its count rises above its entry or it is absent and has any. Migrating a
//! consumer lowers its entry (delete it at zero); W2 exits with the
//! baseline empty.

use std::collections::BTreeMap;
use std::path::Path;

fn workspace_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

/// The textual patterns of hand-done chart handling (design §7).
const PATTERNS: &[&str] = &[
    ".parameter_periodicity()",
    ".periodicity()",
    "geometric_chart_periods",
    "marching_period",
    "has_singular_row",
    ".u_period()",
    ".v_period()",
    "evaluation_wrap",
    "wrap_into_domain",
    "rem_euclid(TAU)",
    "rem_euclid(std::f64::consts::TAU)",
    "rem_euclid(period",
    "rem_euclid(p)",
];

/// Occurrences of the patterns not reached through the chart service's
/// own API (`ChartDescriptor::has_singular_row`, `...::evaluation_wrap_of`).
fn count(text: &str) -> usize {
    PATTERNS
        .iter()
        .map(|p| {
            text.match_indices(p)
                .filter(|(i, _)| !text[..*i].ends_with("ChartDescriptor::"))
                .count()
        })
        .sum()
}

/// The chart service: its module, and the surface representation layer
/// that defines periodicity and the in-domain wrap (the `SurfaceEval`
/// trait and `periodicity`, their impls, the NURBS periodic form, and the
/// old paths' re-exports), which must live with the surface types
/// (chart-service.md §7).
const SERVICE: &[&str] = &[
    "crates/popkern-core/src/geometry/chart/",
    "crates/popkern-core/src/geometry/surface.rs",
    "crates/popkern-core/src/geometry/analytic_surfaces.rs",
    "crates/popkern-core/src/geometry/nurbs_surface/mod.rs",
    "crates/popkern-core/src/geometry/nurbs_surface/periodic_form.rs",
    "crates/popkern-core/src/geometry/nurbs_surface/eval_derivs.rs",
    "crates/popkern-core/src/geometry/mod.rs",
    "crates/popkern-math/src/periodic.rs",
];

fn baseline() -> BTreeMap<String, usize> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("chart_service_baseline.txt");
    let text = std::fs::read_to_string(&path).unwrap_or_default();
    text.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| {
            let (file, n) = l.rsplit_once(' ').expect("`path count`");
            (file.to_string(), n.parse().expect("count"))
        })
        .collect()
}

fn counts() -> BTreeMap<String, usize> {
    let root = workspace_root();
    let mut out = BTreeMap::new();
    for file in popkern_lints::find_rs_files(&root.join("crates")) {
        let rel = file
            .strip_prefix(&root)
            .expect("under root")
            .to_string_lossy()
            .replace('\\', "/");
        let test_file =
            rel.contains("/tests/") || rel.ends_with("_tests.rs") || rel.ends_with("/tests.rs");
        // The chart service, and popkern-math's named period arithmetic
        // (`periodic::wrap_from`, `offset_mod`) that it and curve-angle
        // code call: the homes the patterns move to, not consumers.
        if test_file || SERVICE.iter().any(|s| rel.starts_with(s)) {
            continue;
        }
        let n = count(&std::fs::read_to_string(&file).unwrap_or_default());
        if n > 0 {
            out.insert(rel, n);
        }
    }
    out
}

#[test]
fn chart_artifacts_are_handled_only_by_the_chart_service() {
    let base = baseline();
    let now = counts();
    let mut problems = Vec::new();
    for (file, n) in &now {
        match base.get(file) {
            Some(b) if n <= b => {}
            Some(b) => problems.push(format!(
                "{file}: {n} chart-artifact patterns, baseline {b}: use the chart service"
            )),
            None => problems.push(format!(
                "{file}: {n} chart-artifact patterns, none allowed: use the chart service"
            )),
        }
    }
    for (file, b) in &base {
        let n = now.get(file).copied().unwrap_or(0);
        if n < *b {
            problems.push(format!(
                "{file}: down to {n} from {b}: lower its baseline entry (delete it at 0)"
            ));
        }
    }
    assert!(
        problems.is_empty(),
        "{} chart-service ratchet problems:\n{}",
        problems.len(),
        problems.join("\n")
    );
}

/// Writes the current counts as the baseline (run once at W2a):
/// `CHART_RATCHET_WRITE=1 cargo test -p popkern-structural-tests --test structural_chart_service`.
#[test]
fn write_baseline_when_asked() {
    if std::env::var("CHART_RATCHET_WRITE").is_err() {
        return;
    }
    let mut text = String::from(
        "# Chart-service RATCHET baseline (PK-1228, chart-service.md §7): per-file\n\
         # counts of hand-done chart handling outside the chart service. Only\n\
         # ever lowered; W2 exits with this file empty.\n",
    );
    let total: usize = counts().values().sum();
    for (file, n) in counts() {
        text.push_str(&format!("{file} {n}\n"));
    }
    text.push_str(&format!("# total {total}\n"));
    std::fs::write(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("chart_service_baseline.txt"),
        text,
    )
    .expect("write");
}

/// Pole-specific logic outside the service: the `pole` token (any case,
/// any identifier containing it) in non-comment code. A pole is always a
/// chart artifact (topology has none), so every such token is a consumer
/// special-casing what the service owns; the period patterns above cannot
/// see partial-ratio tests or pole gates.
fn pole_count(text: &str) -> usize {
    code_lines(text)
        .map(|l| l.to_ascii_lowercase().matches("pole").count())
        .sum()
}

/// Raw period arithmetic in consumers (`periodic::modulo`, `offset_mod`,
/// `wrap_from`): the period-pattern ratchet routes every period query
/// through the service, but a consumer reducing parameters by hand still
/// owns a cut. Each site moves to `lift` / `aligned` / `canonical`, or is a
/// curve's own angle (a circle parameter), which moves to the curve's
/// parameter wrap.
fn arithmetic_count(text: &str) -> usize {
    code_lines(text)
        .map(|l| {
            [
                "periodic::modulo(",
                "periodic::offset_mod(",
                "periodic::wrap_from(",
            ]
            .iter()
            .map(|p| l.matches(p).count())
            .sum::<usize>()
        })
        .sum()
}

fn code_lines(text: &str) -> impl Iterator<Item = &str> {
    text.lines()
        .map(str::trim_start)
        .filter(|l| !l.starts_with("//"))
        .map(|l| l.split("//").next().unwrap_or(""))
}

fn counts_by(counter: fn(&str) -> usize) -> BTreeMap<String, usize> {
    let root = workspace_root();
    let mut out = BTreeMap::new();
    for file in popkern_lints::find_rs_files(&root.join("crates")) {
        let rel = file
            .strip_prefix(&root)
            .expect("under root")
            .to_string_lossy()
            .replace('\\', "/");
        let test_file =
            rel.contains("/tests/") || rel.ends_with("_tests.rs") || rel.ends_with("/tests.rs");
        if test_file || SERVICE.iter().any(|s| rel.starts_with(s)) {
            continue;
        }
        let n = counter(&std::fs::read_to_string(&file).unwrap_or_default());
        if n > 0 {
            out.insert(rel, n);
        }
    }
    out
}

fn read_baseline(name: &str) -> BTreeMap<String, usize> {
    std::fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join(name))
        .unwrap_or_default()
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| {
            let (file, n) = l.rsplit_once(' ').expect("`path count`");
            (file.to_string(), n.parse().expect("count"))
        })
        .collect()
}

fn check_ratchet(name: &str, what: &str, counter: fn(&str) -> usize) {
    let (base, now) = (read_baseline(name), counts_by(counter));
    let mut problems = Vec::new();
    for (file, n) in &now {
        match base.get(file) {
            Some(b) if n <= b => {}
            Some(b) => problems.push(format!(
                "{file}: {n} {what}, baseline {b}: use the chart service"
            )),
            None => problems.push(format!(
                "{file}: {n} {what}, none allowed: use the chart service"
            )),
        }
    }
    for (file, b) in &base {
        let n = now.get(file).copied().unwrap_or(0);
        if n < *b {
            problems.push(format!(
                "{file}: down to {n} from {b}: lower its entry in {name} (delete it at 0)"
            ));
        }
    }
    assert!(
        problems.is_empty(),
        "{} {name} problems:\n{}",
        problems.len(),
        problems.join("\n")
    );
}

fn write_ratchet(name: &str, header: &str, counter: fn(&str) -> usize) {
    let counts = counts_by(counter);
    let total: usize = counts.values().sum();
    let mut text = String::from(header);
    for (file, n) in counts {
        text.push_str(&format!("{file} {n}\n"));
    }
    text.push_str(&format!("# total {total}\n"));
    std::fs::write(Path::new(env!("CARGO_MANIFEST_DIR")).join(name), text).expect("write");
}

const POLE_HEADER: &str = "# Chart-service POLE ratchet baseline (PK-1228): per-file counts of the\n\
     # `pole` token in non-comment code outside the service. Only ever\n\
     # lowered; W2 exits with this file empty.\n";
const ARITHMETIC_HEADER: &str = "# Chart-service period-ARITHMETIC ratchet baseline (PK-1228): per-file\n\
     # counts of raw `periodic::` reductions outside the service. Only ever\n\
     # lowered; W2 exits with this file empty.\n";

#[test]
fn pole_logic_lives_only_in_the_chart_service() {
    check_ratchet("chart_pole_baseline.txt", "pole tokens", pole_count);
}

#[test]
fn period_arithmetic_lives_only_in_the_chart_service() {
    check_ratchet(
        "chart_arithmetic_baseline.txt",
        "raw period reductions",
        arithmetic_count,
    );
}

/// `CHART_RATCHET_WRITE=1` writes the pole and arithmetic baselines too.
#[test]
fn write_pole_and_arithmetic_baselines_when_asked() {
    if std::env::var("CHART_RATCHET_WRITE").is_err() {
        return;
    }
    write_ratchet("chart_pole_baseline.txt", POLE_HEADER, pole_count);
    write_ratchet(
        "chart_arithmetic_baseline.txt",
        ARITHMETIC_HEADER,
        arithmetic_count,
    );
}
