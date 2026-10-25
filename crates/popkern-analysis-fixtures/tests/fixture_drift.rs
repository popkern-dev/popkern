//! Each frozen fixture is what the generator builds now: the same STEP
//! entities in the same order, with every number within a relative
//! 1e-9 (an absolute 1e-9 below one). Kernel changes move last digits
//! all the time (a sphere pcurve's v by 2e-11 between 2026-10-10 and
//! 2026-10-25); a failure here means a construction change.
//! Regenerate with `cargo run -p popkern-analysis-fixtures`, inspect,
//! and record the diff in the fixture log
//! (docs/architecture/analysis-integration.md).

use std::path::Path;

/// Relative agreement each STEP number must keep.
const NUMBER_TOLERANCE: f64 = 1e-9;

/// The line with each real number replaced by `#`, and the numbers.
/// STEP reals always carry a decimal point; integers (entity ids,
/// counts) stay in the skeleton.
fn split_numbers(line: &str) -> (String, Vec<f64>) {
    let mut skeleton = String::with_capacity(line.len());
    let mut numbers = Vec::new();
    let bytes = line.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let starts = bytes[i].is_ascii_digit()
            || (bytes[i] == b'-' && bytes.get(i + 1).is_some_and(u8::is_ascii_digit));
        let after_id = i > 0 && bytes[i - 1] == b'#';
        if starts && !after_id {
            let mut j = i + 1;
            while j < bytes.len()
                && (bytes[j].is_ascii_digit()
                    || bytes[j] == b'.'
                    || bytes[j] == b'E'
                    || ((bytes[j] == b'-' || bytes[j] == b'+') && bytes[j - 1] == b'E'))
            {
                j += 1;
            }
            let token = &line[i..j];
            if token.contains('.')
                && let Ok(value) = token.parse::<f64>()
            {
                skeleton.push('#');
                numbers.push(value);
                i = j;
                continue;
            }
            skeleton.push_str(token);
            i = j;
            continue;
        }
        skeleton.push(bytes[i] as char);
        i += 1;
    }
    (skeleton, numbers)
}

/// Why `generated` is not `frozen` up to [`NUMBER_TOLERANCE`], or `None`.
fn drift(generated: &str, frozen: &str) -> Option<String> {
    let (a, b): (Vec<&str>, Vec<&str>) = (generated.lines().collect(), frozen.lines().collect());
    if a.len() != b.len() {
        return Some(format!("{} lines vs {} frozen", a.len(), b.len()));
    }
    for (n, (x, y)) in a.iter().zip(&b).enumerate() {
        let ((sx, nx), (sy, ny)) = (split_numbers(x), split_numbers(y));
        if sx != sy || nx.len() != ny.len() {
            return Some(format!("line {}: `{x}` vs frozen `{y}`", n + 1));
        }
        for (p, q) in nx.iter().zip(&ny) {
            if (p - q).abs() > NUMBER_TOLERANCE * q.abs().max(1.0) {
                return Some(format!("line {}: {p} vs frozen {q}", n + 1));
            }
        }
    }
    None
}

#[test]
fn frozen_fixtures_match_the_generator() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/analysis");
    for fixture in &popkern_analysis_fixtures::FIXTURES {
        let generated = popkern_analysis_fixtures::step_bytes(fixture)
            .unwrap_or_else(|e| panic!("{}: {e}", fixture.name));
        let path = dir.join(format!("{}.step", fixture.name));
        let frozen = std::fs::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let found = drift(
            &String::from_utf8_lossy(&generated),
            &String::from_utf8_lossy(&frozen),
        );
        assert!(
            found.is_none(),
            "{} drifted from the frozen file: {}",
            fixture.name,
            found.unwrap_or_default()
        );
    }
}

#[test]
fn the_drift_check_passes_last_digits_and_fails_construction_changes() {
    let frozen = "#180=CARTESIAN_POINT('',(6.283185307179585E0,-3.398369094345748E-1));";
    let noise = "#180=CARTESIAN_POINT('',(6.283185307179585E0,-3.398369094541238E-1));";
    assert_eq!(drift(noise, frozen), None);
    let moved = "#180=CARTESIAN_POINT('',(6.283185307179585E0,-3.3984E-1));";
    assert!(drift(moved, frozen).is_some());
    let renumbered = "#181=CARTESIAN_POINT('',(6.283185307179585E0,-3.398369094345748E-1));";
    assert!(drift(renumbered, frozen).is_some());
    assert!(drift(&format!("{frozen}\n{frozen}"), frozen).is_some());
}

#[test]
fn generation_is_deterministic() {
    for fixture in &popkern_analysis_fixtures::FIXTURES {
        let a = popkern_analysis_fixtures::step_bytes(fixture).unwrap_or_else(|e| panic!("{e}"));
        let b = popkern_analysis_fixtures::step_bytes(fixture).unwrap_or_else(|e| panic!("{e}"));
        assert!(a == b, "{} is not deterministic", fixture.name);
    }
}
