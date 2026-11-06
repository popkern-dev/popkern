//! The general path as the shortcuts' oracle (plan §4.3a).
//!
//! The touching-only and coincident-containment shortcuts answer a
//! boolean without the cellular stage. Resolution subsumes them: point
//! and line contacts are candidates, and containment is membership's
//! isolated-component case. They stay as performance paths only while
//! the general path agrees with them, so with
//! `POPKERN_SHORTCUT_ORACLE=1` every firing also runs the general path
//! and logs a `ShortcutOracle` shadow line: whether the two bodies
//! agree on membership at a deterministic probe lattice over both
//! operands' boxes, the probes that disagree, and both bodies'
//! topology counts.

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::PointClassification;
use crate::error::{BoolError, BoolOperation};
use crate::intersect::IntersectionResult;

/// Probes per axis of the lattice.
const LATTICE: usize = 9;

/// Whether shortcut firings are checked against the general path.
pub(super) fn oracle_enabled() -> bool {
    std::env::var_os("POPKERN_SHORTCUT_ORACLE").is_some_and(|v| v == "1")
}

/// Run the general path on the boolean `shortcut` answered and log
/// whether the two agree.
pub(super) fn compare(
    name: &str,
    op: BoolOperation,
    (a, b): (&BRep, &BRep),
    ix: &IntersectionResult,
    tolerance: &Tolerance,
    policy: crate::options::NonManifoldPolicy,
    shortcut: &Result<BRep, BoolError>,
) {
    let mut metrics = Vec::new();
    let general =
        super::serve::select_and_serve(a, b, ix, op, tolerance, &mut metrics, &mut None, policy);
    let verdict = match (shortcut, &general) {
        (Ok(s), Ok(g)) => {
            let probes = lattice(ix);
            let (sb, gb) = (
                crate::intersect::compute_face_aabbs(s),
                crate::intersect::compute_face_aabbs(g),
            );
            let mismatched = probes
                .iter()
                .filter(|p| membership(s, &sb, p, tolerance) != membership(g, &gb, p, tolerance))
                .count();
            format!(
                "{} probes={} mismatched={mismatched} short={} general={}",
                if mismatched == 0 { "agree" } else { "disagree" },
                probes.len(),
                counts(s),
                counts(g)
            )
        }
        (Err(_), Err(_)) => "agree both-declined".to_string(),
        (Ok(s), Err(e)) => format!("disagree general-declined short={} ({e})", counts(s)),
        (Err(e), Ok(g)) => format!("disagree short-declined general={} ({e})", counts(g)),
    };
    crate::complex_shadow::append_log_line(&format!("ShortcutOracle {name} {op:?} {verdict}"));
    if verdict.starts_with("disagree") {
        emit_operands(name, op, (a, b));
    }
}

/// With `POPKERN_COMPLEX_ARTIFACT=<dir>`, write a disagreeing boolean's
/// operands to `<dir>` as JSON (`{op, a, b}`), so it reproduces without
/// the model that made it.
fn emit_operands(name: &str, op: BoolOperation, (a, b): (&BRep, &BRep)) {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static SEQUENCE: AtomicUsize = AtomicUsize::new(0);
    let Some(dir) = std::env::var_os("POPKERN_COMPLEX_ARTIFACT") else {
        return;
    };
    let n = SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let path = std::path::Path::new(&dir).join(format!(
        "shortcut-oracle-{name}-{}-{n}.json",
        std::process::id()
    ));
    let body = serde_json::json!({ "op": format!("{op:?}"), "a": a, "b": b });
    if std::fs::write(&path, body.to_string()).is_ok() {
        eprintln!("[shortcut oracle] operands: {}", path.display());
    }
}

/// V/E/F/shells.
fn counts(b: &BRep) -> String {
    format!(
        "V{}E{}F{}S{}",
        b.vertices.len(),
        b.edges.len(),
        b.faces.len(),
        b.back_shell_count()
    )
}

/// Inside (true), outside (false), or no verdict (`None`): on the
/// boundary or unclassifiable. An empty body holds nothing. A probe
/// one body has no verdict on counts as a disagreement only when the
/// other has one.
fn membership(
    body: &BRep,
    boxes: &[(FaceHandle, Aabb)],
    p: &Point3,
    tolerance: &Tolerance,
) -> Option<bool> {
    if body.faces.is_empty() {
        return Some(false);
    }
    match crate::classify::point_in_solid_final(p, body, boxes, tolerance) {
        Ok(PointClassification::Inside) => Some(true),
        Ok(PointClassification::Outside) => Some(false),
        _ => None,
    }
}

/// A lattice over both operands' face boxes, its probes set off the
/// cell centres by irrational fractions so none lands on an axis plane
/// a box face or a flush contact lies on.
fn lattice(ix: &IntersectionResult) -> Vec<Point3> {
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for (_, bx) in ix.aabbs_a.iter().chain(&ix.aabbs_b) {
        for (i, (l, h)) in [
            (bx.min.x, bx.max.x),
            (bx.min.y, bx.max.y),
            (bx.min.z, bx.max.z),
        ]
        .into_iter()
        .enumerate()
        {
            lo[i] = lo[i].min(l);
            hi[i] = hi[i].max(h);
        }
    }
    if (0..3).any(|i| lo[i].partial_cmp(&hi[i]) != Some(std::cmp::Ordering::Less)) {
        return Vec::new();
    }
    let at = |i: usize, k: usize, shift: f64| {
        lo[i] + (hi[i] - lo[i]) * (k as f64 + shift) / LATTICE as f64
    };
    let mut out = Vec::with_capacity(LATTICE.pow(3));
    for x in 0..LATTICE {
        for y in 0..LATTICE {
            for z in 0..LATTICE {
                out.push(Point3::new(
                    at(0, x, 0.381_966),
                    at(1, y, 0.414_214),
                    at(2, z, 0.732_051),
                ));
            }
        }
    }
    out
}

#[cfg(test)]
#[path = "operation_shortcut_oracle_tests.rs"]
mod tests;
