//! PK-304/PK-305 serving-progress pins: the first torus booleans the
//! kernel ever ran.
//!
//! History (2026-08-21): torus × cylinder failed after 54 s ("trim
//! curve end point does not match any vertex"), rotated torus × torus
//! failed after 1050 s. Three root fixes flipped both to completing:
//! the SSI memo (exactly-equal solver inputs solved once — a 16-patch
//! torus re-marched the identical intersection per pair), marched-loop
//! closure welding + seam segmentation (records stay continuous on the
//! principal branch), and the torn-boundary guard in the seam dance
//! (a ±π-shifted boundary whose u-span GREW is a fabricated bridge,
//! never a healed seam edge).
//!
//! These pins ratchet COMPLETION and wall clock only. Strict
//! `validate()` cleanliness of the outputs is the OPEN remainder of
//! PK-304/305 (duplicated cross-v-branch window claims leave open
//! intersection loops; design fork:
//! `docs/architecture/periodic-chart-cell-clipping.md`) — do NOT add
//! an is-invalid assertion here (no either-way pins), and DO extend
//! this file with validity gates when the chart-cell clip lands.
//!
//! Release-only: the ops are seconds-scale; debug wall clock would
//! dominate the sweep.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::time::{Duration, Instant};

fn debug_skip(name: &str) -> bool {
    if cfg!(debug_assertions) {
        eprintln!("torus boolean pin [{name}]: skipped in debug — runs in the CI release perf job");
        return true;
    }
    false
}

#[test]
fn torus_cylinder_union_completes_within_budget() {
    if debug_skip("torus x cylinder") {
        return;
    }
    let tol = Tolerance::default();
    let torus = BRep::torus(1.0, 0.35, tol).expect("torus builds");
    let cyl = BRep::cylinder(0.5, 3.0, tol)
        .expect("cylinder builds")
        .transform(&Transform3::translation(0.8, 0.0, -1.5))
        .expect("translate");
    let t0 = Instant::now();
    let result = boolean(BoolOperation::Union, torus, cyl, tol, None);
    let wall = t0.elapsed();
    let brep = result.expect("torus x cylinder union completes (PK-304)");
    // Completion + structural closure are the ratcheted claims; the
    // carve itself is the OPEN acceptance (current output retains
    // interior walls — validate_all's self-intersection screen
    // honestly flags them; chart-cell campaign continues).
    assert!(
        !brep.faces.is_empty(),
        "union must produce faces, got {}",
        brep.faces.len()
    );
    assert!(
        brep.validate().is_ok(),
        "union must validate structurally: {:?}",
        brep.validate().err()
    );
    // Measured 3.0-3.2 s release (2026-08-21, SSI memo) × 5 ≈ 16 s;
    // was 54 s before the memo — this bound also guards the memo.
    assert!(
        wall < Duration::from_secs(16),
        "torus x cylinder union took {wall:?} (budget 16 s = measured 3.2 s x 5)"
    );
}

/// PK-305 completion ratchet: the rotated torus × torus union
/// COMPLETES (seconds-scale; burned 17.5 minutes then failed before
/// PK-304, then declined transactionally, and completes since the
/// seam-honest torus pcurve authoring + chart-cell containment).
/// `validate()` cleanliness of the output is the OPEN acceptance
/// (current defect: one open intersection loop) — extend with
/// validity gates as the chart-cell campaign closes; never assert the
/// defect itself.
#[test]
fn rotated_torus_torus_union_completes_fast() {
    if debug_skip("torus x torus") {
        return;
    }
    let tol = Tolerance::default();
    let torus = BRep::torus(1.0, 0.35, tol).expect("torus builds");
    let torus_b = BRep::torus(1.0, 0.35, tol)
        .expect("torus builds")
        .transform(&Transform3::rotation_x(0.5).compose(&Transform3::translation(0.7, 0.4, 0.15)))
        .expect("transform");
    let t0 = Instant::now();
    let result = boolean(BoolOperation::Union, torus, torus_b, tol, None);
    let wall = t0.elapsed();
    let brep = result.expect("rotated torus x torus union completes (PK-305)");
    assert!(!brep.faces.is_empty());
    assert!(
        wall < Duration::from_secs(60),
        "must stay fast (was 1050 s before PK-304), took {wall:?}"
    );
}
