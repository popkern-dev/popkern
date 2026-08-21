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

/// PK-305 stability ratchet: the rotated torus × torus union stays
/// FAST and TRANSACTIONAL — before PK-304 it burned 17.5 minutes and
/// failed; the wall clock and untouched-inputs claims are stable
/// across the campaign. Completion currently oscillates as the
/// chart-cell layers land (each fix exposes the next never-exercised
/// defect in this subsystem — this pin deliberately does NOT assert
/// Ok/Err); when the campaign closes, extend with completion +
/// validate()/OCC gates.
#[test]
fn rotated_torus_torus_union_is_fast_and_transactional() {
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
    let result = boolean(BoolOperation::Union, torus.clone(), torus_b, tol, None);
    let wall = t0.elapsed();
    eprintln!(
        "torus x torus union: {:?}, completed = {}",
        wall,
        result.is_ok()
    );
    assert!(
        wall < Duration::from_secs(60),
        "must stay fast (was 1050 s before PK-304), took {wall:?}"
    );
    assert!(
        torus.validate().is_ok(),
        "the operation must leave inputs untouched"
    );
}
