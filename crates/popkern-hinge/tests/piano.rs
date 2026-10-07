//! PK-697: a piano (continuous) hinge is a butt hinge with many
//! knuckles. One comb cut per leaf whatever the count keeps the number
//! of booleans constant; this pins watertightness, the build time and
//! its growth.

use std::time::Instant;

use popkern_core::tolerance::Tolerance;
use popkern_hinge::{HingeComponent, HingeKind, HingeParams, hinge_component};

fn piano(knuckles: u32) -> HingeParams {
    HingeParams {
        knuckles,
        leaf_width: 20.0,
        thickness: 1.5,
        pin_diameter: 2.5,
        clearance: 0.15,
        gap: 0.3,
        ..HingeParams::new(300.0)
    }
}

/// Seconds to build both leaves and the pin (meshing, timed apart, is
/// the display's cost, not the kernel's).
fn build_seconds(knuckles: u32) -> f64 {
    let p = piano(knuckles);
    let mut build = 0.0;
    let mut mesh_s = 0.0;
    for c in [
        HingeComponent::LeafA,
        HingeComponent::LeafB,
        HingeComponent::Pin,
    ] {
        let start = Instant::now();
        let body = hinge_component(HingeKind::Butt, c, &p, Tolerance::default())
            .unwrap_or_else(|e| panic!("{knuckles} knuckles {c:?}: {e}"));
        build += start.elapsed().as_secs_f64();
        let start = Instant::now();
        let mesh = popkern_mesh::tessellate(&body, &popkern_mesh::TessellationOptions::display())
            .expect("mesh");
        mesh_s += start.elapsed().as_secs_f64();
        assert_eq!(
            popkern_mesh::welded_open_edges(&mesh),
            0,
            "{knuckles} {c:?}"
        );
    }
    eprintln!("piano 300 x {knuckles}: build {build:.2} s, display mesh {mesh_s:.2} s");
    build
}

#[test]
fn piano_hinges_build_watertight_within_budget_and_linearly() {
    let counts = [10_u32, 30, 60, 100];
    let times: Vec<f64> = counts.iter().map(|&n| build_seconds(n)).collect();
    // Design budget (§8): 60 knuckles under 2 s natively (test profile).
    // Measured 2026-10-07 after PK-708 (the shared-edge split's
    // bounds prefilter): 0.95 / 1.61 / 2.63 s at 30 / 60 / 100; before,
    // 1.56 / 3.94 / 9.24 s with the merge phase growing ~F^1.9.
    // The absolute budget only means something on a quiet machine: a
    // full parallel suite stretches wall time (measured 2.79 s here).
    if std::env::var_os("POPKERN_TIMING_GATES").is_some() {
        assert!(times[2] < 2.0, "60 knuckles took {:.2} s", times[2]);
    }
    // Near-linear (load stretches both ends alike): 10x the knuckles
    // well under the ~100x a quadratic step would cost.
    assert!(
        times[3] < 15.0 * times[0].max(0.05),
        "100 vs 10 knuckles: {:.2} s vs {:.2} s",
        times[3],
        times[0]
    );
}

/// PK-707: the leaves' coaxial bore through a 300-long knuckle prism
/// serves on the cellular pipeline (ray classification once filtered
/// every cylinder hit beyond |v| = 100, so the cap's bore disc read as
/// kept and the boolean fell back to graph merge).
#[test]
fn a_long_piano_leaf_bores_without_fallback() {
    let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
    for c in [HingeComponent::LeafA, HingeComponent::LeafB] {
        hinge_component(HingeKind::Butt, c, &piano(10), Tolerance::default()).expect("leaf");
    }
    let events = popkern_bool::merge_fallback::take_merge_fallback_events();
    assert!(events.is_empty(), "{events:?}");
}
