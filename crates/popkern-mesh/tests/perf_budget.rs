//! Tessellation CI performance budget gate (PK-277).
//!
//! Absolute wall-clock budgets for end-to-end tessellation, in the
//! PK-190/PK-272 pattern (`popkern-bool/tests/perf_budget.rs`): each
//! budget is `measured × 5 headroom`, measured as best-of-5 release
//! wall clock on the 16-core reference machine with the default rayon
//! pool. A regression past a budget fails the build. Budgets are
//! meaningless in debug builds, so the tests skip there — CI enforces
//! them in the release perf job.
//!
//! These fixtures were introduced with the PK-277 parallel per-face
//! dispatch and the incremental encroachment rounds, so their budgets
//! also guard those two mechanisms against silent regression to
//! serial or all-faces behavior.

use std::time::{Duration, Instant};

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::TessellationOptions;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Skip guard: budgets only bind in release builds.
fn debug_skip(name: &str) -> bool {
    if cfg!(debug_assertions) {
        eprintln!(
            "PK-277 perf budget [{name}]: skipped in debug — enforced by the CI release perf job"
        );
        return true;
    }
    false
}

fn best_of<F: FnMut()>(n: usize, mut op: F) -> Duration {
    let mut best = Duration::MAX;
    for _ in 0..n {
        let t0 = Instant::now();
        op();
        best = best.min(t0.elapsed());
    }
    best
}

/// Assert `measured` is within `budget_ms`, printing the number either
/// way so perf history is greppable from CI logs.
fn assert_budget(name: &str, measured: Duration, budget_ms: u64) {
    let budget = Duration::from_millis(budget_ms);
    eprintln!("PK-277 perf budget [{name}]: best-of-run {measured:?} (budget {budget:?})");
    assert!(
        measured <= budget,
        "PK-277 performance regression: {name} took {measured:?}, budget {budget:?} \
         (= measured baseline × 5 headroom; see file header)."
    );
}

/// Axis-aligned box from `min` to `max` via the prism constructor.
fn axis_box(min: Point3, max: Point3) -> BRep {
    let profile = [
        Point3::new(min.x, min.y, min.z),
        Point3::new(max.x, min.y, min.z),
        Point3::new(max.x, max.y, min.z),
        Point3::new(min.x, max.y, min.z),
    ];
    BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), max.z - min.z, tol()).expect("axis box builds")
}

fn regular_prism(n: usize, radius: f64, height: f64) -> BRep {
    let profile: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
            Point3::new(radius * theta.cos(), radius * theta.sin(), 0.0)
        })
        .collect();
    BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), height, tol()).expect("prism builds")
}

/// Face-count scaling: a 48-gon prism united with a side boss (~55
/// faces), tessellated at a fine chord height — the per-face parallel
/// dispatch is the dominant axis.
#[test]
fn many_face_union_tessellation_within_budget() {
    if debug_skip("~55-face union tessellation") {
        return;
    }
    let prism = regular_prism(48, 1.0, 1.0);
    let boss = axis_box(Point3::new(0.5, -0.25, 0.25), Point3::new(1.5, 0.25, 0.75));
    let solid = popkern_bool::union_ref(&prism, &boss, tol(), None).expect("union succeeds");
    let options = TessellationOptions {
        chord_height: 0.001,
        ..Default::default()
    };
    let measured = best_of(5, || {
        let mesh = popkern_mesh::tessellate(&solid, &options).expect("tessellation succeeds");
        assert!(
            mesh.triangles.len() > 100,
            "mesh must be non-trivial (got {} triangles)",
            mesh.triangles.len()
        );
    });
    // Measured 5.2 ms best-of-5 release, 16-core pool (2026-08-20) × 5 ≈ 26 ms.
    assert_budget("~55-face union tessellation", measured, 26);
}

/// The PK-295 served-subtraction fixture: a rim-filleted cylinder
/// minus a strip-crossing cuboid. Measured to drive the encroachment
/// refinement (2 rounds, 75 + 65 span splits with default options),
/// so this budget guards the incremental re-check — only faces
/// incident to split edges re-run per round.
#[test]
fn encroaching_subtraction_tessellation_within_budget() {
    if debug_skip("encroaching subtraction tessellation") {
        return;
    }
    use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
    let cyl = BRep::cylinder(1.0, 2.0, tol()).expect("cylinder builds");
    let rim: Vec<_> = cyl
        .edges
        .iter()
        .filter(|(_, d)| {
            let he = &cyl.half_edges[d.half_edge];
            let a = cyl.vertices[he.origin].position;
            let b = cyl.vertices[cyl.half_edges[he.twin].origin].position;
            let f1 = cyl.loops[he.loop_].face;
            let f2 = cyl.loops[cyl.half_edges[he.twin].loop_].face;
            f1 != f2 && (a.z - 2.0).abs() < 1e-9 && (b.z - 2.0).abs() < 1e-9
        })
        .map(|(e, _)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.2),
        })
        .collect();
    let part = fillet_edges(&cyl, &rim, &FilletOptions::default(), tol(), None)
        .expect("rim fillet completes")
        .brep;
    let cutter = BRep::unit_cube(tol())
        .expect("cube builds")
        .transform(&popkern_core::geometry::Transform3::scale_uniform(1.8))
        .expect("scales")
        .transform(&popkern_core::geometry::Transform3::translation(
            -0.9, -0.9, 0.6,
        ))
        .expect("translates");
    let solid = popkern_bool::subtraction(part, cutter, tol(), None).expect("subtraction succeeds");
    let options = TessellationOptions::default();
    let measured = best_of(5, || {
        let mesh = popkern_mesh::tessellate(&solid, &options).expect("tessellation succeeds");
        assert!(
            mesh.triangles.len() > 100,
            "mesh must be non-trivial (got {} triangles)",
            mesh.triangles.len()
        );
    });
    // Measured 25.2 ms best-of-5 release, 16-core pool (2026-08-20) × 5 ≈ 126 ms.
    assert_budget("encroaching subtraction tessellation", measured, 126);
}
