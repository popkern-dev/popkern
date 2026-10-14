//! Boolean CI performance budget gate (PK-272).
//!
//! Absolute wall-clock budgets for end-to-end boolean operations, in
//! the PK-190 pattern (`popkern-blend/tests/perf_budget.rs`): each
//! budget is `measured × 5 headroom`, measured as best-of-5 release
//! wall clock on the dev reference machine (16-core linux). A
//! regression past a budget fails the build. Budgets are meaningless
//! in debug builds — every test self-skips outside release. CI runs
//! this file in the release-mode `perf` job; locally use
//! `tools/quality-check.sh --perf`.
//!
//! # Why these cases
//!
//! Before PK-272, boolean end-to-end performance was entirely
//! unmeasured, and every boolean ran an unconditional `validate_all()`
//! (O(n²) surface-surface self-intersection sweep) purely to print a
//! stderr diagnostic — ungated in release. The real correctness gate
//! is `validate::post_check` (which runs either way); the diagnostic
//! is now debug/the `boolean` trace topic-only. The cases below cover the cost
//! regimes that change decision:
//!
//! - **box-box** — minimal all-planar op: pipeline overhead floor.
//! - **~100-face prism** — face-count scaling on all-planar input.
//! - **filleted-operand** — operand carrying fillet strips and corner
//!   patches (cylinder/sphere/NURBS), the normal CAD workflow and the
//!   regime where the ungated `validate_all` cost the most (~1/3 of
//!   the op at bench scale).
//! - **torus-carrying** — operand with a periodic torus fillet strip,
//!   the slowest surface class in the self-intersection sweep (no
//!   analytic torus SSI arm — PK-279).
//!
//! # Before/after at the same sites (PK-272 acceptance row 403)
//!
//! Best-of-5 release wall clock, single-threaded, 2026-08-14 dev
//! reference machine, before → after gating the `validate_all` stderr
//! diagnostic and the compute-and-discard sliver/watertightness
//! diagnostics:
//!
//! | case | before | after |
//! |---|---|---|
//! | box-box union | 9.55 ms | 9.49 ms |
//! | ~100-face prism union | 29.54 ms | 27.91 ms |
//! | filleted-cube pocket subtraction | 26.10 ms | 17.65 ms |
//! | torus-carrying cylinder union | 25.38 ms | 24.24 ms |
//!
//! The scope's "possibly 10–100× on filleted operands" hypothesis was
//! REFUTED at these operand sizes: the largest win is 32% (filleted
//! cube, ~30 curved faces). The sweep is O(n²) in face count, so the
//! gate matters more as parts grow, but at bench scale the boolean
//! pipeline itself dominates.

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::time::{Duration, Instant};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Skip guard: budgets only bind in release builds.
fn debug_skip(name: &str) -> bool {
    if cfg!(debug_assertions) {
        eprintln!(
            "PK-272 perf budget [{name}]: skipped in debug — enforced by the CI release perf job"
        );
        return true;
    }
    false
}

/// Best-of-N wall clock for `op`, with one warmup call.
fn best_of<F: FnMut()>(n: usize, mut op: F) -> Duration {
    op();
    let mut best = Duration::MAX;
    for _ in 0..n {
        let t0 = Instant::now();
        op();
        best = best.min(t0.elapsed());
    }
    best
}

/// Assert `measured` is within `budget_ms`, printing the number either
/// way so CI logs double as a coarse performance history.
fn assert_budget(name: &str, measured: Duration, budget_ms: u64) {
    let budget = Duration::from_millis(budget_ms);
    eprintln!("PK-272 perf budget [{name}]: best-of-run {measured:?} (budget {budget:?})");
    assert!(
        measured <= budget,
        "PK-272 performance regression: {name} took {measured:?}, budget {budget:?} \
         (= measured baseline × 5 headroom; see file header). A jump past this budget \
         means a fast path stopped applying or an ungated diagnostic sweep returned — \
         check log_validation_issues / log_post_processing_diagnostics gating first."
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

/// Minimal all-planar boolean: two overlapping unit cubes. The
/// pipeline-overhead floor (intersection, split, classify, merge,
/// post-check on ~12 faces).
#[test]
fn box_box_union_within_budget() {
    if debug_skip("box-box union") {
        return;
    }
    let a = BRep::unit_cube(tol()).expect("unit cube builds");
    let b = axis_box(Point3::new(0.5, 0.5, 0.5), Point3::new(1.5, 1.5, 1.5));
    let measured = best_of(5, || {
        let r = popkern_bool::union_ref(&a, &b, tol(), None).expect("box-box union succeeds");
        assert!(!r.faces.is_empty(), "union result must be non-empty");
    });
    // Measured 9.49 ms best-of-5 release (2026-08-14) × 5 ≈ 48 ms.
    assert_budget("box-box union", measured, 48);
}

/// Face-count scaling on all-planar input: a 48-gon prism (50 faces)
/// united with a box piercing its side wall — dozens of narrow
/// plane-plane intersections at varied incidence angles.
#[test]
fn many_face_prism_union_within_budget() {
    if debug_skip("~100-face prism union") {
        return;
    }
    let n = 48usize;
    let profile: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
            Point3::new(theta.cos(), theta.sin(), 0.0)
        })
        .collect();
    let prism =
        BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), 1.0, tol()).expect("48-gon prism builds");
    let boss = axis_box(Point3::new(0.5, -0.25, 0.25), Point3::new(1.5, 0.25, 0.75));
    let measured = best_of(5, || {
        let r =
            popkern_bool::union_ref(&prism, &boss, tol(), None).expect("prism-boss union succeeds");
        assert!(
            r.faces.len() > 50,
            "union must keep the prism's facet count (got {})",
            r.faces.len()
        );
    });
    // Measured 27.91 ms best-of-5 release (2026-08-14) × 5 ≈ 140 ms.
    assert_budget("~100-face prism union", measured, 140);
}

/// The normal CAD workflow: boolean on an operand carrying fillet
/// strips. All 12 cube edges filleted (cylinder strips + corner
/// patches), then a pocket subtracted from the top face, clear of the
/// fillets — the cut itself is planar, so this case isolates what the
/// result's curved faces cost in post-op validation.
#[test]
fn filleted_cube_pocket_subtraction_within_budget() {
    if debug_skip("filleted-cube pocket subtraction") {
        return;
    }
    let cube = BRep::unit_cube(tol()).expect("unit cube builds");
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.05),
        })
        .collect();
    let filleted = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None)
        .expect("12-edge cube fillet succeeds")
        .brep;
    let pocket = axis_box(Point3::new(0.4, 0.4, 0.6), Point3::new(0.6, 0.6, 1.4));
    let measured = best_of(5, || {
        let r = popkern_bool::subtraction_ref(&filleted, &pocket, tol(), None)
            .expect("pocket subtraction succeeds");
        assert!(
            r.faces.len() > filleted.faces.len(),
            "subtraction must add pocket walls (got {} faces)",
            r.faces.len()
        );
    });
    // Measured 17.65 ms best-of-5 release (2026-08-14) × 5 ≈ 90 ms.
    // Before PK-272 gating this case measured 26.10 ms — the ungated
    // validate_all diagnostic was ~1/3 of the whole op.
    assert_budget("filleted-cube pocket subtraction", measured, 90);
}

/// Operand carrying a periodic torus fillet strip (cylinder base
/// fillet), united with a box piercing the cylinder wall. Torus faces
/// take the slowest surface-surface path (no analytic torus arm,
/// PK-279); at this face count the sweep cost ~1 ms, but the case
/// pins the slowest surface class so a torus-path regression is
/// caught here first.
#[test]
fn torus_carrying_cylinder_union_within_budget() {
    if debug_skip("torus-carrying cylinder union") {
        return;
    }
    let cyl = BRep::cylinder(1.0, 2.0, tol()).expect("cylinder builds");
    let base_edge = cyl
        .edges
        .keys()
        .find(|&e| {
            let he = cyl.edges[e].half_edge;
            let twin = cyl.half_edges[he].twin;
            let f1 = cyl.loops[cyl.half_edges[he].loop_].face;
            let f2 = cyl.loops[cyl.half_edges[twin].loop_].face;
            if f1 == f2 {
                return false;
            }
            // Base circle: both endpoints at z = 0.
            let a = cyl.vertices[cyl.half_edges[he].origin].position;
            let b = cyl.vertices[cyl.half_edges[twin].origin].position;
            a.z.abs() < 1e-9 && b.z.abs() < 1e-9
        })
        .expect("cylinder has a base edge");
    let specs = vec![FilletSpec {
        edge: base_edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let filleted = fillet_edges(&cyl, &specs, &FilletOptions::default(), tol(), None)
        .expect("cylinder base fillet succeeds")
        .brep;
    let boss = axis_box(Point3::new(0.5, -0.25, 0.75), Point3::new(1.5, 0.25, 1.25));
    let measured = best_of(5, || {
        let r = popkern_bool::union_ref(&filleted, &boss, tol(), None)
            .expect("torus-carrying union succeeds");
        assert!(
            r.faces.len() >= filleted.faces.len(),
            "union must keep the filleted body's faces (got {})",
            r.faces.len()
        );
    });
    // Measured 24.24 ms best-of-5 release (2026-08-14) × 5 ≈ 125 ms.
    assert_budget("torus-carrying cylinder union", measured, 125);
}

/// The 64-slot disc (the six-axis arm's R4 fixture, PK-743): the 64th
/// slot cut against the first. Point-in-face tests visited every hole
/// in full: the clip's hole containment sampled each hole contour, each
/// classification ray hit sampled each inner loop, and every fragment
/// sample walked every edge of the body. Conservative boxes (pcurve
/// forms, edge curves) skip them with the same verdicts. Measured
/// 3.76-3.81 interleaved (slot 0 ~20 ms, slot 63 74-78 ms); with either
/// contour prescreen off 6.9-7.6 on a quiet machine. Heavy contention
/// compresses the ratio (a control read 4.06 at load 30), so the gate
/// sat at 4.6: it can pass under load, never fail because of it.
/// PK-802 then took the walk off untouched faces (identity fragments),
/// off a cap's far slot holes (peeled holes), the quadratic node
/// lookup out of graph build, and rays out of points beyond the
/// closed body's box: 2.97-3.05 (slot 0 ~9.9 ms, slot 63 ~30 ms), gate
/// 3.6. The rest is the serve postprocess and post-check over the
/// whole body (PK-804) and naming lineage (PK-803).
#[test]
fn slotted_disc_late_cut_stays_near_the_first() {
    if debug_skip("64-slot disc, slot 63 vs slot 0") {
        return;
    }
    let tol = tol();
    let slot = |k: usize| {
        let a = std::f64::consts::TAU * k as f64 / 64.0;
        let (ca, sa, w) = (a.cos(), a.sin(), 0.6);
        let pts = [
            Point3::new(22.0 * ca + w * sa, 22.0 * sa - w * ca, -1.0),
            Point3::new(27.0 * ca + w * sa, 27.0 * sa - w * ca, -1.0),
            Point3::new(27.0 * ca - w * sa, 27.0 * sa + w * ca, -1.0),
            Point3::new(22.0 * ca - w * sa, 22.0 * sa + w * ca, -1.0),
        ];
        BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), 4.0, tol).expect("slot prism")
    };
    let cut = |body: &BRep, k: usize| {
        boolean(BoolOperation::Subtraction, body.clone(), slot(k), tol, None).expect("slot cut")
    };
    let disc = BRep::cylinder(30.0, 1.5, tol).expect("disc");
    let mut body = disc.clone();
    for k in 0..63 {
        body = cut(&body, k);
    }
    // Interleaved, so both sides of the ratio see the same contention.
    let (mut first, mut last) = (Duration::MAX, Duration::MAX);
    for _ in 0..7 {
        let t0 = Instant::now();
        cut(&disc, 0);
        first = first.min(t0.elapsed());
        let t0 = Instant::now();
        cut(&body, 63);
        last = last.min(t0.elapsed());
    }
    let ratio = last.as_secs_f64() / first.as_secs_f64().max(1e-9);
    eprintln!(
        "perf budget [slotted disc]: slot 0 {first:?}, slot 63 {last:?}, ratio {ratio:.2} (gate 3.6)"
    );
    assert!(
        ratio <= 3.6,
        "the 64th slot cut costs {ratio:.2}x the first ({last:?} vs {first:?}) > 3.6"
    );
}
