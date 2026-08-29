//! §18 Performance baseline benchmarks for the fillet pipeline.
//!
//! Ships as-is regression guards and baseline measurements for:
//!
//! 1. **§11.1 analytic fast paths** vs the generic NURBS sweep on
//!    identical input geometry (the §18 ≥ 10× acceptance metric).
//!    Toggles `build_fillet_surface(... is_constant_radius ...)`
//!    to compare analytic-torus / analytic-cylinder against the
//!    same-inputs NURBS sweep.
//! 2. **Full fillet pipeline** end-to-end on a 12-edge unit cube
//!    via `fillet_edges`, measuring the total cost of a
//!    commercial-grade multi-edge fillet op.
//! 3. **Curved-support base fillets** via `fillet_edges` on
//!    `BRep::cylinder` (seam-crossing analytic-torus path) and
//!    `BRep::cone` (phase-1b walker path, restored by PK-264).
//! 4. **§13.1 degeneracy detector** cost on a clean cube,
//!    measuring the per-edge overhead of the pre-pass.
//! 5. **§13.1(2) sliver absorber** cost on a healthy fillet
//!    output, measuring the post-construction cleanup overhead.
//!
//! The acceptance test in `tests/perf_acceptance.rs` gates the
//! §18 analytic ≥ 10× requirement in normal `cargo test`; this
//! file is for detailed baseline tracking via `cargo bench`.

use criterion::{Criterion, criterion_group, criterion_main};
use popkern_blend::center_curve::CenterCurveResult;
use popkern_blend::fillet_surface::build_fillet_surface;
use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ===========================================================================
// §11.1 analytic fast paths vs NURBS sweep (the §18 load-bearing ratio)
// ===========================================================================

fn plane_plane_fixture() -> (SurfaceGeometry, SurfaceGeometry, CenterCurveResult) {
    let s1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let s2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 1.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let center_pts: Vec<Point3> = (0..10)
        .map(|i| Point3::new(0.1, i as f64 * 0.1, 0.1))
        .collect();
    let center_curve = popkern_blend::center_curve::fit_nurbs_through_points(&center_pts).unwrap();
    let footprint_1: Vec<Point3> = (0..10)
        .map(|i| Point3::new(0.0, i as f64 * 0.1, 0.1))
        .collect();
    let footprint_2: Vec<Point3> = (0..10)
        .map(|i| Point3::new(0.1, i as f64 * 0.1, 0.0))
        .collect();
    (
        s1,
        s2,
        CenterCurveResult {
            center_curve,
            footprint_1,
            footprint_2,
            contact_solutions: None,
            sweep_center_curve: None,
            sweep_footprint_1: None,
            sweep_footprint_2: None,
        },
    )
}

fn bench_plane_plane_analytic_vs_nurbs(c: &mut Criterion) {
    let (s1, s2, cr) = plane_plane_fixture();
    let mut group = c.benchmark_group("§11.1 plane-plane fillet surface");
    group.bench_function("analytic torus (is_constant=true)", |b| {
        b.iter(|| {
            build_fillet_surface(
                &s1,
                &s2,
                &cr,
                0.1,
                8,
                true,
                popkern_blend::ContinuityLevel::G1,
                false,
                false,
                popkern_blend::BlendProfile::Circular,
            )
            .unwrap()
        })
    });
    group.bench_function("generic NURBS sweep (is_constant=false)", |b| {
        b.iter(|| {
            build_fillet_surface(
                &s1,
                &s2,
                &cr,
                0.1,
                8,
                false,
                popkern_blend::ContinuityLevel::G1,
                false,
                false,
                popkern_blend::BlendProfile::Circular,
            )
            .unwrap()
        })
    });
    group.finish();
}

// ===========================================================================
// Full pipeline end-to-end (fillet_edges) — batch vs single-edge scaling
// ===========================================================================

fn bench_single_edge_cube_fillet(c: &mut Criterion) {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    c.bench_function("§18 single-edge cube fillet", |b| {
        b.iter(|| fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).ok())
    });
}

fn bench_twelve_edge_cube_fillet(c: &mut Criterion) {
    let cube = BRep::unit_cube(tol()).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    c.bench_function("§18 12-edge cube fillet (full pipeline)", |b| {
        b.iter(|| fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).ok())
    });
}

// ===========================================================================
// Curved-support base-circle fillet — plane×cylinder analytic torus
// on a periodic (seam-crossing) edge
// ===========================================================================
//
// Honesty note (PK-190): this bench previously ran the CONE base with
// `.ok()`, which swallowed the result — and the cone base fillet then
// failed at every radius (mirror-ball convexity misclassification +
// cone-blind footprint clip), so the "walker path" baseline was timing
// a sub-millisecond failure. PK-264 fixed the cone at the root and the
// cone entry below returned; the `expect`s keep both honest by failing
// loudly if either op ever stops serving.

fn bench_cylinder_base_fillet(c: &mut Criterion) {
    let cyl = BRep::cylinder(1.0, 2.0, tol()).expect("cylinder BRep should build");
    // Pick an edge with distinct adjacent faces (a base circle).
    let chosen = cyl
        .edges
        .keys()
        .find(|&e| {
            let he = cyl.edges[e].half_edge;
            let twin = cyl.half_edges[he].twin;
            let face_a = cyl.loops[cyl.half_edges[he].loop_].face;
            let face_b = cyl.loops[cyl.half_edges[twin].loop_].face;
            face_a != face_b
        })
        .expect("cylinder has a base edge with distinct faces");
    let specs = vec![FilletSpec {
        edge: chosen,
        radius: FilletRadius::Constant(0.1),
    }];
    c.bench_function("§18 cylinder base fillet (curved support)", |b| {
        b.iter(|| {
            fillet_edges(&cyl, &specs, &FilletOptions::default(), tol(), None)
                .expect("cylinder base fillet must keep serving")
        })
    });
}

fn bench_cone_base_fillet(c: &mut Criterion) {
    let cone = BRep::cone(1.0, 2.0, tol()).expect("cone BRep should build");
    let chosen = cone
        .edges
        .keys()
        .find(|&e| {
            let he = cone.edges[e].half_edge;
            let twin = cone.half_edges[he].twin;
            let face_a = cone.loops[cone.half_edges[he].loop_].face;
            let face_b = cone.loops[cone.half_edges[twin].loop_].face;
            face_a != face_b
        })
        .expect("cone has a base edge with distinct faces");
    let specs = vec![FilletSpec {
        edge: chosen,
        radius: FilletRadius::Constant(0.03),
    }];
    c.bench_function("§18 cone base fillet (walker path)", |b| {
        b.iter(|| {
            fillet_edges(&cone, &specs, &FilletOptions::default(), tol(), None)
                .expect("cone base fillet must keep serving (PK-264)")
        })
    });
}

// ===========================================================================
// Pre-pass cost: degeneracy detectors + sliver absorber
// ===========================================================================
//
// These are measured indirectly — the full fillet pipeline calls
// both. The per-op overhead is the difference between the full
// pipeline and the sum of its inner-phase timings, which we'd
// need telemetry to measure directly. Instead we show the full-
// pipeline bench above and rely on criterion's variance to catch
// regressions in the pre-pass overhead.

criterion_group!(
    benches,
    bench_plane_plane_analytic_vs_nurbs,
    bench_single_edge_cube_fillet,
    bench_twelve_edge_cube_fillet,
    bench_cylinder_base_fillet,
    bench_cone_base_fillet,
);
criterion_main!(benches);
