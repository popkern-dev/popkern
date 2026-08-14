//! Criterion benchmarks for the §14.7 in-place fillet radius
//! modification fast paths against the remove+recreate baseline.
//!
//! Locks in the speedups Round 5d-2/5d-3/5h/5m delivered, so any
//! future regression is caught by `cargo bench -p popkern-fillet`.
//!
//! Coverage:
//!
//! | Case                             | Fast path                          | Baseline           |
//! |----------------------------------|------------------------------------|--------------------|
//! | Open plane-plane Constant        | `modify_open_in_place`             | unfillet+fillet    |
//! | Open plane-plane Variable        | `modify_open_in_place_variable`    | unfillet+fillet    |
//! | Open plane-plane ByPoints        | `modify_open_in_place_bypoints`    | unfillet+fillet    |
//! | Closed plane-cylinder Constant   | `modify_closed_in_place`           | unfillet+fillet    |

use criterion::{Criterion, criterion_group, criterion_main};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_fillet::{
    FilletOptions, FilletRadius, FilletSpec, fillet_edge, fillet_edges, modify_inplace, unfillet,
    unfillet_closed,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Pre-fillet a unit cube edge so the benchmarks operate on the
/// post-fillet brep (which is what `modify_*` takes as input).
fn open_setup() -> (BRep, FaceHandle) {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let opts = FilletOptions {
        propagate_tangent: false,
        vertex_mode: popkern_fillet::types::VertexMode::Mitred,
        ..FilletOptions::default()
    };
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let r1 = fillet_edges(&cube, &specs, &opts, tol(), None).unwrap();
    // Plane-plane fillets on straight edges emit a Cylinder strip
    // (kernel convention — a torus would only arise on a curved
    // edge). The pre-fillet cube has no cylinder faces, so the match
    // uniquely identifies the blend. This bench asserted Torus and
    // panicked on every `cargo test --all-targets` run (PK-284).
    let blend = r1
        .brep
        .faces
        .keys()
        .find(|fh| {
            matches!(
                r1.brep.faces[*fh].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Cylinder { .. })
            )
        })
        .expect("post-fillet brep should contain the Cylinder blend strip");
    (r1.brep, blend)
}

/// Pre-fillet a closed plane-cylinder hole rim so the closed
/// benchmarks operate on the resulting annulus face.
fn closed_setup() -> Option<(BRep, FaceHandle)> {
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;
    use popkern_core::geometry::surface::SurfaceGeometry;
    let cube = BRep::unit_cube(tol()).unwrap();
    let cyl = BRep::cylinder(0.3, 2.0, tol()).unwrap();
    let translated = cyl
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .ok()?;
    let subtracted = boolean(BoolOperation::Subtraction, cube, translated, tol(), None).ok()?;
    let edge = subtracted.edges.keys().find(|e| {
        let he = &subtracted.half_edges[subtracted.edges[*e].half_edge];
        let fa = subtracted.loops[he.loop_].face;
        let fb = subtracted.loops[subtracted.half_edges[he.twin].loop_].face;
        let sa = subtracted.faces[fa].surface.as_ref().map(|s| &s.geometry);
        let sb = subtracted.faces[fb].surface.as_ref().map(|s| &s.geometry);
        matches!(
            (sa, sb),
            (
                Some(SurfaceGeometry::Plane { .. }),
                Some(SurfaceGeometry::Cylinder { .. })
            ) | (
                Some(SurfaceGeometry::Cylinder { .. }),
                Some(SurfaceGeometry::Plane { .. })
            )
        )
    })?;
    let r1 = fillet_edge(&subtracted, edge, 0.05, tol(), None).ok()?;
    let annulus = r1.brep.faces.keys().find(|fh| {
        let f = &r1.brep.faces[*fh];
        if f.inner_loops.len() != 1 {
            return false;
        }
        let single = |lh: popkern_core::topology::LoopHandle| -> bool {
            let Some(l) = r1.brep.loops.get(lh) else {
                return false;
            };
            let start = l.half_edge;
            let Some(he) = r1.brep.half_edges.get(start) else {
                return false;
            };
            he.next == start && he.prev == start
        };
        single(f.outer_loop) && single(f.inner_loops[0])
    })?;
    Some((r1.brep, annulus))
}

/// Locate the restored sharp edge in an unfilleted open brep by
/// matching against the captured plan endpoints.
fn find_open_restored_edge(
    unfilleted: &BRep,
    plan: &unfillet::UnfilletPlan,
) -> Option<popkern_core::topology::EdgeHandle> {
    let p0 = plan.restored_edge_start;
    let p1 = plan.restored_edge_end;
    let close = |a: &popkern_core::geometry::Point3, b: &popkern_core::geometry::Point3| -> bool {
        let dx = a.x - b.x;
        let dy = a.y - b.y;
        let dz = a.z - b.z;
        (dx * dx + dy * dy + dz * dz).sqrt() < 1e-6
    };
    unfilleted.edges.keys().find(|e| {
        let he = &unfilleted.half_edges[unfilleted.edges[*e].half_edge];
        let twin = &unfilleted.half_edges[he.twin];
        let v0 = unfilleted.vertices[he.origin].position;
        let v1 = unfilleted.vertices[twin.origin].position;
        (close(&v0, &p0) && close(&v1, &p1)) || (close(&v0, &p1) && close(&v1, &p0))
    })
}

fn bench_open_constant(c: &mut Criterion) {
    let (brep, blend) = open_setup();
    let mut group = c.benchmark_group("open/constant");

    group.bench_function("in_place", |b| {
        b.iter(|| {
            modify_inplace::modify_open_in_place(&brep, blend, 0.15).expect("constant in-place");
        });
    });

    group.bench_function("remove_recreate", |b| {
        b.iter(|| {
            let unfilleted = unfillet::unfillet_face(&brep, blend).unwrap();
            let plan = unfillet::compute_unfillet_plan(&brep, blend).unwrap();
            let restored = find_open_restored_edge(&unfilleted, &plan).unwrap();
            fillet_edge(&unfilleted, restored, 0.15, tol(), None).unwrap();
        });
    });

    group.finish();
}

fn bench_open_variable(c: &mut Criterion) {
    let (brep, blend) = open_setup();
    let mut group = c.benchmark_group("open/variable");

    group.bench_function("in_place", |b| {
        b.iter(|| {
            modify_inplace::modify_open_in_place_variable(&brep, blend, 0.05, 0.20)
                .expect("variable in-place");
        });
    });

    group.bench_function("remove_recreate", |b| {
        b.iter(|| {
            let unfilleted = unfillet::unfillet_face(&brep, blend).unwrap();
            let plan = unfillet::compute_unfillet_plan(&brep, blend).unwrap();
            let restored = find_open_restored_edge(&unfilleted, &plan).unwrap();
            let specs = vec![FilletSpec {
                edge: restored,
                radius: FilletRadius::Variable {
                    start: 0.05,
                    end: 0.20,
                },
            }];
            fillet_edges(&unfilleted, &specs, &FilletOptions::default(), tol(), None).unwrap();
        });
    });

    group.finish();
}

fn bench_open_bypoints(c: &mut Criterion) {
    let (brep, blend) = open_setup();
    let pts = vec![(0.0, 0.05), (0.5, 0.10), (1.0, 0.20)];
    let mut group = c.benchmark_group("open/bypoints");

    group.bench_function("in_place", |b| {
        b.iter(|| {
            modify_inplace::modify_open_in_place_bypoints(&brep, blend, &pts)
                .expect("bypoints in-place");
        });
    });

    group.bench_function("remove_recreate", |b| {
        b.iter(|| {
            let unfilleted = unfillet::unfillet_face(&brep, blend).unwrap();
            let plan = unfillet::compute_unfillet_plan(&brep, blend).unwrap();
            let restored = find_open_restored_edge(&unfilleted, &plan).unwrap();
            let specs = vec![FilletSpec {
                edge: restored,
                radius: FilletRadius::ByPoints(pts.clone()),
            }];
            fillet_edges(&unfilleted, &specs, &FilletOptions::default(), tol(), None).unwrap();
        });
    });

    group.finish();
}

fn bench_closed_constant(c: &mut Criterion) {
    let Some((brep, annulus)) = closed_setup() else {
        return;
    };
    let mut group = c.benchmark_group("closed/constant");

    group.bench_function("in_place", |b| {
        b.iter(|| {
            modify_inplace::modify_closed_in_place(&brep, annulus, 0.08).expect("closed in-place");
        });
    });

    group.bench_function("remove_recreate", |b| {
        use popkern_core::geometry::surface::SurfaceGeometry;
        b.iter(|| {
            let unfilleted = unfillet_closed::unfillet_closed_annulus(&brep, annulus).unwrap();
            let rim = unfilleted
                .edges
                .keys()
                .find(|e| {
                    let he = &unfilleted.half_edges[unfilleted.edges[*e].half_edge];
                    if he.origin != unfilleted.half_edges[he.twin].origin {
                        return false;
                    }
                    let fa = unfilleted.loops[he.loop_].face;
                    let fb = unfilleted.loops[unfilleted.half_edges[he.twin].loop_].face;
                    let sa = unfilleted.faces[fa].surface.as_ref().map(|s| &s.geometry);
                    let sb = unfilleted.faces[fb].surface.as_ref().map(|s| &s.geometry);
                    matches!(
                        (sa, sb),
                        (
                            Some(SurfaceGeometry::Plane { .. }),
                            Some(SurfaceGeometry::Cylinder { .. })
                        ) | (
                            Some(SurfaceGeometry::Cylinder { .. }),
                            Some(SurfaceGeometry::Plane { .. })
                        )
                    )
                })
                .unwrap();
            fillet_edge(&unfilleted, rim, 0.08, tol(), None).unwrap();
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_open_constant,
    bench_open_variable,
    bench_open_bypoints,
    bench_closed_constant,
);
criterion_main!(benches);
