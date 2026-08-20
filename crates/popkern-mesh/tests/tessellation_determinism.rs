//! Tessellation determinism pins (PK-277).
//!
//! The parallel per-face dispatch contract: `tessellate` output is
//! bit-identical call-to-call and at any thread count. The fixture is
//! the PK-295 served subtraction (rim-filleted cylinder minus a
//! strip-crossing cuboid) at a fine chord height — heavy enough to
//! drive encroachment rounds, the legacy fallback, boundary sewing,
//! and the crack zipper. This pin caught two real defects on landing:
//! hash-ordered Bowyer–Watson cavity refill (thread-dependent triangle
//! order) and hash-ordered crack-zipper candidate ties (call-to-call
//! triangle-count drift, 27304 vs 27333 on this fixture).

use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_mesh::{TessellationOptions, TriMesh};

fn served_subtraction() -> BRep {
    let tol = Tolerance::default();
    let cyl = BRep::cylinder(1.0, 2.0, tol).expect("cylinder builds");
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
    let part = fillet_edges(&cyl, &rim, &FilletOptions::default(), tol, None)
        .expect("rim fillet completes")
        .brep;
    let cutter = BRep::unit_cube(tol)
        .expect("cube builds")
        .transform(&Transform3::scale_uniform(1.8))
        .expect("scales")
        .transform(&Transform3::translation(-0.9, -0.9, 0.6))
        .expect("translates");
    popkern_bool::subtraction(part, cutter, tol, None).expect("subtraction succeeds")
}

fn options() -> TessellationOptions {
    TessellationOptions {
        chord_height: 0.0005,
        ..Default::default()
    }
}

fn assert_meshes_identical(lhs: &TriMesh, rhs: &TriMesh, label: &str) {
    assert_eq!(
        lhs.triangles, rhs.triangles,
        "{label}: triangle lists differ"
    );
    assert_eq!(
        lhs.vertices.len(),
        rhs.vertices.len(),
        "{label}: vertex counts differ"
    );
    for (i, (l, r)) in lhs.vertices.iter().zip(rhs.vertices.iter()).enumerate() {
        let bits =
            |p: &popkern_core::geometry::Point3| (p.x.to_bits(), p.y.to_bits(), p.z.to_bits());
        assert_eq!(
            bits(&l.position),
            bits(&r.position),
            "{label}: vertex {i} position differs"
        );
    }
}

/// Repeat calls on the same solid must be bit-identical (pins the
/// crack-zipper canonical candidate order).
#[test]
fn repeat_tessellation_is_bit_identical() {
    let solid = served_subtraction();
    let m1 = popkern_mesh::tessellate(&solid, &options()).expect("first tessellation");
    let m2 = popkern_mesh::tessellate(&solid, &options()).expect("second tessellation");
    assert_meshes_identical(&m1, &m2, "repeat call");
}

/// 1-thread and N-thread pools must produce bit-identical meshes
/// (pins the ordered parallel dispatch and the canonical Delaunay
/// cavity fill).
#[test]
fn thread_count_does_not_change_served_mesh() {
    let solid = served_subtraction();
    let single = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .expect("1-thread pool builds")
        .install(|| popkern_mesh::tessellate(&solid, &options()).expect("tessellation"));
    let many = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .expect("8-thread pool builds")
        .install(|| popkern_mesh::tessellate(&solid, &options()).expect("tessellation"));
    assert_meshes_identical(&single, &many, "1-thread vs 8-thread");
}
