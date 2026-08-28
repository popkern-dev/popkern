//! PK-339: sphere-cap and full-sphere revolution charts weld by
//! construction. The full-wrap revolution path used to decline any
//! pole chart whose meridian seam carried interior stations (curved
//! seams chord-refine — "cross chain not at constant v"), dropping
//! sphere caps and full spheres to uniform-grid sampling whose pole
//! fans duplicate (welded census: plain sphere 27 open edges,
//! sphere_capped_rod 16, PK-326 apex cap 5). The meridian-ladder
//! extension takes the seam chain's stations AS the latitude rows —
//! the face boundary equals the ring exactly — with a two-apex arm
//! for full spheres.
//!
//! The census here is STRICT: every welded edge must be shared by
//! EXACTLY two triangles (c == 1 misses the ×4 pole-fan overlaps the
//! uniform-grid path emitted).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::{TessellationOptions, tessellate};
use popkern_shell::ShellOptions;
use popkern_shell::offset::offset_solid;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Welded edges not shared by exactly two triangles (quantized
/// positions — the TriMesh never dedups vertices, PK-338 lesson).
fn non_manifold_welded_edges(brep: &BRep) -> usize {
    let mesh = tessellate(brep, &TessellationOptions::default()).expect("tessellate");
    let q = |x: f64| (x * 1.0e7).round() as i64;
    let key = |i: usize| {
        let p = mesh.vertices[i].position;
        (q(p.x), q(p.y), q(p.z))
    };
    let mut counts = std::collections::HashMap::new();
    for t in &mesh.triangles {
        let w = [key(t[0]), key(t[1]), key(t[2])];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0_i32) += 1;
        }
    }
    counts.values().filter(|&&c| c != 2).count()
}

#[test]
fn sphere_class_charts_weld_by_construction() {
    let sphere = BRep::sphere(0.4, tol()).unwrap();
    assert_eq!(
        non_manifold_welded_edges(&sphere),
        0,
        "full sphere (bipolar revolution arm)"
    );

    let rod = popkern_examples::offset_corpus::sphere_capped_rod(tol())
        .unwrap()
        .brep;
    assert_eq!(
        non_manifold_welded_edges(&rod),
        0,
        "sphere_capped_rod caps (meridian-ladder arm)"
    );

    let cone = BRep::cone(1.0, 1.5, tol()).unwrap();
    let capped = offset_solid(&cone, 0.1, &ShellOptions::default(), tol()).unwrap();
    assert_eq!(
        non_manifold_welded_edges(&capped.brep),
        0,
        "PK-326 apex sphere cap"
    );
}

#[test]
fn sphere_class_charts_need_no_zipper_repair() {
    for (name, brep) in [
        ("sphere", BRep::sphere(0.4, tol()).unwrap()),
        (
            "sphere_capped_rod",
            popkern_examples::offset_corpus::sphere_capped_rod(tol())
                .unwrap()
                .brep,
        ),
    ] {
        let mesh = tessellate(&brep, &TessellationOptions::default()).expect("tessellate");
        assert_eq!(
            mesh.zipped_t_junctions, 0,
            "{name}: welds must close by construction, not zipper repair"
        );
    }
}
