//! Blend×blend corner seam acceptance (PK-152 stage 1).
//!
//! Two fillets on adjacent cube edges meeting at a valence-2 blend
//! corner (third edge sharp) must close by trimming both strips to
//! their mutual intersection curve. Design:
//! `docs/architecture/blend-blend-intersection-trimming.md`.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "support/mesh_gauges.rs"]
mod mesh_gauges;
use mesh_gauges::welded_open_edges;

fn default_tol() -> Tolerance {
    Tolerance::default()
}

/// Two cube edges sharing exactly one vertex.
fn adjacent_edge_pair(brep: &BRep) -> (EdgeHandle, EdgeHandle) {
    for (e1, d1) in brep.edges.iter() {
        for (e2, d2) in brep.edges.iter() {
            if e1 == e2 {
                continue;
            }
            let ends = |d: &popkern_core::topology::Edge| {
                let h = &brep.half_edges[d.half_edge];
                (h.origin, brep.half_edges[h.twin].origin)
            };
            let (a1, b1) = ends(d1);
            let (a2, b2) = ends(d2);
            if a1 == a2 || a1 == b2 || b1 == a2 || b1 == b2 {
                return (e1, e2);
            }
        }
    }
    panic!("cube has adjacent edges");
}

fn pair_body(r1: f64, r2: f64) -> BRep {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let (e1, e2) = adjacent_edge_pair(&brep);
    fillet_edges(
        &brep,
        &[
            FilletSpec {
                edge: e1,
                radius: FilletRadius::Constant(r1),
            },
            FilletSpec {
                edge: e2,
                radius: FilletRadius::Constant(r2),
            },
        ],
        &FilletOptions::default(),
        default_tol(),
        None,
    )
    .expect("adjacent pair completes")
    .brep
}

/// Current measured state: the pair-corner body validates but ~5
/// corner edges run through the old cube corner at (√2−1)·r off both
/// strips, and the welded gauge reads ~58 opens. This test PINS that
/// measurement so the defect cannot silently grow while stage 1 is
/// built; the completion gates below flip when the seam lands.
#[test]
fn equal_radius_pair_corner_pinned_measurement() {
    let out = pair_body(0.08, 0.08);
    out.validate()
        .expect("pair body validates (measured today)");
    let worst = out
        .edges
        .keys()
        .map(|eh| out.measure_edge_deviation(eh))
        .fold(0.0_f64, f64::max);
    let stale = (2.0_f64.sqrt() - 1.0) * 0.08;
    assert!(
        worst <= stale * 1.1,
        "corner-edge staleness must not grow past the (√2−1)·r class: {worst:.3e}"
    );
    assert!(
        worst >= stale * 0.9,
        "corner staleness vanished ({worst:.3e}) — the seam landed? Flip the \
         completion gates below and retire this pin."
    );
}

#[test]
#[ignore = "PK-152 stage 1: flip when the blend-x-blend corner seam lands — equal-radius \
            perpendicular pair closes welded-watertight with every edge on its surfaces"]
fn equal_radius_pair_corner_closes_watertight() {
    let out = pair_body(0.08, 0.08);
    out.validate().expect("pair body validates");
    assert_eq!(welded_open_edges(&out), 0, "pair corner must weld shut");
    let worst = out
        .edges
        .keys()
        .map(|eh| out.measure_edge_deviation(eh))
        .fold(0.0_f64, f64::max);
    assert!(
        worst < 1e-6,
        "every edge must lie on its faces' surfaces: worst {worst:.3e}"
    );
}

#[test]
#[ignore = "PK-152 stage 2: flip when the marching-intersector seam serves mismatched radii"]
fn mixed_radius_pair_corner_closes_watertight() {
    let out = pair_body(0.06, 0.10);
    out.validate().expect("mixed pair body validates");
    assert_eq!(
        welded_open_edges(&out),
        0,
        "mixed pair corner must weld shut"
    );
}
