//! Blend×blend corner seam acceptance (PK-152).
//!
//! Two fillets on adjacent cube edges meeting at a valence-2 blend
//! corner (third edge sharp) must close by trimming both strips to
//! their mutual intersection curve. Design:
//! `docs/architecture/blend-blend-intersection-trimming.md`.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};
use popkern_fillet::{BlendResult, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

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

fn pair_result(r1: f64, r2: f64) -> BlendResult {
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
}

/// PK-152 stage 1 HARD gate: the equal-radius perpendicular pair
/// corner closes welded-watertight with every edge on its surfaces.
/// The blend-seam-corner settlement retires the phantom flap and
/// re-authors the shared edge as the exact intersection ellipse
/// (formerly the `(√2−1)·r` stale-corner class: worst deviation
/// 3.314e-2 and 58 welded opens at r = 0.08).
#[test]
fn equal_radius_pair_corner_closes_watertight() {
    let result = pair_result(0.08, 0.08);
    let report = &result.repair_report.blend_seam_corner;
    assert_eq!(report.seams_built, 1, "pair corner must build one seam");
    assert!(
        report.declines.is_empty(),
        "no seam declines expected: {:?}",
        report.declines
    );
    let out = result.brep;
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

/// PK-152 stage 2 HARD gate: the mismatched-radius pair corner
/// closes welded-watertight. The seam is the marched cylinder×
/// cylinder quartic from the rail crossing to the smaller strip's
/// rail exit; the near anchor moves onto the exit point and the Δr
/// mismatch stub becomes the larger strip's exact cylinder×plane
/// ellipse arc (formerly 58 welded opens, worst 4.142e-2).
#[test]
fn mixed_radius_pair_corner_closes_watertight() {
    let result = pair_result(0.06, 0.10);
    let report = &result.repair_report.blend_seam_corner;
    assert_eq!(report.seams_built, 1, "mixed corner must build one seam");
    assert!(
        report.declines.is_empty(),
        "no seam declines expected: {:?}",
        report.declines
    );
    let out = result.brep;
    out.validate().expect("mixed pair body validates");
    assert_eq!(
        welded_open_edges(&out),
        0,
        "mixed pair corner must weld shut"
    );
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
