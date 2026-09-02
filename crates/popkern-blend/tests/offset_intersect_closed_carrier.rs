//! PK-421: the offset-intersect centre-curve path must serve a valid
//! sweep for an edge-span window on a CLOSED carrier.
//!
//! A prism boss through a cylinder wall, centred on the wall's seam,
//! leaves plane × cylinder chain pieces that are short exact arcs
//! ending at the seam vertex. Their offset carrier is the full offset
//! circle, and the piece's window straddles that carrier's parameter
//! seam. With the spine walker forced off (the production route when
//! the walker declines), the sweep channel must still section only
//! the piece's own span — measured before the fix: the seam-crossing
//! piece's sections ran the long way round the carrier and the sweep
//! twisted (`SelfIntersectingBlend`, normal field flipped between
//! adjacent samples).

use popkern_blend::dispatch_overrides::WALKER_DISABLE;
use popkern_blend::{BlendFailure, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

/// The seam-straddling boss-on-wall union and its plane × cylinder
/// chain pieces (edges on the wall between the caps).
fn boss_on_wall_chain(tol: Tolerance) -> (BRep, Vec<EdgeHandle>) {
    let cyl = BRep::cylinder(1.0, 2.0, tol).expect("cylinder");
    let outline = [
        Point3::new(0.6, -0.25, 0.75),
        Point3::new(1.3, -0.25, 0.75),
        Point3::new(1.3, 0.25, 0.75),
        Point3::new(0.6, 0.25, 0.75),
    ];
    let boss = BRep::prism(&outline, Vec3::new(0.0, 0.0, 1.0), 0.5, tol).expect("prism");
    let united = popkern_bool::union(cyl, boss, tol, None).expect("union");
    let kind = |brep: &BRep, f| {
        brep.faces[f]
            .surface
            .as_ref()
            .map(|s| matches!(s.geometry, SurfaceGeometry::Cylinder { .. }))
    };
    let mut chain = Vec::new();
    for (e, d) in united.edges.iter() {
        let he = &united.half_edges[d.half_edge];
        let f1 = united.loops[he.loop_].face;
        let f2 = united.loops[united.half_edges[he.twin].loop_].face;
        let a = united.vertices[he.origin].position;
        if f1 != f2
            && kind(&united, f1) != kind(&united, f2)
            && kind(&united, f1).is_some()
            && kind(&united, f2).is_some()
            && a.z > 0.5
            && a.z < 1.5
        {
            chain.push(e);
        }
    }
    chain.sort();
    (united, chain)
}

/// Every chain piece's geometry computes on the offset-intersect
/// path: no piece declines with a self-intersecting (twisted) sweep.
#[test]
fn a_seam_straddling_arc_piece_sweeps_untwisted_on_the_offset_intersect_path() {
    let tol = Tolerance::default();
    let (united, chain) = boss_on_wall_chain(tol);
    assert_eq!(chain.len(), 6, "six chain pieces");
    let _forced = WALKER_DISABLE.set_guard(true);
    let specs: Vec<FilletSpec> = chain
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.05),
        })
        .collect();
    let twisted: Vec<String> =
        match fillet_edges(&united, &specs, &FilletOptions::default(), tol, None) {
            Ok(_) => Vec::new(),
            Err(popkern_blend::BlendError::PartialFailure { failures, .. }) => failures
                .iter()
                .filter_map(|(edge, failure)| match failure {
                    BlendFailure::SelfIntersectingBlend { reason, .. } => {
                        Some(format!("{edge:?}: {reason}"))
                    }
                    _ => None,
                })
                .collect(),
            Err(other) => panic!("unexpected decline: {other:?}"),
        };
    assert!(
        twisted.is_empty(),
        "twisted sweeps on the offset-intersect path: {twisted:?}"
    );
}
