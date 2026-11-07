use super::*;
use crate::topology::FaceSide;
use crate::topology::radial_order::book_test_helpers::{DEG, book, tol};

/// Pages linked out of order come back in angular order about the
/// spine, starting at its first fin, whichever way each fin runs.
#[test]
fn a_book_orders_its_pages_by_angle_about_the_spine() {
    let angles = [0.0, 220.0 * DEG, 100.0 * DEG, 300.0 * DEG];
    let (brep, spine, fins) = book(&angles, &[true, false, true, false]);
    let order = brep.radial_order(spine, &tol()).expect("ordered");
    assert_eq!(order, vec![fins[0], fins[2], fins[1], fins[3]]);
}

/// Ordering the ring relinks it, and every wedge then spans two
/// adjacent pages, the side of each facing it read off the fin sense.
#[test]
fn ordering_the_ring_relinks_it_in_order() {
    let angles = [0.0, 200.0 * DEG, 90.0 * DEG];
    let (mut brep, spine, fins) = book(&angles, &[true, true, false]);
    brep.order_ring(spine, &tol()).expect("ordered");
    assert_eq!(brep.radial_fins(spine), vec![fins[0], fins[2], fins[1]]);
    let wedges = brep.edge_wedges(spine).expect("decidable");
    assert_eq!(wedges.len(), 3);
    // The wedge after the forward first page lies on its front; the
    // reversed page at 90° bounds it with its front too.
    assert_eq!((wedges[0].from, wedges[0].to), (fins[0], fins[2]));
    assert_eq!(wedges[0].from_side, FaceSide::Front);
    assert_eq!(wedges[0].to_side, FaceSide::Front);
}

/// Two pages in one plane leave the spine together and nothing tells
/// them apart: a typed decline naming the tie, the ring unchanged.
#[test]
fn coincident_pages_decline_with_their_evidence() {
    let angles = [0.0, 90.0 * DEG, 90.0 * DEG];
    let (mut brep, spine, fins) = book(&angles, &[true, true, true]);
    let before = brep.radial_fins(spine);
    match brep.order_ring(spine, &tol()) {
        Err(RadialOrderError::Declined {
            edge,
            fins: listed,
            decline:
                RadialOrderDecline::Ambiguous {
                    fins: tied,
                    evidence,
                    ..
                },
        }) => {
            assert_eq!(edge, spine);
            assert_eq!(listed, fins);
            let mut tied = tied;
            tied.sort_unstable();
            assert_eq!(tied, vec![1, 2]);
            assert_eq!(evidence.len(), 2);
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(brep.radial_fins(spine), before);
}

/// Edges of one or two fins are their own order; an unknown edge and
/// a fin on a surfaceless face are typed errors.
#[test]
fn short_rings_unknown_edges_and_bare_faces() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    for edge in cube.edges.keys() {
        assert_eq!(
            cube.radial_order(edge, &tol()).expect("two fins"),
            cube.radial_fins(edge)
        );
    }
    let (mut brep, spine, fins) = book(&[0.0, 1.0, 2.0], &[true, true, true]);
    let mut other = brep.clone();
    other.edges.remove(spine);
    assert_eq!(
        other.radial_order(spine, &tol()),
        Err(RadialOrderError::UnknownEdge(spine))
    );
    let face = brep.fin_face(fins[1]).expect("face");
    brep.faces[face].surface = None;
    assert_eq!(
        brep.radial_order(spine, &tol()),
        Err(RadialOrderError::NoSurface {
            edge: spine,
            fin: fins[1]
        })
    );
}

/// The geometric sense reading agrees with topology on every fin of a
/// cylinder's closed rims and seam, and on a cube's edges: the rule a
/// closed edge of more than two fins is read by.
#[test]
fn the_pcurve_sense_agrees_with_topology() {
    for body in [
        BRep::cylinder(0.5, 2.0, tol()).expect("cylinder"),
        BRep::unit_cube(tol()).expect("cube"),
    ] {
        for (edge, _) in &body.edges {
            for fin in body.radial_fins(edge) {
                let topology = body.fin_is_forward(fin);
                assert!(topology.is_some());
                assert_eq!(body.fin_sense_by_pcurve(fin), topology, "{edge:?} {fin:?}");
                assert_eq!(body.fin_sense(fin), topology);
            }
        }
    }
}
