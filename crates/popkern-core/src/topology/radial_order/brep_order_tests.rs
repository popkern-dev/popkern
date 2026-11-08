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

/// An edge whose first fin runs its stored curve backward (an operator
/// re-anchored the ring without reversing the curve) is ordered about
/// that fin's direction (PK-1118): the cycle comes back reversed, as
/// seen from the other end, and every wedge spans the same two pages
/// with the same sides. Ordered about the curve instead, the ring read
/// one solid's material against another's.
#[test]
fn an_edge_whose_first_fin_runs_backward_orders_about_that_fin() {
    let angles = [0.0, 220.0 * DEG, 100.0 * DEG, 300.0 * DEG];
    // The spine's stored curve runs x from 0 to 1, as page 0 does.
    let line = |b: &mut BRep, spine| {
        b.edges[spine].curve = Some(crate::geometry::curve::CurveOnEdge {
            geometry: crate::geometry::curve::CurveGeometry::Line {
                origin: crate::geometry::Point3::new(0.0, 0.0, 0.0),
                direction: crate::geometry::Vec3::new(1.0, 0.0, 0.0),
            },
            t_start: 0.0,
            t_end: 1.0,
        });
    };
    let (mut forward, spine, fins) = book(&angles, &[true, false, true, false]);
    line(&mut forward, spine);
    forward.order_ring(spine, &tol()).expect("ordered");
    let sides = |b: &BRep| {
        let mut out: Vec<_> = b
            .edge_wedges(spine)
            .expect("decidable")
            .into_iter()
            .map(|w| {
                // Each wedge as its two pages and the side each shows it,
                // whichever way round it was walked.
                let (a, b) = ((w.from, w.from_side), (w.to, w.to_side));
                if a.0 < b.0 { (a, b) } else { (b, a) }
            })
            .collect();
        out.sort_by_key(|((f, _), (g, _))| (*f, *g));
        out
    };
    let want = sides(&forward);
    let (mut backward, _, _) = book(&angles, &[true, false, true, false]);
    line(&mut backward, spine);
    // The reversed page made first, the curve as it was.
    backward.edges[spine].half_edge = fins[1];
    let order = backward.radial_order(spine, &tol()).expect("ordered");
    assert_eq!(order, vec![fins[1], fins[2], fins[0], fins[3]]);
    backward.order_ring(spine, &tol()).expect("ordered");
    assert_eq!(sides(&backward), want);
}
