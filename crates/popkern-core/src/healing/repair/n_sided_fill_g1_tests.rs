//! Gates for the G¹ hole fill.
//!
//! These run against a real `BRep` rather than synthetic point rings,
//! because the whole job of this module is reading frames out of
//! topology: curves, neighbouring surfaces, and edge orientation that
//! hole detection does not normalise.

use super::*;
use crate::geometry::NSidedFillError;
use crate::healing::repair::hole_fill::DetectedHole;
use crate::tolerance::Tolerance;

/// The four edges bounding one face of a unit cube, in loop order —
/// the shape a deleted face leaves behind.
fn cube_face_boundary(brep: &BRep) -> Vec<EdgeHandle> {
    let (face, _) = brep.faces.iter().next().expect("cube has faces");
    let loop_handle = brep.faces[face].outer_loop;
    let start = brep.loops[loop_handle].half_edge;
    let mut edges = Vec::new();
    let mut he = start;
    loop {
        let data = &brep.half_edges[he];
        edges.push(data.edge);
        he = data.next;
        if he == start {
            break;
        }
    }
    edges
}

fn cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).expect("unit cube fixture")
}

#[test]
fn fills_a_cube_face_opening_with_a_g1_cap() {
    let brep = cube();
    let boundary_edges = cube_face_boundary(&brep);
    assert_eq!(boundary_edges.len(), 4, "a cube face is bounded by 4 edges");
    let hole = DetectedHole {
        num_sides: boundary_edges.len(),
        boundary_edges,
    };

    let fill = fill_n_sided_hole_g1(&brep, &hole, 1e-6)
        .expect("a square opening with planar surroundings must fill");

    assert_eq!(fill.surfaces.len(), 4, "one patch per side");
    assert!(
        fill.worst_boundary_gap < 1e-9,
        "the cap must interpolate the opening it was built from, worst gap {:.3e}",
        fill.worst_boundary_gap
    );
}

#[test]
fn reversed_boundary_edges_are_reoriented_not_rejected() {
    // Hole detection returns edges in adjacency order but each edge
    // keeps its own stored direction, so about half of them run
    // backwards relative to the loop. Reversing the list must still
    // fill: if orientation were not normalised the ring would read as
    // open and decline on a boundary that is in fact closed.
    let brep = cube();
    let mut boundary_edges = cube_face_boundary(&brep);
    boundary_edges.reverse();
    let hole = DetectedHole {
        num_sides: boundary_edges.len(),
        boundary_edges,
    };
    let fill = fill_n_sided_hole_g1(&brep, &hole, 1e-6)
        .expect("a reversed boundary list is still a closed loop");
    assert!(fill.worst_boundary_gap < 1e-9);
}

#[test]
fn a_hole_with_too_many_sides_declines_through_the_primitive() {
    // 13 sides is past the bi-4 cap's range; the decline must come
    // through as the primitive's own typed error, not a generic one.
    let brep = cube();
    let edges = cube_face_boundary(&brep);
    let mut boundary_edges = Vec::new();
    while boundary_edges.len() < 13 {
        boundary_edges.extend(edges.iter().copied());
    }
    boundary_edges.truncate(13);
    let hole = DetectedHole {
        num_sides: boundary_edges.len(),
        boundary_edges,
    };
    match fill_n_sided_hole_g1(&brep, &hole, 1e-6) {
        Err(G1HoleFillError::Fill(NSidedFillError::UnsupportedSideCount { got })) => {
            assert_eq!(got, 13);
        }
        // Repeating edges can also break the loop-order check first;
        // either decline is honest, but it must be typed and specific.
        Err(G1HoleFillError::BoundaryOutOfOrder { .. }) => {}
        other => panic!("an over-long boundary must decline specifically, got {other:?}"),
    }
}

#[test]
fn a_face_without_surface_geometry_declines_naming_the_edge() {
    let mut brep = cube();
    let boundary_edges = cube_face_boundary(&brep);
    // Strip the surface from every face so no tangent plane exists.
    let faces: Vec<_> = brep.faces.iter().map(|(f, _)| f).collect();
    for face in faces {
        brep.faces[face].surface = None;
    }
    let hole = DetectedHole {
        num_sides: boundary_edges.len(),
        boundary_edges,
    };
    match fill_n_sided_hole_g1(&brep, &hole, 1e-6) {
        Err(G1HoleFillError::MissingNeighbourSurface { index }) => {
            assert_eq!(index, 0, "the first edge is where it should be noticed");
        }
        other => panic!("a surfaceless neighbour must decline explicitly, got {other:?}"),
    }
}

/// The G¹ cap is a coverage gain, not only a continuity gain.
///
/// Measured on a unit-cube face opening: the incumbent Coons fan
/// returns `None` — every one of its four sub-patches fails to
/// construct, because each is built with a *degenerate zero-length
/// top edge* at the centroid — while the G¹ cap fills the same
/// opening and interpolates it exactly. So the two strategies are not
/// interchangeable: on this class the G¹ path is the only one that
/// produces anything at all.
///
/// This gate pins that relationship. If the Coons fan is later fixed
/// to serve this opening, this test should be updated to compare the
/// two rather than record the gap.
#[test]
fn g1_cap_serves_an_opening_the_coons_fan_declines() {
    let brep = cube();
    let boundary_edges = cube_face_boundary(&brep);
    let hole = DetectedHole {
        num_sides: boundary_edges.len(),
        boundary_edges,
    };
    let g1 = fill_n_sided_hole_g1(&brep, &hole, 1e-6);
    assert!(g1.is_ok(), "G1 cap must serve this opening: {:?}", g1.err());
    assert_eq!(
        g1.expect("checked above").surfaces.len(),
        4,
        "one patch per side"
    );

    let coons = super::super::n_sided_fill::fill_n_sided_hole(&brep, &hole, 1e-6);
    assert!(
        coons.is_none(),
        "recorded gap: the Coons fan declines this opening (its sub-patches carry a \
         zero-length top edge). If it now serves it, update this gate to a comparison."
    );
}
