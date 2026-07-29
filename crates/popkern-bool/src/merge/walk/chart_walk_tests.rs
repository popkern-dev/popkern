//! Convention pins for the chart-subdivision walk (PK-104).
//!
//! These fix the rotation rule's handedness by construction: the
//! boundary-T alternation case (where the legacy 3D angular pick
//! measurably chose straight-through over the trim) and the orbit
//! partition property. The full-corpus behavior is covered by the
//! boolean and part suites plus the breached-bore regression gate;
//! these tests exist so a sign error in the angle formula fails a
//! unit test, not a fixture hunt.

use super::chart_subdivision::build_face_subdivision;
use super::chart_walk::{canonical_loops, compare_walks, walk_chart_orbits};
use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph, NodeSource, Operand};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;

fn xy_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

/// A unit square face split by a vertical trim through the midpoints
/// of the bottom and top edges. The subdivision has six boundary
/// half-edges (bottom/top split at the trim feet) plus both trim
/// directions, and the correct partition is exactly two four-loops.
fn square_with_transversal_trim() -> (IntersectionGraph, FaceRef, Vec<usize>) {
    let mut g = IntersectionGraph::new();
    let face = FaceRef {
        operand: Operand::A,
        face: popkern_core::topology::FaceHandle::default(),
    };
    let eh = popkern_core::topology::EdgeHandle::default();
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    let n = [
        g.add_node(p(0.0, 0.0), NodeSource::Intersection), // 0 SW
        g.add_node(p(0.5, 0.0), NodeSource::Intersection), // 1 bottom T
        g.add_node(p(1.0, 0.0), NodeSource::Intersection), // 2 SE
        g.add_node(p(1.0, 1.0), NodeSource::Intersection), // 3 NE
        g.add_node(p(0.5, 1.0), NodeSource::Intersection), // 4 top T
        g.add_node(p(0.0, 1.0), NodeSource::Intersection), // 5 NW
    ];
    // Boundary CCW: only the face-side directions are registered.
    let mut hes = Vec::new();
    for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)] {
        let (fwd, _) = g.add_half_edge_pair(n[a], n[b], EdgeSource::OriginalA(eh), None);
        g.set_half_edge_face(fwd, face);
        hes.push(fwd);
    }
    // The trim: both directions registered (cross-registration).
    let (t_up, t_down) = g.add_half_edge_pair(
        n[1],
        n[4],
        EdgeSource::Shared {
            face_a: popkern_core::topology::FaceHandle::default(),
            face_b: popkern_core::topology::FaceHandle::default(),
        },
        None,
    );
    g.set_half_edge_face(t_up, face);
    g.set_half_edge_face(t_down, face);
    hes.push(t_up);
    hes.push(t_down);
    (g, face, hes)
}

#[test]
fn transversal_trim_yields_two_four_loops() {
    let (g, face, hes) = square_with_transversal_trim();
    let sub = build_face_subdivision(
        &g,
        &face,
        &hes,
        &xy_plane(),
        true,
        &[],
        &BRep::new(popkern_core::tolerance::Tolerance::default()),
    );
    assert!(
        sub.complete(),
        "substrate must pass I1–I5: {:?}",
        sub.report.issues
    );
    let orbits = walk_chart_orbits(&sub, &g).expect("orbits");
    assert_eq!(orbits.len(), 2, "one trim splits the square in two");
    assert!(orbits.iter().all(|o| o.len() == 4));
    // The two fragments consume opposite trim directions and never
    // share a half-edge (partition property).
    let all: Vec<usize> = orbits.iter().flatten().copied().collect();
    let set: std::collections::BTreeSet<usize> = all.iter().copied().collect();
    assert_eq!(all.len(), set.len(), "orbits partition the half-edges");
}

/// The boundary-T alternation pin: arriving at the trim foot along
/// the bottom-left boundary, the material-left walk MUST turn onto
/// the trim (the legacy 3D angular pick measurably continued
/// straight past it — the mega-loop class). With material on the
/// left and orientation=true, the west fragment is
/// [SW→T, T→top, top→NW, NW→SW].
#[test]
fn walk_turns_onto_trim_at_boundary_t() {
    let (g, face, hes) = square_with_transversal_trim();
    let sub = build_face_subdivision(
        &g,
        &face,
        &hes,
        &xy_plane(),
        true,
        &[],
        &BRep::new(popkern_core::tolerance::Tolerance::default()),
    );
    let orbits = walk_chart_orbits(&sub, &g).expect("orbits");
    // hes[0] = bottom-left (SW→T), hes[6] = trim up (T→top).
    let west = orbits
        .iter()
        .find(|o| o.contains(&hes[0]))
        .expect("west fragment");
    assert!(
        west.contains(&hes[6]),
        "the walk arriving at the T must turn onto the trim (up), got {west:?}"
    );
    assert!(
        !west.contains(&hes[1]),
        "the west fragment must NOT continue straight past the trim"
    );
}

/// Orientation=false mirrors the turn: the same graph walked with a
/// flipped chart pairs the bottom-left boundary with the OTHER trim
/// direction's region.
#[test]
fn orientation_false_mirrors_the_partition() {
    let (g, face, hes) = square_with_transversal_trim();
    let sub_t = build_face_subdivision(
        &g,
        &face,
        &hes,
        &xy_plane(),
        true,
        &[],
        &BRep::new(popkern_core::tolerance::Tolerance::default()),
    );
    let sub_f = build_face_subdivision(
        &g,
        &face,
        &hes,
        &xy_plane(),
        false,
        &[],
        &BRep::new(popkern_core::tolerance::Tolerance::default()),
    );
    // σ self-calibrates from the stored boundary winding, so BOTH
    // orientations resolve to the same (correct) handedness — the
    // flag is only the fallback when no boundary exists.
    assert_eq!(sub_t.sigma, sub_f.sigma, "σ is boundary-calibrated");
    let o_t = walk_chart_orbits(&sub_t, &g).expect("orbits");
    let o_f = walk_chart_orbits(&sub_f, &g).expect("orbits");
    assert_eq!(canonical_loops(&o_t), canonical_loops(&o_f));
}

#[test]
fn compare_walks_reports_divergence_classes() {
    let legacy = vec![vec![1usize, 2, 3], vec![7, 8]];
    let orbits = vec![vec![2usize, 3, 1], vec![9, 10]];
    let diff = compare_walks(&legacy, &orbits);
    assert_eq!(diff.matched, 1, "rotation-invariant loop match");
    assert_eq!(diff.legacy_only, vec![vec![7, 8]]);
    assert_eq!(diff.orbit_only, vec![vec![9, 10]]);
    assert_eq!(diff.contested(), 0, "disjoint extras are not contested");
}
