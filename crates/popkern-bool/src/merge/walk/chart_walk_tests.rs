//! Convention pins for the chart-subdivision walk (PK-104).
//!
//! These fix the rotation rule's handedness by construction: the
//! boundary-T alternation case (where the 3D-Newell angular pick
//! measurably chose straight-through over the trim) and the orbit
//! partition property. The full-corpus behavior is covered by the
//! boolean and part suites plus the breached-bore regression gate;
//! these tests exist so a sign error in the angle formula fails a
//! unit test, not a fixture hunt.

use std::collections::BTreeSet;

use super::chart_subdivision::build_face_subdivision;
use super::chart_walk::walk_chart_orbits;
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
    square_with_transversal_trim_wound(true)
}

/// The same square, its boundary registered CCW in the chart or CW
/// (as a face with `orientation == false` winds under the PK-772
/// contract). CW reverses each boundary half-edge and the cycle;
/// `hes[0]` is then bottom-left T→SW and `hes[1]` bottom-right SE→T.
fn square_with_transversal_trim_wound(ccw: bool) -> (IntersectionGraph, FaceRef, Vec<usize>) {
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
    // Only the face-side directions are registered.
    let mut hes = Vec::new();
    for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)] {
        let (a, b) = if ccw { (a, b) } else { (b, a) };
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
        &popkern_core::tolerance::Tolerance::default(),
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
/// the trim (the 3D-Newell angular pick measurably continued
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
        &popkern_core::tolerance::Tolerance::default(),
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

/// Orientation=false mirrors the turn. Under the PK-772 contract such
/// a face's boundary winds CW in the chart, σ is the flag (−1), and
/// the walk, material on the left in the V-mirrored frame, partitions
/// the square into the same two halves: the west half takes T→SW and
/// turns onto the trim's DOWN direction, never past it to SE→T.
#[test]
fn orientation_false_mirrors_the_partition() {
    let (g, face, hes) = square_with_transversal_trim_wound(false);
    let sub = build_face_subdivision(
        &g,
        &face,
        &hes,
        &xy_plane(),
        false,
        &[],
        &BRep::new(popkern_core::tolerance::Tolerance::default()),
        &popkern_core::tolerance::Tolerance::default(),
    );
    // No stored ring reads on this synthetic face, so σ is the flag.
    assert_eq!(sub.sigma, -1.0, "σ is the orientation sign");
    assert!(sub.complete(), "I1–I5: {:?}", sub.report.issues);
    let orbits = walk_chart_orbits(&sub, &g).expect("orbits");
    assert_eq!(orbits.len(), 2, "one trim splits the square in two");
    assert!(orbits.iter().all(|o| o.len() == 4));
    // hes[0] = T→SW, hes[1] = SE→T, hes[7] = trim down (top→T).
    let west = orbits
        .iter()
        .find(|o| o.contains(&hes[0]))
        .expect("west fragment");
    assert!(west.contains(&hes[7]), "west takes the trim down: {west:?}");
    assert!(!west.contains(&hes[1]), "west never runs past the trim");
}

#[test]
fn compare_walks_reports_divergence_classes() {
    let newell = vec![vec![1usize, 2, 3], vec![7, 8]];
    let orbits = vec![vec![2usize, 3, 1], vec![9, 10]];
    let diff = compare_walks(&newell, &orbits);
    assert_eq!(diff.matched, 1, "rotation-invariant loop match");
    assert_eq!(diff.newell_only, vec![vec![7, 8]]);
    assert_eq!(diff.orbit_only, vec![vec![9, 10]]);
    assert_eq!(diff.contested(), 0, "disjoint extras are not contested");
}

/// Canonical form for loop-set comparison: each loop rotated so its
/// smallest half-edge index leads (loops are cyclic; direction is
/// preserved because direction is part of the contract). Test-only
/// since the dual-walk probe retired (PK-108); the convention pins
/// use it to compare partitions.
pub(super) fn canonical_loops(loops: &[Vec<usize>]) -> BTreeSet<Vec<usize>> {
    loops
        .iter()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let min_pos = l
                .iter()
                .enumerate()
                .min_by_key(|&(_, &v)| v)
                .map(|(i, _)| i)
                .unwrap_or(0);
            let mut rotated = Vec::with_capacity(l.len());
            rotated.extend_from_slice(&l[min_pos..]);
            rotated.extend_from_slice(&l[..min_pos]);
            rotated
        })
        .collect()
}

/// Divergence between two walk partitions of one face (test-only
/// since the dual-walk probe retired — PK-108).
#[derive(Debug)]
pub(super) struct DualWalkDiff {
    /// Loops present in both partitions (canonical form).
    pub matched: usize,
    /// 3D-Newell loops the orbit partition does not contain.
    pub newell_only: Vec<Vec<usize>>,
    /// Orbits the 3D-Newell walk never emitted. Orbits consuming only
    /// half-edges the 3D-Newell walk left unwalked are benign extras (it
    /// skips some cross-registered seeds); orbits sharing half-edges
    /// with differing 3D-Newell loops are true routing divergences.
    pub orbit_only: Vec<Vec<usize>>,
}

impl DualWalkDiff {
    /// Orbit-only loops that share half-edges with Newell-only loops —
    /// the genuinely divergent routing (as opposed to benign extras).
    pub fn contested(&self) -> usize {
        let newell_hes: BTreeSet<usize> = self.newell_only.iter().flatten().copied().collect();
        self.orbit_only
            .iter()
            .filter(|o| o.iter().any(|h| newell_hes.contains(h)))
            .count()
    }
}

/// Compare two walk partitions (test-only since PK-108).
pub(super) fn compare_walks(newell: &[Vec<usize>], orbits: &[Vec<usize>]) -> DualWalkDiff {
    let l = canonical_loops(newell);
    let o = canonical_loops(orbits);
    let matched = l.intersection(&o).count();
    DualWalkDiff {
        matched,
        newell_only: l.difference(&o).cloned().collect(),
        orbit_only: o.difference(&l).cloned().collect(),
    }
}
