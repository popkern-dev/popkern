//! The partition corpus for many bodies at one place (design 7.5): k
//! bodies meeting at one vertex, bodies sharing a corner and its three
//! faces, and the same regions whatever order the bodies come in.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{BodySet, SetExpr, partition};
use crate::error::BoolError;
use crate::options::BoolOptions;

fn block(lo: [f64; 3], hi: [f64; 3]) -> BRep {
    BRep::prism(
        &[
            Point3::new(lo[0], lo[1], lo[2]),
            Point3::new(hi[0], lo[1], lo[2]),
            Point3::new(hi[0], hi[1], lo[2]),
            Point3::new(lo[0], hi[1], lo[2]),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        hi[2] - lo[2],
        Tolerance::default(),
    )
    .expect("block")
}

/// A planar body's volume: the divergence theorem over its faces' fan
/// triangles (exact on convex planar faces).
fn volume(b: &BRep) -> f64 {
    let mut v = 0.0;
    for (fh, _) in &b.faces {
        let ring = b.face_outer_loop_positions(fh);
        let Some(p0) = ring.first() else { continue };
        let o = Vec3::new(p0.x, p0.y, p0.z);
        for w in ring[1..].windows(2) {
            let p1 = Vec3::new(w[0].x, w[0].y, w[0].z);
            let p2 = Vec3::new(w[1].x, w[1].y, w[1].z);
            v += o.dot(&p1.cross(&p2)) / 6.0;
        }
    }
    v.abs()
}

/// The eight unit octant cubes around the origin, body `k` at the
/// octant whose sign bits are `k`'s.
fn octants() -> Vec<BRep> {
    (0..8)
        .map(|k: usize| {
            let lo = [0, 1, 2].map(|axis| if k & (1 << axis) != 0 { 0.0 } else { -1.0 });
            block(lo, lo.map(|c| c + 1.0))
        })
        .collect()
}

/// Eight cubes meet at one vertex, each pair of neighbours sharing a
/// face: eight regions, twelve interfaces, and their union the 2-cube
/// with six faces.
#[test]
fn eight_octant_cubes_meet_at_one_vertex() {
    let cubes = octants();
    let refs: Vec<&BRep> = cubes.iter().collect();
    let p = partition(&refs, Tolerance::default()).expect("partition");
    let mut got: Vec<BodySet> = p.regions().iter().map(|r| r.bodies).collect();
    got.sort();
    let mut want: Vec<BodySet> = (0..8).map(|k| BodySet::of(&[k])).collect();
    want.sort();
    assert_eq!(got, want);
    let union = p
        .select(&SetExpr::Any, &BoolOptions::default())
        .expect("union");
    assert!(
        union.validate_all().is_empty(),
        "{:?}",
        union.validate_all()
    );
    assert_eq!(union.faces.len(), 6);
    assert!((volume(&union) - 8.0).abs() < 1e-9);
    let cells = p.region_solids(&SetExpr::Any).expect("region solids");
    assert_eq!(cells.solids.len(), 8);
    assert_eq!(cells.interfaces.len(), 12);
    for s in &cells.solids {
        assert!((volume(&s.body) - 1.0).abs() < 1e-9, "{:?}", s.bodies);
    }
}

/// Two octants touching only at the origin: their union is two cubes
/// at one vertex, a non-manifold answer with the contact where it is.
#[test]
fn opposite_octants_touch_at_one_vertex() {
    let cubes = octants();
    let refs: Vec<&BRep> = cubes.iter().collect();
    let p = partition(&refs, Tolerance::default()).expect("partition");
    let corners = SetExpr::Body(0).or(SetExpr::Body(7));
    match p.select(&corners, &BoolOptions::default()) {
        Err(BoolError::NonManifoldResult { contacts }) => {
            assert!(!contacts.is_empty());
            assert!(
                contacts
                    .iter()
                    .all(|c| c.location.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1e-9),
                "{contacts:?}"
            );
        }
        other => panic!("expected a vertex contact, got {other:?}"),
    }
}

/// Three cubes nested at one corner share its three faces: each
/// region a shell of the next, and the shared corner one vertex.
#[test]
fn nested_cubes_share_a_corner_and_its_faces() {
    let cubes: Vec<BRep> = [1.0, 0.75, 0.5]
        .iter()
        .map(|&s| block([0.0; 3], [s; 3]))
        .collect();
    let refs: Vec<&BRep> = cubes.iter().collect();
    let p = partition(&refs, Tolerance::default()).expect("partition");
    let mut got: Vec<BodySet> = p.regions().iter().map(|r| r.bodies).collect();
    got.sort();
    let mut want = vec![
        BodySet::of(&[0]),
        BodySet::of(&[0, 1]),
        BodySet::of(&[0, 1, 2]),
    ];
    want.sort();
    assert_eq!(got, want);
    let cells = p.region_solids(&SetExpr::Any).expect("region solids");
    let vol = |m: &[usize]| {
        cells
            .solids
            .iter()
            .find(|s| s.bodies == BodySet::of(m))
            .map(|s| volume(&s.body))
            .expect("region")
    };
    assert!((vol(&[0, 1, 2]) - 0.125).abs() < 1e-9);
    assert!((vol(&[0, 1]) - (0.421875 - 0.125)).abs() < 1e-9);
    assert!((vol(&[0]) - (1.0 - 0.421875)).abs() < 1e-9);
    let union = p
        .select(&SetExpr::Any, &BoolOptions::default())
        .expect("union");
    assert_eq!(union.faces.len(), 6);
}

/// Every permutation of `bodies` cuts the same regions: the same
/// memberships (mapped back through the permutation), each with the
/// same volume, and the union with the same faces and volume.
fn permutation_invariant(bodies: &[BRep]) {
    let n = bodies.len();
    let census = |order: &[usize]| {
        let refs: Vec<&BRep> = order.iter().map(|&i| &bodies[i]).collect();
        let p = partition(&refs, Tolerance::default()).expect("partition");
        let cells = p.region_solids(&SetExpr::Any).expect("region solids");
        let mut regions: Vec<(BodySet, i64)> = cells
            .solids
            .iter()
            .map(|s| {
                let original: Vec<usize> = s.bodies.iter().map(|k| order[k]).collect();
                (
                    BodySet::of(&original),
                    (volume(&s.body) * 1e9).round() as i64,
                )
            })
            .collect();
        regions.sort();
        let union = p
            .select(&SetExpr::Any, &BoolOptions::default())
            .expect("union");
        (
            regions,
            union.faces.len(),
            (volume(&union) * 1e9).round() as i64,
        )
    };
    let identity: Vec<usize> = (0..n).collect();
    let reference = census(&identity);
    let mut order = identity.clone();
    // Heap's algorithm over every order.
    let mut c = vec![0; n];
    let mut i = 0;
    while i < n {
        if c[i] < i {
            order.swap(if i % 2 == 0 { 0 } else { c[i] }, i);
            assert_eq!(census(&order), reference, "order {order:?}");
            c[i] += 1;
            i = 0;
        } else {
            c[i] = 0;
            i += 1;
        }
    }
}

/// Three cubes on diagonal offsets, every pair and the triple
/// overlapping.
#[test]
fn diagonal_cubes_are_permutation_invariant() {
    permutation_invariant(&[
        block([0.0; 3], [1.0; 3]),
        block([0.5, 0.3, 0.2], [1.5, 1.3, 1.2]),
        block([0.25, 0.6, 0.45], [1.25, 1.6, 1.45]),
    ]);
}

/// Four blocks on one floor: coplanar bottoms between every pair, a
/// pinched shared boundary in the union's bottom skin.
#[test]
fn coplanar_blocks_are_permutation_invariant() {
    permutation_invariant(&[
        block([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        block([0.5, 0.0, 0.0], [1.5, 1.0, 0.8]),
        block([0.25, 0.5, 0.0], [1.25, 1.5, 1.2]),
        block([0.1, 0.2, 0.0], [0.9, 0.7, 0.5]),
    ]);
}
