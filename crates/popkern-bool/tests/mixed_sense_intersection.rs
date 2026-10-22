//! PK-911: a disc intersected with an inward-shelled boss (`BRep::prism`
//! primitives are inward shells) in either operand order. The disc-first
//! order once served an invalid body: the cellular merge wired the
//! mixed-sense pair's record edges one way from both sides, and the
//! graph-merge fallback served record edges whose pcurves ran against
//! their curves (recorded edge tolerance 4.0). The served body must now
//! validate and hold the true slab volume in both orders.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Volume of an all-planar body with straight edges: the divergence
/// theorem over each face's vertex fan.
fn polyhedral_volume(brep: &BRep) -> f64 {
    signed_volume(brep).abs()
}

/// [`polyhedral_volume`] signed by the body's winding: positive for an
/// outward shell, negative for an inward one.
fn signed_volume(brep: &BRep) -> f64 {
    let mut volume = 0.0;
    for face in brep.faces.values() {
        for l in std::iter::once(face.outer_loop).chain(face.inner_loops.iter().copied()) {
            let start = brep.loops[l].half_edge;
            let mut points = Vec::new();
            let mut h = start;
            loop {
                points.push(brep.vertices[brep.half_edges[h].origin].position);
                h = brep.half_edges[h].next;
                if h == start {
                    break;
                }
            }
            let o = points[0];
            for w in points[1..].windows(2) {
                let (a, b) = (w[0], w[1]);
                volume += (o.x * (a.y * b.z - a.z * b.y) - o.y * (a.x * b.z - a.z * b.x)
                    + o.z * (a.x * b.y - a.y * b.x))
                    / 6.0;
            }
        }
    }
    volume
}

fn boss(z0: f64, tol: Tolerance) -> BRep {
    let pts = [
        Point3::new(-2.0, -2.0, z0),
        Point3::new(2.0, -2.0, z0),
        Point3::new(2.0, 2.0, z0),
        Point3::new(-2.0, 2.0, z0),
    ];
    BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), 2.0, tol).expect("boss")
}

#[test]
fn a_disc_intersected_with_an_inward_boss_validates_in_either_order() {
    let tol = Tolerance::default();
    let disc = BRep::cylinder(30.0, 1.5, tol).expect("disc");
    // A boss crossing the top plane by 0.1, and one by 1.0.
    for (z0, depth) in [(1.4, 0.1), (0.5, 1.0)] {
        let boss = boss(z0, tol);
        let want = 16.0 * depth;
        for (name, a, b) in [("disc ∩ boss", &disc, &boss), ("boss ∩ disc", &boss, &disc)] {
            let r = boolean(BoolOperation::Intersection, a.clone(), b.clone(), tol, None)
                .unwrap_or_else(|e| panic!("{name} (z0 {z0}): {e}"));
            assert!(
                r.validate().is_ok(),
                "{name} (z0 {z0}): {:?}",
                r.validate().err()
            );
            let v = polyhedral_volume(&r);
            assert!(
                (v - want).abs() < 1e-9 * want.max(1.0) + 1e-9,
                "{name} (z0 {z0}): volume {v} vs {want}"
            );
        }
    }
}

/// Every operation in both operand orders on the same mixed-sense pairs:
/// each serves a valid body, a union has the same shape either way, and
/// the all-planar results hold their true volumes.
#[test]
fn every_operation_on_the_mixed_sense_pair_validates_in_either_order() {
    let tol = Tolerance::default();
    let disc = BRep::cylinder(30.0, 1.5, tol).expect("disc");
    let counts = |r: &BRep| (r.faces.len(), r.edges.len(), r.vertices.len());
    for (z0, depth) in [(1.4, 0.1), (0.5, 1.0)] {
        let boss = boss(z0, tol);
        let run = |op, a: &BRep, b: &BRep| {
            let r = boolean(op, a.clone(), b.clone(), tol, None)
                .unwrap_or_else(|e| panic!("{op:?} (z0 {z0}): {e}"));
            assert!(
                r.validate().is_ok(),
                "{op:?} (z0 {z0}): {:?}",
                r.validate().err()
            );
            r
        };
        for op in [BoolOperation::Union, BoolOperation::Intersection] {
            let (ab, ba) = (run(op, &disc, &boss), run(op, &boss, &disc));
            assert_eq!(counts(&ab), counts(&ba), "{op:?} (z0 {z0}) is order-free");
        }
        run(BoolOperation::Subtraction, &disc, &boss);
        let stub = run(BoolOperation::Subtraction, &boss, &disc);
        let want = 16.0 * (2.0 - depth);
        let v = polyhedral_volume(&stub);
        assert!(
            (v - want).abs() < 1e-9 * want,
            "boss − disc (z0 {z0}): {v} vs {want}"
        );
    }
}

/// PK-916: a mixed-sense union whose tool's far faces the scope
/// carries (an inward 16-gon boss of radius 10 on an outward 2 x 2
/// block) serves with the first operand's sense whole in either order:
/// the carried far faces reverse with the tool's flipped fragments.
/// Carried unreversed, the scoped union validated with the block's part
/// outward and the boss's inward (signed volume -606.7, true 617.9).
#[test]
fn a_scoped_mixed_sense_union_keeps_one_sense_in_either_order() {
    let tol = Tolerance::default();
    let square = [
        Point3::new(-1.0, -1.0, 0.0),
        Point3::new(1.0, -1.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(-1.0, 1.0, 0.0),
    ];
    let mut block = BRep::prism(&square, Vec3::new(0.0, 0.0, 1.0), 1.5, tol).expect("block");
    block.reverse_body_orientation();
    let boss_area = 16.0 * 0.5 * 100.0 * (std::f64::consts::TAU / 16.0).sin();
    for z0 in [1.4, 0.5] {
        let pts: Vec<Point3> = (0..16)
            .map(|i| {
                let t = std::f64::consts::TAU * i as f64 / 16.0;
                Point3::new(10.0 * t.cos(), 10.0 * t.sin(), z0)
            })
            .collect();
        let boss = BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), 2.0, tol).expect("boss");
        let want = 6.0 + 2.0 * boss_area - 4.0 * (1.5 - z0);
        for (name, a, b) in [
            ("block ∪ boss", &block, &boss),
            ("boss ∪ block", &boss, &block),
        ] {
            let r = boolean(BoolOperation::Union, a.clone(), b.clone(), tol, None)
                .unwrap_or_else(|e| panic!("{name} (z0 {z0}): {e}"));
            assert!(
                r.validate().is_ok(),
                "{name} (z0 {z0}): {:?}",
                r.validate().err()
            );
            let (v, sense) = (signed_volume(&r), signed_volume(a).signum());
            assert!(
                (v - sense * want).abs() < 1e-9 * want,
                "{name} (z0 {z0}): signed volume {v} vs {}",
                sense * want
            );
        }
    }
}
