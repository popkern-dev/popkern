//! PK-913: a mixed-sense pair (an outward disc, an inward prism boss)
//! assembles on the cellular path itself, with no rescue retry, for
//! every operation in both orders and both PK-911 crossing depths,
//! scoped as the serve scopes it and unscoped; a same-sense pair is
//! left alone.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BoolOperation;

fn boss(z0: f64, tol: Tolerance) -> BRep {
    let pts = [
        Point3::new(-2.0, -2.0, z0),
        Point3::new(2.0, -2.0, z0),
        Point3::new(2.0, 2.0, z0),
        Point3::new(-2.0, 2.0, z0),
    ];
    BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), 2.0, tol).expect("boss")
}

fn merge(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    scoped: bool,
) -> Result<BRep, String> {
    let hits = crate::intersect::find_intersections(a, b, &tol).map_err(|e| format!("{e:?}"))?;
    let scopes = crate::scope::serving_scopes(a, b, &hits, op, &tol);
    crate::cellular::boolean::cellular_graph_merge(a, b, &hits, op, &tol, scoped.then_some(&scopes))
        .map(|(body, _)| body)
        .map_err(|d| format!("{d:?}").chars().take(300).collect())
}

#[test]
fn a_mixed_sense_pair_assembles_on_the_cellular_path() {
    let tol = Tolerance::default();
    let disc = BRep::cylinder(30.0, 1.5, tol).expect("disc");
    let mut failures = Vec::new();
    for z0 in [1.4, 0.5] {
        let boss = boss(z0, tol);
        for op in [
            BoolOperation::Union,
            BoolOperation::Intersection,
            BoolOperation::Subtraction,
        ] {
            for (order, a, b) in [("disc, boss", &disc, &boss), ("boss, disc", &boss, &disc)] {
                for scoped in [false, true] {
                    let label = format!("{op:?} ({order}) z0 {z0} scoped {scoped}");
                    match merge(op, a, b, tol, scoped) {
                        Ok(body) if body.validate().is_ok() => {}
                        Ok(body) => failures.push(format!("{label}: {:?}", body.validate().err())),
                        Err(e) => failures.push(format!("{label}: {e}")),
                    }
                }
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// An inward 16-gon prism, radius 10, from `z0` up 2: its sides and
/// top lie clear of the 2 x 2 block at the origin, so a scoped union
/// carries them as far faces (PK-816).
fn wide_boss(z0: f64, tol: Tolerance) -> BRep {
    let pts: Vec<Point3> = (0..16)
        .map(|i| {
            let t = std::f64::consts::TAU * i as f64 / 16.0;
            Point3::new(10.0 * t.cos(), 10.0 * t.sin(), z0)
        })
        .collect();
    BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), 2.0, tol).expect("boss")
}

/// Volume of an all-planar body: the divergence theorem over each
/// face's vertex fans, signed by winding.
fn planar_volume(brep: &BRep) -> f64 {
    let mut volume = 0.0;
    for face in brep.faces.values() {
        for l in std::iter::once(face.outer_loop).chain(face.inner_loops.iter().copied()) {
            let start = brep.loops[l].half_edge;
            let mut p = Vec::new();
            let mut h = start;
            loop {
                p.push(brep.vertices[brep.half_edges[h].origin].position);
                h = brep.half_edges[h].next;
                if h == start {
                    break;
                }
            }
            for w in p[1..].windows(2) {
                let (o, a, b) = (p[0], w[0], w[1]);
                volume += (o.x * (a.y * b.z - a.z * b.y) - o.y * (a.x * b.z - a.z * b.x)
                    + o.z * (a.x * b.y - a.y * b.x))
                    / 6.0;
            }
        }
    }
    volume
}

/// PK-916: with B's far faces carried (a scoped union), the flip
/// reverses them in the same pass, so a mixed-sense pair serves on
/// the cellular path with A's sense whole and the unscoped volume.
#[test]
fn a_mixed_sense_scoped_union_reverses_the_carried_far_faces() {
    let tol = Tolerance::default();
    let square = [
        Point3::new(-1.0, -1.0, 0.0),
        Point3::new(1.0, -1.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(-1.0, 1.0, 0.0),
    ];
    let mut block = BRep::prism(&square, Vec3::new(0.0, 0.0, 1.0), 1.5, tol).expect("block");
    block.reverse_body_orientation();
    let op = BoolOperation::Union;
    for z0 in [1.4, 0.5] {
        let boss = wide_boss(z0, tol);
        let hits = crate::intersect::find_intersections(&block, &boss, &tol).expect("hits");
        let scopes = crate::scope::serving_scopes(&block, &boss, &hits, op, &tol);
        let far_b = scopes.carried(op).1.map(|s| s.far.len());
        assert!(
            far_b.is_some_and(|n| n >= 8),
            "z0 {z0}: B far faces {far_b:?}"
        );
        let scoped =
            merge(op, &block, &boss, tol, true).unwrap_or_else(|e| panic!("z0 {z0} scoped: {e}"));
        let unscoped = merge(op, &block, &boss, tol, false)
            .unwrap_or_else(|e| panic!("z0 {z0} unscoped: {e}"));
        assert!(
            scoped.validate().is_ok(),
            "z0 {z0}: {:?}",
            scoped.validate().err()
        );
        // The block's sense (outward, positive volume) is kept whole.
        let boss_area = 16.0 * 0.5 * 100.0 * (std::f64::consts::TAU / 16.0).sin();
        let want = 6.0 + 2.0 * boss_area - 4.0 * (1.5 - z0);
        let (vs, vu) = (planar_volume(&scoped), planar_volume(&unscoped));
        assert!(
            (vu - want).abs() <= 1e-9 * want,
            "z0 {z0}: unscoped {vu} vs {want}"
        );
        assert!(
            (vs - want).abs() <= 1e-9 * want,
            "z0 {z0}: scoped {vs} vs {want}"
        );
    }
}
