//! The 64-slot disc's last cut (PK-813): three faces active, the 252
//! other slot walls far, and the unscoped result honours every premise;
//! a scope that leaves out the faces the tool touches is refused.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{active_scope, scope_violation};
use crate::error::BoolOperation;
use crate::merge::graph::Operand;

fn slot(k: usize, tol: Tolerance) -> BRep {
    let a = std::f64::consts::TAU * k as f64 / 64.0;
    let (ca, sa, w) = (a.cos(), a.sin(), 0.6);
    let pts = [
        Point3::new(22.0 * ca + w * sa, 22.0 * sa - w * ca, -1.0),
        Point3::new(27.0 * ca + w * sa, 27.0 * sa - w * ca, -1.0),
        Point3::new(27.0 * ca - w * sa, 27.0 * sa + w * ca, -1.0),
        Point3::new(22.0 * ca - w * sa, 22.0 * sa + w * ca, -1.0),
    ];
    BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), 4.0, tol).expect("slot prism")
}

#[test]
fn the_last_slot_cut_touches_three_faces_and_carries_the_rest() {
    let tol = Tolerance::default();
    let mut body = BRep::cylinder(30.0, 1.5, tol).expect("disc");
    for k in 0..63 {
        body = crate::boolean(BoolOperation::Subtraction, body, slot(k, tol), tol, None)
            .expect("slot");
    }
    let tool = slot(63, tol);
    let hits = crate::intersect::find_intersections(&body, &tool, &tol).expect("intersections");
    let scope = active_scope(&body, &hits.aabbs_a, &hits.aabbs_b, 10.0 * tol.linear());
    assert_eq!(
        (scope.active.len(), scope.far.len(), scope.frontier.len()),
        (3, 252, 504),
        "top, bottom and outer wall active; the other slots' walls far"
    );
    // Stage 2b: the top and bottom faces' 63 earlier slot holes each lie
    // away from the tool and are carried like far faces.
    assert_eq!(scope.far_holes.values().map(Vec::len).sum::<usize>(), 126);
    let result =
        crate::boolean(BoolOperation::Subtraction, body.clone(), tool, tol, None).expect("cut");
    assert_eq!(
        scope_violation(
            BoolOperation::Subtraction,
            Operand::A,
            &body,
            &hits,
            &scope,
            &result
        ),
        None
    );
    // Negative control: a pad that empties the tool box calls the faces
    // it touches far.
    let blind = active_scope(&body, &hits.aabbs_a, &hits.aabbs_b, -3.0);
    assert!(
        scope_violation(
            BoolOperation::Subtraction,
            Operand::A,
            &body,
            &hits,
            &blind,
            &result
        )
        .is_some(),
        "a scope missing the touched faces must be refused"
    );
}

/// Stage 2 (PK-814): the scoped merge — graph from the three active
/// faces, the 252 far faces carried — assembles the same body as the
/// unscoped merge: counts, and every face's surface kind with its loop
/// vertex positions bit for bit. A scope that calls a touched face far
/// cannot stitch its frontier and declines (the caller reruns unscoped).
#[test]
fn the_scoped_merge_assembles_the_unscoped_body() {
    let tol = Tolerance::default();
    let mut body = BRep::cylinder(30.0, 1.5, tol).expect("disc");
    for k in 0..63 {
        body = crate::boolean(BoolOperation::Subtraction, body, slot(k, tol), tol, None)
            .expect("slot");
    }
    let tool = slot(63, tol);
    let hits = crate::intersect::find_intersections(&body, &tool, &tol).expect("intersections");
    let scope = active_scope(&body, &hits.aabbs_a, &hits.aabbs_b, 10.0 * tol.linear());
    let merge_by = |s: Option<&super::ActiveScope>, source| {
        let s = s.map(|a| super::OperandScopes {
            a: Some(a.clone()),
            b: None,
        });
        crate::cellular::boolean::cellular_graph_merge(
            &body,
            &tool,
            &hits,
            BoolOperation::Subtraction,
            &tol,
            s.as_ref(),
            source,
        )
    };
    let merge = |s| merge_by(s, crate::complex::serve::SelectionSource::Cellular);
    let (unscoped, _) = merge(None).expect("unscoped merge");
    let (scoped, carried) = merge(Some(&scope)).expect("scoped merge");
    let counts = |b: &BRep| (b.faces.len(), b.edges.len(), b.vertices.len());
    assert_eq!(counts(&scoped), counts(&unscoped));
    assert!(scoped.validate().is_ok());
    let keys = |b: &BRep| {
        let mut k: Vec<_> = b.faces.keys().map(|f| super::face_key(b, f)).collect();
        k.sort();
        k
    };
    assert_eq!(keys(&scoped), keys(&unscoped), "same faces, bit for bit");
    assert!(
        carried.edges.len() >= scope.frontier.len(),
        "the far walls' edges are carried: {}",
        carried.edges.len()
    );
    // Negative control: the blind scope's far faces include the ones the
    // tool cuts; their frontier edges never come through unsplit.
    let blind = active_scope(&body, &hits.aabbs_a, &hits.aabbs_b, -3.0);
    let declined = merge(Some(&blind))
        .map(|(b, _)| b.validate().is_err() || counts(&b) != counts(&unscoped))
        .unwrap_or(true);
    assert!(declined, "a blind scope must not assemble the cut body");
    // The cell complex decides the scoped merge alike (PK-1037): the
    // frontier edges it lifts bound open shells, and the far material
    // is carried. The blind scope leaves the top face without its rim:
    // its slot ring is the face's only boundary, both ways round, and
    // the complex declines instead of guessing.
    let complex = crate::complex::serve::SelectionSource::Complex;
    let (by_complex, _) = merge_by(Some(&scope), complex).expect("the complex serves the scope");
    assert_eq!(
        keys(&by_complex),
        keys(&unscoped),
        "same faces, bit for bit"
    );
    assert!(
        merge_by(Some(&blind), complex).is_err(),
        "the complex declines the blind scope"
    );
}

/// A far hole's rim bounds the scope even when the face inside the
/// hole is active: the rim's half-edge on the holed face never enters
/// the graph, the inner face's does (measured: a pipped die's dimple
/// rims, whose dimple faces stayed active, read as free edges in the
/// cell complex).
#[test]
fn a_far_hole_rim_bounds_the_scope() {
    let tol = Tolerance::default();
    let block = BRep::prism(
        &[
            Point3::new(-1.0, -1.0, 0.0),
            Point3::new(1.0, -1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(-1.0, 1.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        tol,
    )
    .expect("block");
    let pocket = BRep::cylinder(0.3, 1.0, tol)
        .expect("pocket")
        .transform(&popkern_core::geometry::transform::Transform3::translation(
            0.0, 0.0, 0.5,
        ))
        .expect("placed");
    let body =
        crate::boolean(BoolOperation::Subtraction, block, pocket, tol, None).expect("pocket cut");
    let (top, hole) = body
        .faces
        .iter()
        .find_map(|(f, face)| face.inner_loops.first().map(|&l| (f, l)))
        .expect("the top face carries the pocket's rim as a hole");
    let mut scope = super::ActiveScope {
        active: body.faces.keys().collect(),
        ..Default::default()
    };
    let rim: Vec<_> = body
        .edges
        .keys()
        .filter(|&e| {
            let he = body.edges[e].half_edge;
            [he, body.half_edges[he].mate().unwrap_or_default()]
                .iter()
                .any(|&h| body.half_edges[h].loop_ == hole)
        })
        .collect();
    assert!(!rim.is_empty());
    // Every face active, no far hole: nothing bounds the scope.
    assert!(body.edges.keys().all(|e| !scope.bounds(&body, e)));
    scope.far_holes.insert(top, vec![hole]);
    let bounding: Vec<_> = body
        .edges
        .keys()
        .filter(|&e| scope.bounds(&body, e))
        .collect();
    assert_eq!(bounding, rim);
}

/// The 64-slot disc as the TOOL of a union with a small inward boss
/// (PK-816 mirrored): the tool's scope carries its top face's 63 slot
/// holes (its bottom face is far, carried whole) like the target's would, reversed with its far faces
/// (the boss winds opposite to the disc), and both classifiers assemble
/// the unscoped body bit for bit. Without the holes the tool's active
/// top face brought every slot rim into the graph: the cell complex
/// lifted 252 rim edges the far region never needed (boss ∪ disc 1.3-1.46x
/// disc ∪ boss on the perf gate).
#[test]
fn the_disc_as_tool_carries_its_far_holes() {
    let tol = Tolerance::default();
    let mut disc = BRep::cylinder(30.0, 1.5, tol).expect("disc");
    for k in 0..63 {
        disc = crate::boolean(BoolOperation::Subtraction, disc, slot(k, tol), tol, None)
            .expect("slot");
    }
    let boss = BRep::prism(
        &[
            Point3::new(-2.0, -2.0, 1.4),
            Point3::new(2.0, -2.0, 1.4),
            Point3::new(2.0, 2.0, 1.4),
            Point3::new(-2.0, 2.0, 1.4),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        2.0,
        tol,
    )
    .expect("boss");
    let op = BoolOperation::Union;
    let hits = crate::intersect::find_intersections(&boss, &disc, &tol).expect("intersections");
    let scopes = super::serving_scopes(&boss, &disc, &hits, op, &tol);
    let carried = scopes
        .b
        .as_ref()
        .map_or(0, |s| s.far_holes.values().map(Vec::len).sum());
    assert_eq!(carried, 63, "the tool's slot holes are carried");
    let merge = |s: Option<&super::OperandScopes>, source| {
        crate::cellular::boolean::cellular_graph_merge(&boss, &disc, &hits, op, &tol, s, source)
    };
    let keys = |b: &BRep| {
        let mut k: Vec<_> = b.faces.keys().map(|f| super::face_key(b, f)).collect();
        k.sort();
        k
    };
    use crate::complex::serve::SelectionSource;
    let (unscoped, _) = merge(None, SelectionSource::Cellular).expect("unscoped merge");
    for source in [SelectionSource::Cellular, SelectionSource::Complex] {
        let (scoped, _) = merge(Some(&scopes), source).expect("scoped merge");
        assert!(scoped.validate().is_ok(), "{source:?}");
        assert_eq!(keys(&scoped), keys(&unscoped), "{source:?}: same faces");
    }
}
