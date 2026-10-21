//! Stage 4 (PK-816): the tool's scope, multi-solid operands. A small
//! target against a slotted disc scopes the DISC as the tool: a union
//! carries its far faces, a subtraction and an intersection drop them,
//! and in every case the scoped merge assembles the unscoped body.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{OperandScopes, active_scope, scope_violation, serving_scopes};
use crate::error::BoolOperation;
use crate::merge::graph::Operand;

fn prism(lo: Point3, hi: Point3, tol: Tolerance) -> BRep {
    let pts = [
        Point3::new(lo.x, lo.y, lo.z),
        Point3::new(hi.x, lo.y, lo.z),
        Point3::new(hi.x, hi.y, lo.z),
        Point3::new(lo.x, hi.y, lo.z),
    ];
    BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), hi.z - lo.z, tol).expect("prism")
}

/// A disc with twelve radial slots between 22 and 27 from the axis,
/// around its first quarter turn: 48 slot walls, all far from the
/// origin and from the rim at +x.
fn slotted(tol: Tolerance) -> BRep {
    let mut body = BRep::cylinder(30.0, 1.5, tol).expect("disc");
    for k in 1..=12 {
        let a = std::f64::consts::TAU * k as f64 / 64.0;
        let (ca, sa, w) = (a.cos(), a.sin(), 0.6);
        let pts = [
            Point3::new(22.0 * ca + w * sa, 22.0 * sa - w * ca, -1.0),
            Point3::new(27.0 * ca + w * sa, 27.0 * sa - w * ca, -1.0),
            Point3::new(27.0 * ca - w * sa, 27.0 * sa + w * ca, -1.0),
            Point3::new(22.0 * ca - w * sa, 22.0 * sa + w * ca, -1.0),
        ];
        let slot = BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), 4.0, tol).expect("slot");
        body = crate::boolean(BoolOperation::Subtraction, body, slot, tol, None).expect("slot");
    }
    body
}

fn keys(b: &BRep) -> Vec<(String, Vec<[u64; 3]>)> {
    let mut k: Vec<_> = b.faces.keys().map(|f| super::face_key(b, f)).collect();
    k.sort();
    k
}

/// Scoped against unscoped merge of `a op b`; returns the scopes and the
/// scoped run's carried structure.
fn scoped_matches_unscoped(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
) -> (OperandScopes, crate::serve_region::Carried, BRep) {
    let hits = crate::intersect::find_intersections(a, b, &tol).expect("intersections");
    let scopes = serving_scopes(a, b, &hits, op, &tol);

    let merge = |s: Option<&OperandScopes>| {
        crate::cellular::boolean::cellular_graph_merge(a, b, &hits, op, &tol, s)
    };
    let (unscoped, _) = merge(None).expect("unscoped merge");
    let (scoped, carried) = merge(Some(&scopes)).expect("scoped merge");
    let counts = |b: &BRep| {
        (
            b.faces.len(),
            b.edges.len(),
            b.vertices.len(),
            b.shells.len(),
        )
    };
    assert_eq!(counts(&scoped), counts(&unscoped), "{op:?}");
    assert!(scoped.validate().is_ok(), "{op:?}");
    assert_eq!(
        keys(&scoped),
        keys(&unscoped),
        "{op:?}: same faces, bit for bit"
    );
    (scopes, carried, scoped)
}

#[test]
fn a_small_target_scopes_the_tool_carrying_or_dropping_its_far_faces() {
    // Compare mode serves unscoped by design (it checks the premises).
    if super::compare_enabled() {
        return;
    }
    let tol = Tolerance::default();
    let body = slotted(tol);
    // On the top face, away from every slot: the tool's top face and
    // outer wall (whose box spans the disc) are active, its 48 slot
    // walls and bottom face are not.
    let boss = prism(
        Point3::new(-2.0, -2.0, 1.4),
        Point3::new(2.0, 2.0, 3.4),
        tol,
    );
    for op in [BoolOperation::Union, BoolOperation::Subtraction] {
        let (scopes, carried, result) = scoped_matches_unscoped(op, &boss, &body, tol);
        assert!(scopes.a.is_none(), "the boss has one far face");
        let far = scopes.b.as_ref().map_or(0, |s| s.far.len());
        assert_eq!(far, 49, "{op:?}: the slot walls and bottom are far");
        // A union keeps the disc's far faces and their structure stands
        // as in the disc: the serve tail leaves them alone. A
        // subtraction keeps none of the tool's outside material.
        let region = crate::serve_region::ServeRegion::of(&result, &carried);
        let untouched = result
            .faces
            .keys()
            .filter(|&f| carried.b_faces.contains_key(f) && !region.admits_face(f))
            .count();
        match op {
            BoolOperation::Union => assert!(untouched >= far - 2, "{untouched} of {far}"),
            _ => assert!(carried.b_faces.is_empty() && carried.b_edges.is_empty()),
        }
    }
    // Through the rim at +x, clear of every slot: an intersection drops
    // the far faces of both operands in either order.
    let rim = prism(
        Point3::new(29.0, -1.0, -2.0),
        Point3::new(31.0, 1.0, 2.0),
        tol,
    );
    let (scopes, ..) = scoped_matches_unscoped(BoolOperation::Intersection, &body, &rim, tol);
    assert!(
        scopes
            .a
            .as_ref()
            .is_some_and(|s| s.far.len() >= 45 && s.far_holes.is_empty())
    );
    let (scopes, carried, _) =
        scoped_matches_unscoped(BoolOperation::Intersection, &rim, &body, tol);
    assert!(scopes.b.is_some() && carried.b_faces.is_empty());
}

#[test]
fn the_tools_premises_hold_and_a_blind_tool_scope_is_refused() {
    let tol = Tolerance::default();
    let body = slotted(tol);
    let boss = prism(
        Point3::new(-2.0, -2.0, 1.4),
        Point3::new(2.0, 2.0, 3.4),
        tol,
    );
    let hits = crate::intersect::find_intersections(&boss, &body, &tol).expect("intersections");
    let pad = 10.0 * tol.linear();
    let scope = active_scope(&body, &hits.aabbs_b, &hits.aabbs_a, pad);
    // Negative control: a pad that empties the boss's box calls the
    // disc's top face far, which the boss changes.
    let blind = active_scope(&body, &hits.aabbs_b, &hits.aabbs_a, -3.0);
    for op in [BoolOperation::Union, BoolOperation::Subtraction] {
        let result = crate::boolean(op, boss.clone(), body.clone(), tol, None).expect("boolean");
        assert_eq!(
            scope_violation(op, Operand::B, &body, &hits, &scope, &result),
            None,
            "{op:?}"
        );
        assert!(
            scope_violation(op, Operand::B, &body, &hits, &blind, &result).is_some(),
            "{op:?}: a tool scope missing the touched face must be refused"
        );
    }
}

/// A body of two solids cut at one of them: the other solid's faces are
/// far, carried into the one assembly shell, and the shells come apart
/// again as two solids, as unscoped.
#[test]
fn a_two_solid_body_cut_at_one_solid_keeps_the_other_whole() {
    if super::compare_enabled() {
        return;
    }
    let tol = Tolerance::default();
    let body = slotted(tol);
    let far_block = prism(
        Point3::new(100.0, 0.0, 0.0),
        Point3::new(104.0, 4.0, 4.0),
        tol,
    );
    let two = crate::boolean(BoolOperation::Union, body.clone(), far_block, tol, None)
        .expect("two solids");
    assert_eq!((two.solids.len(), two.shells.len()), (2, 2));
    let pocket = prism(
        Point3::new(101.0, 1.0, -1.0),
        Point3::new(103.0, 3.0, 5.0),
        tol,
    );
    let (scopes, _, result) =
        scoped_matches_unscoped(BoolOperation::Subtraction, &two, &pocket, tol);
    let far = scopes.a.as_ref().map_or(0, |s| s.far.len());
    assert!(far >= body.faces.len(), "the whole disc is far ({far})");
    assert_eq!((result.solids.len(), result.shells.len()), (2, 2));
}
