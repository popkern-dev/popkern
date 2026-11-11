//! Contacts are named (PK-1160, naming scheme 6).
//!
//! A ball resting on a cube, united under the general policy, touches
//! its top at a point inside both faces: the glued body holds the
//! contact as a vertex loop of each face, one vertex with no edges.
//! Scheme 5 named the pieces before the glue and left that vertex
//! unnamed (and an edge split at a contact, the piece the split made);
//! a reference to the contact had nothing to hold. Scheme 6 names a
//! general result after the glue and names an edgeless contact by the
//! faces it lies on, so the name survives the ball moving across the
//! face. A tree recorded under scheme 5 replays under it, unchanged.

use popkern_bool::{
    BoolOperation, BoolOptions, NonManifoldPolicy, boolean_bodies, boolean_with_options,
};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};

const GENERAL: BoolOptions = BoolOptions {
    non_manifold: NonManifoldPolicy::General,
};

/// A ball of radius `r` turned so no pole or seam lies at its bottom,
/// its bottom at `(x, y, z)`.
fn ball(r: f64, at: (f64, f64, f64), tol: Tolerance) -> BRep {
    let turn =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0), 1.0)
            .expect("axis");
    BRep::sphere(r, tol)
        .and_then(|b| b.transform(&turn))
        .and_then(|b| b.transform(&Transform3::translation(at.0, at.1, at.2 + r)))
        .expect("ball")
}

/// `body` with every face, edge and vertex named by operation `op`, in
/// slot order, as a feature tree names a primitive.
fn named(mut body: BRep, op: u64) -> BRep {
    use popkern_core::naming::EntityRole;
    let mint = |pool: &mut popkern_core::naming::NamePool, i: usize| {
        pool.intern(
            EntityName {
                operation: OperationId(op),
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: i as u32,
                },
                ancestors: Vec::new(),
            },
            None,
        )
    };
    let faces: Vec<_> = body.faces.keys().collect();
    for (i, f) in faces.into_iter().enumerate() {
        body.faces[f].name = Some(mint(&mut body.name_pool, i));
    }
    let edges: Vec<_> = body.edges.keys().collect();
    for (i, e) in edges.into_iter().enumerate() {
        body.edges[e].name = Some(mint(&mut body.name_pool, 100 + i));
    }
    let vertices: Vec<_> = body.vertices.keys().collect();
    for (i, v) in vertices.into_iter().enumerate() {
        body.vertices[v].name = Some(mint(&mut body.name_pool, 200 + i));
    }
    body
}

fn ball_on_cube(x: f64, y: f64, scheme: u32) -> BRep {
    let tol = Tolerance::default();
    let cube = named(BRep::unit_cube(tol).expect("cube"), 1);
    let mut ctx = TrackingContext::new(OperationId(7)).with_naming_scheme(scheme);
    boolean_with_options(
        BoolOperation::Union,
        cube,
        named(ball(0.25, (x, y, 1.0), tol), 2),
        tol,
        Some(&mut ctx),
        &GENERAL,
    )
    .expect("general union")
}

/// The vertex of the body's vertex loops (the contact).
fn contact(body: &BRep) -> VertexHandle {
    let loops = body.vertex_loops();
    assert!(!loops.is_empty(), "no vertex loop");
    assert!(
        loops.iter().all(|&(_, v)| v == loops[0].1),
        "one contact vertex"
    );
    loops[0].1
}

/// A vertex's name with its ancestors spelled out, comparable across
/// bodies (pool ids are not).
fn spelled(body: &BRep, v: VertexHandle) -> Option<(RoleCategory, u32, Vec<EntityName>)> {
    let name = body.name_pool.get(body.vertices[v].name?)?;
    let ancestors = name
        .ancestors
        .iter()
        .map(|&id| body.name_pool.get(id).cloned())
        .collect::<Option<Vec<_>>>()?;
    Some((name.role.category, name.role.index, ancestors))
}

fn unnamed(body: &BRep) -> (usize, usize, usize) {
    (
        body.vertices.values().filter(|v| v.name.is_none()).count(),
        body.edges.values().filter(|e| e.name.is_none()).count(),
        body.faces.values().filter(|f| f.name.is_none()).count(),
    )
}

#[test]
fn a_point_contact_is_named_by_the_faces_it_lies_on() {
    let body = ball_on_cube(0.5, 0.5, 6);
    assert!(body.validate_all().is_empty());
    assert_eq!(unnamed(&body), (0, 0, 0), "every entity named");
    let v = contact(&body);
    let (category, index, ancestors) = spelled(&body, v).expect("named contact");
    assert_eq!((category, index), (RoleCategory::Contact, 0));
    // Its ancestors are the names of the faces holding it: the cube's
    // top and the ball's face, nothing else.
    let mut holding: Vec<EntityName> = body
        .vertex_loops()
        .into_iter()
        .map(|(l, _)| body.loops[l].face)
        .filter_map(|f| body.name_pool.get(body.faces[f].name?).cloned())
        .collect();
    holding.sort_by_key(|n| format!("{n:?}"));
    let mut got = ancestors.clone();
    got.sort_by_key(|n| format!("{n:?}"));
    assert_eq!(got, holding);
    assert_eq!(got.len(), 2);

    // The ball moved across the top: the same contact, the same name.
    for (x, y) in [(0.3, 0.6), (0.7, 0.35)] {
        let moved = ball_on_cube(x, y, 6);
        assert_eq!(
            spelled(&moved, contact(&moved)),
            Some((category, index, ancestors.clone())),
            "ball at ({x}, {y})"
        );
    }
}

#[test]
fn a_tree_recorded_under_scheme_five_names_as_it_did() {
    // Scheme 5 named the pieces before the glue: the contact vertex
    // has no name, everything else does, as before.
    let body = ball_on_cube(0.5, 0.5, 5);
    assert!(body.vertices[contact(&body)].name.is_none());
    assert_eq!(unnamed(&body), (1, 0, 0));
}

#[test]
fn a_contact_on_an_edge_names_the_edge_pieces_too() {
    // Two balls touching where one's seam passes: the seam is split at
    // the contact, and the piece the split made is named with the rest
    // (scheme 5 left it unnamed), as is the contact, a junction of the
    // seam's pieces. Through the bodies entry point, as history calls.
    let tol = Tolerance::default();
    let turn =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0), 1.0)
            .expect("axis");
    let one = || {
        BRep::sphere(1.0, tol)
            .and_then(|b| b.transform(&turn))
            .expect("ball")
    };
    let other = one()
        .transform(&Transform3::translation(2.0, 0.0, 0.0))
        .expect("moved");
    for (scheme, want) in [(6, (0, 0, 0)), (5, (1, 1, 0))] {
        let mut ctx = TrackingContext::new(OperationId(3)).with_naming_scheme(scheme);
        let bodies = boolean_bodies(
            BoolOperation::Union,
            one(),
            other.clone(),
            tol,
            Some(&mut ctx),
            &GENERAL,
        )
        .expect("general union");
        assert_eq!(bodies.len(), 1);
        assert!(bodies[0].validate_all().is_empty());
        assert_eq!(unnamed(&bodies[0]), want, "scheme {scheme}");
    }
}
