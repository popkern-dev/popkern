//! Recorded contacts survive later operators (PK-986). Two cubes
//! touching along one edge are served as pieces meeting at recorded
//! contacts (`NonManifoldPolicy::SplitToManifold`); later booleans
//! under the default policy carry that deliberate touch (recorded
//! again, `validate_all` clean, volumes exact) and still decline a
//! contact no operand recorded. A fillet and a rigid motion keep the
//! records as they are. The kernel mesh never welds across a contact:
//! each face meshes its own vertices, so the two copies of a contact
//! edge stay two (the welded census counts the edge used four times,
//! not open).
use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_bool::{
    BoolError, BoolOperation, BoolOptions, NonManifoldPolicy, boolean, boolean_with_options,
};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn cube(x: f64, y: f64, z: f64, s: f64) -> BRep {
    let tol = Tolerance::default();
    BRep::unit_cube(tol)
        .and_then(|c| c.transform(&Transform3::scale(s, s, s)))
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube")
}

const SPLIT: BoolOptions = BoolOptions {
    non_manifold: NonManifoldPolicy::SplitToManifold,
};

/// Unit cubes at the origin and at (1, 1, 0): one edge contact on the
/// line x = y = 1.
fn touching() -> BRep {
    let tol = Tolerance::default();
    let body = boolean_with_options(
        BoolOperation::Union,
        cube(0.0, 0.0, 0.0, 1.0),
        cube(1.0, 1.0, 0.0, 1.0),
        tol,
        None,
        &SPLIT,
    )
    .expect("touching cubes serve as pieces");
    assert!(!body.topology_contacts.is_empty());
    body
}

/// Served with contacts, `validate_all` clean, the expected volume.
fn served(what: &str, body: &BRep, volume: f64) {
    let tol = Tolerance::default();
    assert!(!body.topology_contacts.is_empty(), "{what}: contacts kept");
    let issues = body.validate_all();
    assert!(issues.is_empty(), "{what}: {issues:?}");
    let v = popkern_measure::body_volume(body, &tol).expect("volume");
    assert!((v - volume).abs() < 1e-6, "{what}: volume {v} vs {volume}");
}

#[test]
fn booleans_carry_the_contacts_their_operands_recorded() {
    let tol = Tolerance::default();
    let cases = [
        (
            "union with a far cube",
            BoolOperation::Union,
            cube(5.0, 0.0, 0.0, 1.0),
            3.0,
        ),
        // x -0.5..0.1 over piece A: 0.216 less the 0.036 inside A.
        (
            "union overlapping a piece",
            BoolOperation::Union,
            cube(-0.5, 0.2, 0.2, 0.6),
            2.18,
        ),
        (
            "a pocket in a piece",
            BoolOperation::Subtraction,
            cube(0.2, 0.2, 0.5, 0.3),
            1.973,
        ),
        // 0.2 x 0.2 x 0.4 out of each piece, the contact cut in two.
        (
            "a cut across the contact",
            BoolOperation::Subtraction,
            cube(0.8, 0.8, 0.4, 0.4),
            1.968,
        ),
        (
            "inside a big box",
            BoolOperation::Intersection,
            cube(-1.0, -1.0, -1.0, 4.0),
            2.0,
        ),
    ];
    for (what, op, tool, volume) in cases {
        let body =
            boolean(op, touching(), tool, tol, None).unwrap_or_else(|e| panic!("{what}: {e}"));
        served(what, &body, volume);
    }
    // Both operands carry contacts.
    let other = touching()
        .transform(&Transform3::translation(0.0, 0.0, 5.0))
        .expect("moved");
    let body = boolean(BoolOperation::Union, touching(), other, tol, None).expect("two");
    served("two touching pairs", &body, 4.0);
}

#[test]
fn a_contact_no_operand_recorded_still_declines() {
    let tol = Tolerance::default();
    // A third cube touching piece A along its edge x = y = 0: new.
    let third = || cube(-1.0, -1.0, 0.0, 1.0);
    match boolean(BoolOperation::Union, touching(), third(), tol, None) {
        Err(BoolError::NonManifoldResult { contacts }) => {
            assert!(!contacts.is_empty());
            assert!(
                contacts
                    .iter()
                    .all(|c| c.location.x.abs() < 1e-6 && c.location.y.abs() < 1e-6),
                "only the new contact is reported: {contacts:?}"
            );
        }
        other => panic!("a new edge contact must decline: {other:?}"),
    }
    // Asked for, it serves: three pieces, every contact recorded.
    let body = boolean_with_options(BoolOperation::Union, touching(), third(), tol, None, &SPLIT)
        .expect("split serves");
    served("three pieces", &body, 3.0);
    assert_eq!(body.solid_regions().len(), 3);
}

#[test]
fn a_fillet_and_a_motion_keep_the_records() {
    let tol = Tolerance::default();
    let base = touching();
    // The z edge of piece A at x = y = 0, away from the contact.
    let edge = base
        .edges
        .iter()
        .find(|(_, e)| {
            let he = &base.half_edges[e.half_edge];
            let a = base.vertices[he.origin].position;
            let b = base.vertices[base.half_edges[he.mate().unwrap_or_default()].origin].position;
            [a, b].iter().all(|p| p.x.abs() < 1e-9 && p.y.abs() < 1e-9)
        })
        .map(|(h, _)| h)
        .expect("edge");
    let spec = FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    };
    let filleted = fillet_edges(&base, &[spec], &FilletOptions::default(), tol, None)
        .expect("fillet")
        .brep;
    assert_eq!(filleted.topology_contacts, base.topology_contacts);
    served(
        "filleted",
        &filleted,
        2.0 - 0.01 * (1.0 - std::f64::consts::FRAC_PI_4),
    );
    let moved = base
        .transform(&Transform3::translation(1.0, 2.0, 3.0))
        .expect("moved");
    assert_eq!(moved.topology_contacts.len(), base.topology_contacts.len());
    served("moved", &moved, 2.0);
}

#[test]
fn the_mesh_keeps_both_sides_of_a_contact() {
    let body = touching();
    let mesh = popkern_mesh::tessellate(&body, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    assert_eq!(popkern_mesh::welded_open_edges(&mesh), 0);
    // No mesh vertex is shared between faces of different pieces.
    let solid_of = |f| body.face_back_shell(f).and_then(|s| body.shell_region(s));
    let mut owner = std::collections::HashMap::new();
    for (t, tri) in mesh.triangles.iter().enumerate() {
        let solid = solid_of(mesh.face_map[t]);
        for &v in tri {
            assert_eq!(*owner.entry(v).or_insert(solid), solid, "vertex {v} welded");
        }
    }
}

#[test]
fn partition_operations_carry_recorded_contacts() {
    let tol = Tolerance::default();
    let base = touching();
    let far = cube(5.0, 0.0, 0.0, 1.0);
    let any = popkern_bool::nary_boolean(
        &[&base, &far],
        &popkern_bool::SetExpr::Any,
        tol,
        &BoolOptions::default(),
    );
    match &any {
        Ok(body) => served("nary any", body, 3.0),
        Err(e) => panic!("nary any: {e}"),
    }
    // A tool above z = 0.5 cuts each piece in two: four pieces, the
    // pieces of different cubes still meeting along the contact line.
    let cut = popkern_bool::split(&base, &[&cube(-1.0, -1.0, 0.5, 4.0)], tol).expect("split");
    assert_eq!(cut.pieces.len(), 4);
    for piece in &cut.pieces {
        let v = popkern_measure::body_volume(&piece.body, &tol).expect("volume");
        assert!((v - 0.5).abs() < 1e-9, "{v}");
    }
}
