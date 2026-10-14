//! Naming scheme 4 (PK-803): a boolean's carried entities — the
//! target's faces, edges and vertices left unchanged and clear of the
//! tool — keep the target's names; entities the tool changes still mint
//! names whose ancestry reaches what they came from.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, OperationId};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn slot(k: usize, n: usize, tol: Tolerance) -> BRep {
    let a = std::f64::consts::TAU * k as f64 / n as f64;
    let (ca, sa, w) = (a.cos(), a.sin(), 0.6);
    let pts = [
        Point3::new(22.0 * ca + w * sa, 22.0 * sa - w * ca, -1.0),
        Point3::new(27.0 * ca + w * sa, 27.0 * sa - w * ca, -1.0),
        Point3::new(27.0 * ca - w * sa, 27.0 * sa + w * ca, -1.0),
        Point3::new(22.0 * ca - w * sa, 22.0 * sa + w * ca, -1.0),
    ];
    BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), 4.0, tol).expect("slot prism")
}

/// The disc after each of `n` slot cuts under naming `scheme`.
fn slotted(n: usize, scheme: u32) -> Vec<BRep> {
    let tol = Tolerance::default();
    let mut bodies = vec![BRep::cylinder(30.0, 1.5, tol).expect("disc")];
    for k in 0..n {
        let mut ctx = TrackingContext::new(OperationId(100 + k as u64)).with_naming_scheme(scheme);
        let target = bodies.last().cloned().expect("a body");
        let cut = boolean(
            BoolOperation::Subtraction,
            target,
            slot(k, n, tol),
            tol,
            Some(&mut ctx),
        )
        .expect("slot cut serves");
        bodies.push(cut);
    }
    bodies
}

/// The name of the face whose outer-loop vertex centroid lies nearest
/// `probe` — located by geometry, not handle.
fn face_name_near(body: &BRep, probe: Point3) -> Option<EntityName> {
    body.faces
        .iter()
        .filter_map(|(fh, f)| {
            let ps = body.face_outer_loop_positions(fh);
            let n = ps.len().max(1) as f64;
            let c = ps.iter().fold(Point3::new(0.0, 0.0, 0.0), |acc, p| {
                Point3::new(acc.x + p.x / n, acc.y + p.y / n, acc.z + p.z / n)
            });
            Some((c.distance_to(&probe), f.name?))
        })
        .min_by(|a, b| a.0.total_cmp(&b.0))
        .and_then(|(_, id)| body.name_pool.get(id).cloned())
}

/// The first slot's outer radial wall (vertex centroid at x = 27).
fn first_wall_corner() -> Point3 {
    Point3::new(27.0, 0.0, 0.0)
}

#[test]
fn carried_faces_keep_their_names_through_later_cuts() {
    let bodies = slotted(6, 4);
    let named_at_cut = face_name_near(&bodies[1], first_wall_corner()).expect("wall named");
    for (k, body) in bodies.iter().enumerate().skip(2) {
        assert_eq!(
            face_name_near(body, first_wall_corner()).as_ref(),
            Some(&named_at_cut),
            "the first slot's wall keeps its name through cut {k}"
        );
    }
}

#[test]
fn scheme_three_renames_the_same_face_every_cut() {
    // Control: under scheme 3 the same wall is minted a new Modified name
    // by every later boolean — the behaviour scheme 4 retires.
    let bodies = slotted(3, 3);
    let first = face_name_near(&bodies[1], first_wall_corner()).expect("wall named");
    let later = face_name_near(&bodies[3], first_wall_corner()).expect("wall named");
    assert_ne!(first, later);
}

#[test]
fn touched_faces_still_descend_from_their_inputs() {
    // The disc's top cap is cut by every slot: it is never carried, and
    // its latest name still descends from the cap the first cut saw.
    let bodies = slotted(4, 4);
    let top = Point3::new(30.0, 0.0, 0.75);
    let first = face_name_near(&bodies[1], top).expect("cap named");
    let last_body = bodies.last().expect("a body");
    let last_id = last_body
        .faces
        .values()
        .filter_map(|f| f.name)
        .find(|id| last_body.name_pool.get(*id) == face_name_near(last_body, top).as_ref())
        .expect("cap id");
    assert_ne!(last_body.name_pool.get(last_id), Some(&first));
    assert!(
        last_body.name_pool.descends(last_id, &first),
        "the cut cap's name reaches the cap it was cut from"
    );
}

#[test]
fn carried_entities_record_no_new_names() {
    // The pool grows by what the cut changed, not by the body: a later
    // cut adds far fewer names than the body has entities.
    let bodies = slotted(8, 4);
    let (before, after) = (&bodies[7], &bodies[8]);
    let grown = after.name_pool.names().count() - before.name_pool.names().count();
    let entities = after.faces.len() + after.edges.len() + after.vertices.len();
    assert!(
        grown * 4 < entities,
        "a slot cut minted {grown} names on a body of {entities} entities"
    );
}
