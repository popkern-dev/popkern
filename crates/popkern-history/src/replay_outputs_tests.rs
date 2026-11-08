use super::*;
use popkern_core::geometry::transform::Transform3;
use popkern_core::naming::{EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;

fn name(index: u32) -> EntityName {
    EntityName {
        operation: OperationId(7),
        role: EntityRole {
            category: RoleCategory::Generated,
            index,
        },
        ancestors: Vec::new(),
    }
}

/// A unit cube at `x` whose faces are named `first..first + 6`, every
/// piece's pool interning the same names in the same order (as pieces
/// taken from one named body share its pool).
fn piece(x: f64, first: u32) -> BRep {
    let mut cube = BRep::unit_cube(Tolerance::default())
        .expect("cube")
        .transform(&Transform3::translation(x, 0.0, 0.0))
        .expect("placed");
    let ids: Vec<NameId> = (0..18)
        .map(|i| cube.name_pool.intern(name(i), None))
        .collect();
    let faces: Vec<_> = cube.faces.keys().collect();
    for (k, f) in faces.into_iter().enumerate() {
        cube.faces[f].name = Some(ids[first as usize + k]);
    }
    cube
}

fn x_of(b: &BRep) -> f64 {
    b.vertices
        .values()
        .map(|v| v.position.x)
        .fold(f64::MAX, f64::min)
}

fn counter(start: u32) -> impl FnMut() -> BodyId {
    let mut next = start;
    move || {
        next += 1;
        BodyId(next - 1)
    }
}

/// First serve: the own body takes the first piece, the others new
/// bodies, each keyed by a name only its piece carries.
#[test]
fn a_first_serve_allocates_and_keys_every_piece() {
    let pieces = vec![piece(0.0, 0), piece(2.0, 6), piece(4.0, 12)];
    let out = assign(BodyId(0), &[], pieces, &mut counter(5));
    let bodies: Vec<_> = out.served.iter().map(|(b, p)| (b.0, x_of(p))).collect();
    assert_eq!(bodies, vec![(0, 0.0), (5, 2.0), (6, 4.0)]);
    assert_eq!(out.outputs.len(), 3);
    assert_eq!(out.outputs[1].key, Some(name(6)));
    assert!(out.absent.is_empty());
}

/// One piece never recorded several: nothing is recorded.
#[test]
fn a_single_piece_records_nothing() {
    let out = assign(BodyId(0), &[], vec![piece(0.0, 0)], &mut counter(5));
    assert_eq!(out.served.len(), 1);
    assert!(out.outputs.is_empty());
}

/// Reordered pieces keep their bodies by key; a gone piece leaves its
/// body absent with its key, and comes back to it.
#[test]
fn pieces_keep_their_bodies_by_key_through_reorder_and_loss() {
    let first = assign(
        BodyId(0),
        &[],
        vec![piece(0.0, 0), piece(2.0, 6), piece(4.0, 12)],
        &mut counter(5),
    );
    let reordered = assign(
        BodyId(0),
        &first.outputs,
        vec![piece(4.0, 12), piece(0.0, 0), piece(2.0, 6)],
        &mut counter(9),
    );
    let bodies: Vec<_> = reordered
        .served
        .iter()
        .map(|(b, p)| (b.0, x_of(p)))
        .collect();
    assert_eq!(bodies, vec![(0, 0.0), (5, 2.0), (6, 4.0)]);
    let lost = assign(
        BodyId(0),
        &first.outputs,
        vec![piece(4.0, 12), piece(0.0, 0)],
        &mut counter(9),
    );
    assert_eq!(lost.absent, vec![BodyId(5)]);
    assert_eq!(lost.outputs, first.outputs, "the gone body keeps its key");
    let back = assign(
        BodyId(0),
        &lost.outputs,
        vec![piece(2.0, 6), piece(0.0, 0), piece(4.0, 12)],
        &mut counter(9),
    );
    assert!(back.served.iter().any(|(b, p)| b.0 == 5 && x_of(p) == 2.0));
}

/// The own body's piece gone: it takes the first unclaimed piece, so
/// features after it keep a body; a new piece gets a new body.
#[test]
fn the_own_body_always_takes_a_piece_and_new_pieces_get_new_bodies() {
    let first = assign(
        BodyId(0),
        &[],
        vec![piece(0.0, 0), piece(2.0, 6)],
        &mut counter(5),
    );
    let out = assign(
        BodyId(0),
        &first.outputs,
        vec![piece(2.0, 6), piece(8.0, 12), piece(9.0, 0).clone()],
        &mut counter(9),
    );
    // Piece at 9 carries the own body's key; 8 is new.
    let bodies: Vec<_> = out.served.iter().map(|(b, p)| (b.0, x_of(p))).collect();
    assert_eq!(bodies, vec![(0, 9.0), (5, 2.0), (9, 8.0)]);
    let orphaned = assign(
        BodyId(0),
        &first.outputs,
        vec![piece(2.0, 6), piece(8.0, 12)],
        &mut counter(9),
    );
    let bodies: Vec<_> = orphaned
        .served
        .iter()
        .map(|(b, p)| (b.0, x_of(p)))
        .collect();
    assert_eq!(bodies, vec![(0, 8.0), (5, 2.0)]);
}
