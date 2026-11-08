use super::*;
use crate::{BoolOperation, BoolOptions, NonManifoldPolicy};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    let tol = Tolerance::default();
    BRep::unit_cube(tol)
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube")
}

fn two_apart() -> BRep {
    crate::concatenate_solids(cube_at(0.0, 0.0, 0.0), cube_at(3.0, 0.0, 0.0)).expect("two")
}

#[test]
fn each_piece_keeps_its_own_topology_and_names() {
    let both = two_apart();
    let pieces = piece_bodies(&both).expect("pieces");
    assert_eq!(pieces.len(), 2);
    for piece in &pieces {
        assert!(piece.validate().is_ok(), "{:?}", piece.validate());
        assert_eq!(
            (
                piece.solid_regions().len(),
                piece.faces.len(),
                piece.edges.len()
            ),
            (1, 6, 12)
        );
        assert_eq!(piece.vertices.len(), 8);
        for (h, face) in &piece.faces {
            assert_eq!(face.name, both.faces[h].name, "a kept face keeps its name");
        }
    }
    let xs: Vec<f64> = pieces
        .iter()
        .map(|p| {
            p.vertices
                .values()
                .map(|v| v.position.x)
                .fold(f64::MAX, f64::min)
        })
        .collect();
    assert!((xs[0] - xs[1]).abs() > 2.0, "two different pieces: {xs:?}");
}

#[test]
fn keeping_every_piece_or_none() {
    let both = two_apart();
    let all = retain_pieces(&both, &piece_handles(&both)).expect("all");
    assert_eq!(
        (
            all.solid_regions().len(),
            all.faces.len(),
            all.vertices.len()
        ),
        (2, 12, 16)
    );
    let none = retain_pieces(&both, &[]).expect("none");
    assert!(none.faces.is_empty() && none.vertices.is_empty());
}

#[test]
fn touching_split_bodies_are_each_valid_alone() {
    // Two cubes touching along one edge, served split: two bodies,
    // each a valid cube with its own copy of the contact edge.
    let tol = Tolerance::default();
    let options = BoolOptions {
        non_manifold: NonManifoldPolicy::SplitToManifold,
    };
    let bodies = crate::boolean_bodies(
        BoolOperation::Union,
        cube_at(0.0, 0.0, 0.0),
        cube_at(1.0, 1.0, 0.0),
        tol,
        None,
        &options,
    )
    .expect("edge contact serves as bodies");
    assert_eq!(bodies.len(), 2);
    for one in &bodies {
        assert_eq!(one.solid_regions().len(), 1);
        assert!(one.validate_all().is_empty(), "{:?}", one.validate_all());
    }
}

#[test]
fn solids_sharing_a_vertex_are_not_pieces() {
    let mut both = two_apart();
    let handles = piece_handles(&both);
    // Re-point one half-edge of the second solid at a vertex of the
    // first: the solids now share that vertex.
    let first = retain_pieces(&both, &handles[..1]).expect("first");
    let borrowed = first.vertices.keys().next().expect("a vertex");
    let second = retain_pieces(&both, &handles[1..]).expect("second");
    let victim = second.half_edges.keys().next().expect("a half-edge");
    both.half_edges[victim].origin = borrowed;
    let error = retain_pieces(&both, &handles[..1]).expect_err("shared vertex");
    assert!(matches!(error, BoolError::DegenerateInput(ref m) if m.contains("vertex")));
}
