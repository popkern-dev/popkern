use crate::constructors::BRepBuilder;
use crate::geometry::Point3;
use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, Vertex, VertexUse};

/// A polyline wire through `points` (open), built by the builder.
fn polyline(points: &[(f64, f64, f64)]) -> BRep {
    let mut b = BRepBuilder::new(Tolerance::default());
    let vs: Vec<_> = points
        .iter()
        .map(|&(x, y, z)| b.add_vertex(Point3::new(x, y, z)))
        .collect();
    for w in vs.windows(2) {
        b.add_edge(w[0], w[1]);
    }
    b.build_unchecked()
}

/// A built wire is one shell of the outside holding its edges; each
/// inner vertex has two wire uses, each end one.
#[test]
fn a_built_wire_lies_in_one_shell_of_the_outside() {
    let wire = polyline(&[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0)]);
    assert_eq!(wire.shells.len(), 1);
    let (shell, data) = wire.shells.iter().next().expect("a shell");
    assert_eq!(data.region, wire.outside);
    assert_eq!(data.wire_edges.len(), 2);
    for e in wire.wire_edges() {
        assert_eq!(wire.wire_edge_shell(e), Some(shell));
    }
    assert!(wire.validate_all().is_empty(), "{:?}", wire.validate_all());
    let uses: Vec<usize> = wire.vertex_uses_all().values().map(Vec::len).collect();
    assert_eq!(uses, vec![1, 2, 1]);
    assert!(
        wire.vertex_uses_all()
            .values()
            .flatten()
            .all(|u| matches!(u, VertexUse::Wire(_)))
    );
    assert_eq!(wire.back_shell_count(), 0);
}

/// Two wires apart in one body are two shells; appended into another
/// body they keep their shells.
#[test]
fn disjoint_wires_are_shells_of_their_own() {
    let mut a = polyline(&[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]);
    let b = polyline(&[(0.0, 5.0, 0.0), (1.0, 5.0, 0.0), (2.0, 5.0, 0.0)]);
    a.append_body(&b);
    assert_eq!(a.shells.len(), 2);
    assert!(a.validate_all().is_empty(), "{:?}", a.validate_all());
    assert!(a.settle_wire_edges().is_empty(), "appended wires stay held");
}

/// A wire edge no shell holds fails validation (the negative control),
/// and settling holds it again.
#[test]
fn an_unheld_wire_edge_fails_validation_until_settled() {
    let mut wire = polyline(&[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]);
    wire.shells.clear();
    assert!(!wire.validate_all().is_empty());
    assert_eq!(wire.settle_wire_edges().len(), 1);
    assert!(wire.validate_all().is_empty(), "{:?}", wire.validate_all());
}

/// Settling drops a held edge that is gone, and the shell it emptied.
#[test]
fn settling_drops_a_removed_wire_edge_and_its_shell() {
    let mut wire = polyline(&[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]);
    let edge = wire.wire_edges()[0];
    wire.half_edges.retain(|_, fin| fin.edge != edge);
    wire.edges.remove(edge);
    assert!(wire.settle_wire_edges().is_empty());
    assert!(wire.shells.is_empty());
}

/// A wire body written before wires were held loads with its shell.
#[test]
fn a_wire_written_before_wire_shells_loads_held() {
    let mut wire = polyline(&[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]);
    wire.shells.clear();
    let json = serde_json::to_string(&wire).expect("serializes");
    let loaded: BRep = serde_json::from_str(&json).expect("loads");
    assert_eq!(loaded.shells.len(), 1);
    assert!(
        loaded.validate_all().is_empty(),
        "{:?}",
        loaded.validate_all()
    );
}

/// Every vertex of a cube has one use, a fan.
#[test]
fn a_cubes_vertices_each_have_one_fan() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    for uses in cube.vertex_uses_all().values() {
        assert!(matches!(uses.as_slice(), [VertexUse::Fan(_)]), "{uses:?}");
    }
    let v = cube.vertices.keys().next().expect("a vertex");
    assert_eq!(cube.vertex_uses(v), cube.vertex_uses_all()[&v]);
    assert!(!cube.is_nonmanifold_vertex(v));
}

/// Two cubes sharing a corner vertex: that vertex has two fans.
#[test]
fn cubes_sharing_a_corner_make_a_two_fan_vertex() {
    let mut a = BRep::unit_cube(Tolerance::default()).expect("cube");
    let b = BRep::unit_cube(Tolerance::default())
        .and_then(|b| b.transform(&Transform3::translation(1.0, 1.0, 1.0)))
        .expect("moved cube");
    a.append_body(&b);
    let at = |p: Point3| {
        a.vertices
            .iter()
            .filter(|(_, v)| v.position.distance_to(&p) < 1e-12)
            .map(|(h, _)| h)
            .collect::<Vec<_>>()
    };
    let corner = at(Point3::new(1.0, 1.0, 1.0));
    assert_eq!(corner.len(), 2);
    let (keep, gone) = (corner[0], corner[1]);
    for fin in a.half_edges.values_mut() {
        if fin.origin == gone {
            fin.origin = keep;
        }
    }
    a.vertices.remove(gone);
    let uses = a.vertex_uses(keep);
    assert_eq!(uses.len(), 2, "{uses:?}");
    assert!(uses.iter().all(|u| matches!(u, VertexUse::Fan(_))));
    assert!(a.is_nonmanifold_vertex(keep));
}

/// An acorn is a vertex held by a shell; an orphan vertex is not one
/// and fails validation.
#[test]
fn an_acorn_is_held_and_an_orphan_fails() {
    let mut body = BRep::new(Tolerance::default());
    let v = body.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    assert!(
        !body.validate_all().is_empty(),
        "an orphan vertex is refused"
    );
    let outside = body.outside;
    let shell = body.add_acorn_shell(outside, v).expect("an acorn");
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    assert_eq!(body.vertex_uses(v), vec![VertexUse::Isolated(shell)]);
    assert_eq!(body.acorns(), vec![(shell, v)]);
    assert!(body.add_acorn_shell(outside, v).is_none(), "held once");
}
