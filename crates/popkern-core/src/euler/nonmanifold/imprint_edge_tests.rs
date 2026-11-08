use super::*;
use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::tolerance::Tolerance;
use crate::topology::BodyKind;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn assert_valid(body: &BRep, what: &str) {
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{what}: {errors:?}");
}

fn top_face(body: &BRep) -> FaceHandle {
    body.faces
        .keys()
        .find(|&f| {
            body.face_point_trim_containment(f, &Point3::new(0.5, 0.5, 1.0))
                .is_some_and(|p| p.distance < 1e-9)
        })
        .expect("top face")
}

fn line(p: Point3, q: Point3) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: p,
            direction: q - p,
        },
        t_start: 0.0,
        t_end: 1.0,
    }
}

/// The vertex splitting the cube's edge through `p` there.
fn split_at(body: &mut BRep, p: Point3) -> VertexHandle {
    let edge = body
        .edges
        .keys()
        .find(|&e| {
            body.edges[e].curve.as_ref().is_some_and(|c| {
                c.geometry
                    .closest_point(&p)
                    .is_some_and(|t| c.geometry.point_at(t).distance_to(&p) < 1e-9)
            }) && body.edge_vertices(e).is_some_and(|(a, b)| {
                let (a, b) = (body.vertices[a].position, body.vertices[b].position);
                (a.distance_to(&p) + p.distance_to(&b) - a.distance_to(&b)).abs() < 1e-9
            })
        })
        .expect("edge through the point");
    body.split_edge_at_point(edge, p, &tol())
        .expect("split")
        .vertex
}

/// An edge across the top face between points of two opposite edges
/// splits the face in two; killing it merges them again. Both valid.
#[test]
fn an_edge_across_a_face_splits_it_and_goes_again() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let (p, q) = (Point3::new(0.5, 0.0, 1.0), Point3::new(0.5, 1.0, 1.0));
    let (a, b) = (split_at(&mut cube, p), split_at(&mut cube, q));
    let top = top_face(&cube);
    let faces = cube.faces.len();
    let made = cube.imprint_edge(top, a, b, line(p, q)).expect("imprint");
    assert!(made.new_face.is_some());
    assert_eq!(cube.faces.len(), faces + 1);
    assert_eq!(cube.edge_vertices(made.edge), Some((a, b)));
    assert_eq!(cube.body_kind(), BodyKind::Solid);
    assert_valid(&cube, "imprinted across");
    cube.kill_imprinted_edge(made.edge).expect("kill");
    assert_eq!(cube.faces.len(), faces);
    assert_valid(&cube, "killed again");
}

/// An edge from the boundary to a vertex loop is a spur of the outer
/// loop, the vertex loop absorbed; killing it stands the far end as a
/// vertex loop again.
#[test]
fn an_edge_to_a_vertex_loop_is_a_spur() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let p = Point3::new(0.5, 0.0, 1.0);
    let c = Point3::new(0.5, 0.5, 1.0);
    let a = split_at(&mut cube, p);
    let top = top_face(&cube);
    let lone = cube.make_vertex_loop(top, c).expect("vertex loop");
    let loops = cube.loops.len();
    // From the vertex loop, so the canonical fin leaves it.
    let made = cube
        .imprint_edge(top, lone.vertex, a, line(c, p))
        .expect("imprint");
    assert_eq!(made.new_face, None);
    assert_eq!(cube.loops.len(), loops - 1);
    assert!(cube.vertex_loops().is_empty());
    let canonical = cube.edges[made.edge].half_edge;
    assert_eq!(cube.half_edges[canonical].origin, lone.vertex);
    assert_valid(&cube, "spur");
    cube.kill_imprinted_edge(made.edge).expect("kill");
    assert_eq!(cube.loops.len(), loops);
    assert_eq!(cube.vertex_loops().len(), 1);
    assert_valid(&cube, "spur killed");
}

/// An edge between two vertex loops is a slit: one inner loop of its
/// two fins. Killing it stands both ends as vertex loops again.
#[test]
fn an_edge_between_vertex_loops_is_a_slit() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let (p, q) = (Point3::new(0.3, 0.5, 1.0), Point3::new(0.7, 0.5, 1.0));
    let top = top_face(&cube);
    let a = cube.make_vertex_loop(top, p).expect("a");
    let b = cube.make_vertex_loop(top, q).expect("b");
    let made = cube
        .imprint_edge(top, a.vertex, b.vertex, line(p, q))
        .expect("imprint");
    assert_eq!(cube.faces[top].inner_loops.len(), 1);
    assert!(cube.vertex_loops().is_empty());
    assert_valid(&cube, "slit");
    cube.kill_imprinted_edge(made.edge).expect("kill");
    assert!(cube.faces[top].inner_loops.is_empty());
    assert_eq!(cube.vertex_loops().len(), 2);
    assert_valid(&cube, "slit killed");
}

/// A vertex off the face, and one vertex for both ends, decline with
/// the body unchanged; a face edge is not an imprinted one to kill
/// unless its faces share a surface.
#[test]
fn imprints_decline_off_the_face() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let top = top_face(&cube);
    let bottom = cube
        .vertices
        .iter()
        .find(|(_, v)| v.position.z.abs() < 1e-9)
        .map(|(h, _)| h)
        .expect("bottom vertex");
    let up = cube
        .vertices
        .iter()
        .find(|(_, v)| (v.position.z - 1.0).abs() < 1e-9)
        .map(|(h, _)| h)
        .expect("top vertex");
    let before = super::super::NonManifoldCounts::of(&cube);
    let seg = line(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0));
    assert!(matches!(
        cube.imprint_edge(top, bottom, up, seg.clone()),
        Err(EulerError::NotImprintable { .. })
    ));
    assert!(matches!(
        cube.imprint_edge(top, up, up, seg),
        Err(EulerError::NotImprintable { .. })
    ));
    let cube_edge = cube.edges.keys().next().expect("edge");
    assert!(matches!(
        cube.kill_imprinted_edge(cube_edge),
        Err(EulerError::NotAnImprintedEdge { .. })
    ));
    assert_eq!(super::super::NonManifoldCounts::of(&cube), before);
}
