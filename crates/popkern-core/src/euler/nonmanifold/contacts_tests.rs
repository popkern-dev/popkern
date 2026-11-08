use super::*;
use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::BodyKind;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn assert_valid(body: &BRep, what: &str) {
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{what}: {errors:?}");
}

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube")
}

/// Two cubes, one at the origin and one moved by `(x, y, z)`, glued at
/// every coincident vertex and edge: the general body they make.
fn glued(x: f64, y: f64, z: f64) -> BRep {
    let mut body = cube_at(0.0, 0.0, 0.0);
    body.append_body(&cube_at(x, y, z));
    body.split_and_glue_reference()
}

impl BRep {
    /// Glue every pair of coincident vertices, then every pair of
    /// coincident edges, by hand.
    fn split_and_glue_reference(mut self) -> BRep {
        loop {
            let keys: Vec<_> = self.vertices.keys().collect();
            let pair = keys.iter().enumerate().find_map(|(i, &a)| {
                keys[i + 1..]
                    .iter()
                    .copied()
                    .find(|&b| {
                        self.vertices[a]
                            .position
                            .distance_to(&self.vertices[b].position)
                            < 1e-9
                    })
                    .map(|b| (a, b))
            });
            let Some((a, b)) = pair else { break };
            self.glue_vertices(a, b).expect("vertex");
        }
        loop {
            let keys: Vec<_> = self.edges.keys().collect();
            let same = |e: crate::topology::EdgeHandle, f: crate::topology::EdgeHandle| {
                let (Some((a, b)), Some((c, d))) = (self.edge_vertices(e), self.edge_vertices(f))
                else {
                    return false;
                };
                (a == c && b == d) || (a == d && b == c)
            };
            let pair = keys.iter().enumerate().find_map(|(i, &a)| {
                keys[i + 1..]
                    .iter()
                    .copied()
                    .find(|&b| same(a, b))
                    .map(|b| (a, b))
            });
            let Some((a, b)) = pair else { break };
            self.glue_edges(a, b).expect("edge");
        }
        self
    }
}

/// `bodies` appended into one body and glued where they touch.
fn reglued(bodies: &[BRep]) -> BRep {
    let mut body = BRep::new(tol());
    for b in bodies {
        body.append_body(b);
    }
    body.glue_touching().expect("glue");
    body
}

/// Two cubes glued along an edge split into the two cubes as separate
/// bodies, each a valid solid holding nothing of the other; glued again
/// where they touch they are the general body.
#[test]
fn an_edge_contact_splits_into_bodies_and_glues_again() {
    let general = glued(1.0, 1.0, 0.0);
    assert_eq!(general.body_kind(), BodyKind::General);
    let bodies = general.split_to_manifold().expect("split");
    assert_eq!(bodies.len(), 2);
    for b in &bodies {
        assert_eq!(b.body_kind(), BodyKind::Solid);
        assert_eq!((b.faces.len(), b.edges.len(), b.vertices.len()), (6, 12, 8));
        assert!(b.contacts().is_empty());
        assert_valid(b, "piece");
    }
    let again = reglued(&bodies);
    assert_eq!(
        NonManifoldCounts::of(&again),
        NonManifoldCounts::of(&general)
    );
    assert_eq!(again.body_kind(), BodyKind::General);
    assert_valid(&again, "glued again");
}

/// Two cubes sharing a corner split at the vertex alone.
#[test]
fn a_vertex_contact_splits_into_bodies_and_glues_again() {
    let general = glued(1.0, 1.0, 1.0);
    let bodies = general.split_to_manifold().expect("split");
    assert_eq!(bodies.len(), 2);
    for b in &bodies {
        assert_valid(b, "piece");
    }
    assert_eq!(
        NonManifoldCounts::of(&reglued(&bodies)),
        NonManifoldCounts::of(&general)
    );
}

/// Taken apart in place, the pieces stay in one body touching without
/// sharing: an intermediate, which validation rejects (negative
/// control), coming apart into the same bodies.
#[test]
fn pieces_left_in_one_body_do_not_validate() {
    let general = glued(1.0, 1.0, 0.0);
    let mut pieces = general.clone();
    pieces.separate_at_contacts().expect("separate");
    assert_eq!(pieces.solid_regions().len(), 2);
    assert!(!pieces.validate_all().is_empty(), "touching lumps");
    assert_eq!(pieces.piece_bodies().expect("bodies").len(), 2);
}

/// A body stored with recorded contacts (written before they retired)
/// loads, the list ignored.
#[test]
fn an_old_document_with_recorded_contacts_loads() {
    let cube = cube_at(0.0, 0.0, 0.0);
    let mut value: serde_json::Value = serde_json::to_value(&cube).expect("save");
    value["topology_contacts"] = serde_json::json!([{
        "kind": "Vertex",
        "location": {"x": 1.0, "y": 1.0, "z": 1.0},
        "entities": []
    }]);
    let loaded: BRep = serde_json::from_value(value).expect("load");
    assert_eq!(NonManifoldCounts::of(&loaded), NonManifoldCounts::of(&cube));
    assert_valid(&loaded, "loaded");
}

/// A face between two solids has no manifold pieces: the split
/// declines and leaves the body as it was. (Gluing a facing pair's
/// edges alone meets coincident faces, a tie radial order declines;
/// `glue_faces` glues the pair in one step.)
#[test]
fn bodies_without_manifold_pieces_do_not_split() {
    let mut stacked = cube_at(0.0, 0.0, 0.0);
    stacked.append_body(&cube_at(0.0, 0.0, 1.0));
    loop {
        let keys: Vec<_> = stacked.vertices.keys().collect();
        let pair = keys.iter().enumerate().find_map(|(i, &a)| {
            keys[i + 1..]
                .iter()
                .copied()
                .find(|&b| {
                    stacked.vertices[a]
                        .position
                        .distance_to(&stacked.vertices[b].position)
                        < 1e-9
                })
                .map(|b| (a, b))
        });
        let Some((a, b)) = pair else { break };
        stacked.glue_vertices(a, b).expect("vertex");
    }
    let at_one: Vec<_> = stacked
        .faces
        .keys()
        .filter(|&f| {
            stacked
                .loop_cycle(stacked.faces[f].outer_loop)
                .iter()
                .all(|&h| {
                    (stacked.vertices[stacked.half_edges[h].origin].position.z - 1.0).abs() < 1e-9
                })
        })
        .collect();
    stacked.glue_faces(at_one[0], at_one[1]).expect("faces");
    let before = NonManifoldCounts::of(&stacked);
    assert!(matches!(
        stacked.split_to_manifold(),
        Err(EulerError::NotSplittable {
            reason: "a face lies between two solids"
        })
    ));
    assert!(stacked.separate_at_contacts().is_err());
    assert_eq!(NonManifoldCounts::of(&stacked), before);
}
