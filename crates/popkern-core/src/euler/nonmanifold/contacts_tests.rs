use super::*;
use crate::geometry::Point3;
use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::{BodyKind, TopologyContactKind};

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

/// Two cubes glued along an edge split into the two cubes meeting at
/// one recorded edge contact and two vertex contacts, which validate;
/// gluing the contacts gives the general body back.
#[test]
fn an_edge_contact_splits_and_glues_again() {
    let general = glued(1.0, 1.0, 0.0);
    assert_eq!(general.body_kind(), BodyKind::General);
    let general_counts = NonManifoldCounts::of(&general);
    let mut pieces = general.clone();
    let contacts = pieces.split_to_manifold().expect("split");
    let kinds: Vec<_> = contacts.iter().map(|c| c.kind).collect();
    assert_eq!(
        kinds
            .iter()
            .filter(|&&k| k == TopologyContactKind::Edge)
            .count(),
        1
    );
    assert_eq!(
        kinds
            .iter()
            .filter(|&&k| k == TopologyContactKind::Vertex)
            .count(),
        2
    );
    assert_eq!(pieces.body_kind(), BodyKind::Solid);
    assert_eq!(pieces.topology_contacts.len(), 3);
    assert_valid(&pieces, "pieces with recorded contacts");
    let apart = {
        let mut b = cube_at(0.0, 0.0, 0.0);
        b.append_body(&cube_at(1.0, 1.0, 0.0));
        NonManifoldCounts::of(&b)
    };
    assert_eq!(NonManifoldCounts::of(&pieces), apart);

    let report = pieces.glue_contacts().expect("glue");
    assert_eq!((report.glued, report.kept.len()), (3, 0));
    assert!(pieces.topology_contacts.is_empty());
    assert_eq!(NonManifoldCounts::of(&pieces), general_counts);
    assert_eq!(pieces.body_kind(), BodyKind::General);
    assert_valid(&pieces, "glued again");
}

/// Two cubes sharing a corner split at the vertex alone.
#[test]
fn a_vertex_contact_splits_and_glues_again() {
    let general = glued(1.0, 1.0, 1.0);
    let mut pieces = general.clone();
    let contacts = pieces.split_to_manifold().expect("split");
    assert_eq!(contacts.len(), 1);
    assert_eq!(contacts[0].kind, TopologyContactKind::Vertex);
    assert_eq!(contacts[0].location, Point3::new(1.0, 1.0, 1.0));
    assert_valid(&pieces, "pieces");
    pieces.glue_contacts().expect("glue");
    assert_eq!(
        NonManifoldCounts::of(&pieces),
        NonManifoldCounts::of(&general)
    );
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
    assert_eq!(
        stacked.split_to_manifold(),
        Err(EulerError::NotSplittable {
            reason: "a face lies between two solids"
        })
    );
    assert_eq!(NonManifoldCounts::of(&stacked), before);
}
