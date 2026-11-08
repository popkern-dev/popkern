use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, TopologyContactKind};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Two cubes glued on an edge: one edge contact and its two end
/// vertices; split into separate bodies, neither touches itself.
#[test]
fn a_general_body_lists_its_contacts_and_its_pieces_none() {
    let mut pieces = BRep::unit_cube(tol()).expect("cube");
    pieces.append_body(
        &BRep::unit_cube(tol())
            .and_then(|c| c.transform(&Transform3::translation(1.0, 1.0, 0.0)))
            .expect("cube"),
    );
    // Glue by hand.
    let mut general = pieces.clone();
    loop {
        let keys: Vec<_> = general.vertices.keys().collect();
        let pair = keys.iter().enumerate().find_map(|(i, &a)| {
            keys[i + 1..]
                .iter()
                .copied()
                .find(|&b| {
                    general.vertices[a]
                        .position
                        .distance_to(&general.vertices[b].position)
                        < 1e-9
                })
                .map(|b| (a, b))
        });
        let Some((a, b)) = pair else { break };
        general.glue_vertices(a, b).expect("vertex");
    }
    let shared: Vec<_> = general
        .edges
        .keys()
        .filter(|&e| {
            general.edge_vertices(e).is_some_and(|(a, b)| {
                [a, b].iter().all(|&v| {
                    let p = general.vertices[v].position;
                    (p.x - 1.0).abs() < 1e-9 && (p.y - 1.0).abs() < 1e-9
                })
            })
        })
        .collect();
    general.glue_edges(shared[0], shared[1]).expect("edge");
    let split = general.split_to_manifold().expect("split");
    let kinds = |b: &BRep| {
        let c = b.contacts();
        (
            c.iter()
                .filter(|c| c.kind == TopologyContactKind::Edge)
                .count(),
            c.iter()
                .filter(|c| c.kind == TopologyContactKind::Vertex)
                .count(),
        )
    };
    assert_eq!(kinds(&general), (1, 2));
    // Separate bodies touch no part of themselves.
    assert!(split.iter().all(|b| kinds(b) == (0, 0)));
    assert!(BRep::unit_cube(tol()).expect("cube").contacts().is_empty());
}
