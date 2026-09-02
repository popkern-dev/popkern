//! PK-407: every edge and vertex a served fillet, chamfer, or shell
//! leaves behind carries a stable name, so a query after a blend
//! records names (D3) and the document can spell them. Before this,
//! a filleted cube had 20 of 24 edges unnamed and any dress-up after
//! a fillet failed with "carry no stable name".

use popkern_core::naming::RoleCategory;
use popkern_script::prelude::*;

fn census(m: &mut Model) -> (usize, usize, usize, usize) {
    let brep = m.to_brep().expect("brep");
    let unnamed_edges = brep.edges.values().filter(|e| e.name.is_none()).count();
    let unnamed_vertices = brep.vertices.values().filter(|v| v.name.is_none()).count();
    (
        brep.edges.len(),
        unnamed_edges,
        brep.vertices.len(),
        unnamed_vertices,
    )
}

fn box_body(m: &mut Model, height: f64) -> Body {
    m.profile(popkern_script::Plane::XY)
        .rect(1.0, 1.0)
        .sketch("base")
        .expect("sketch")
        .extrude(m, End::Blind(height.into()))
        .expect("extrude")
}

#[test]
fn every_edge_and_vertex_after_a_fillet_is_named() {
    let mut m = Model::new();
    let body = box_body(&mut m, 0.75);
    body.fillet(&mut m, edges().straight().direction(Dir::Z), 0.1)
        .expect("fillet");
    let (edges_total, unnamed_edges, vertices_total, unnamed_vertices) = census(&mut m);
    assert_eq!(edges_total, 24);
    assert_eq!(unnamed_edges, 0, "unnamed edges after fillet");
    assert_eq!(vertices_total, 16);
    assert_eq!(unnamed_vertices, 0, "unnamed vertices after fillet");
}

#[test]
fn shortened_edges_keep_their_ancestor_and_new_edges_are_lateral() {
    let mut m = Model::new();
    let body = box_body(&mut m, 0.75);
    body.fillet(&mut m, edges().straight().direction(Dir::Z), 0.1)
        .expect("fillet");
    let brep = m.to_brep().expect("brep");
    let mut modified = 0;
    let mut lateral = 0;
    for edge in brep.edges.values() {
        let name = edge
            .name
            .and_then(|nid| brep.name_pool.get(nid))
            .expect("named");
        match name.role.category {
            RoleCategory::Modified => {
                assert_eq!(
                    name.ancestors.len(),
                    1,
                    "modified edge carries its ancestor"
                );
                modified += 1;
            }
            RoleCategory::Lateral => lateral += 1,
            other => panic!("unexpected edge role {other:?}"),
        }
    }
    // Eight top/bottom edges were shortened by the corner fillets;
    // the four vertical edges were replaced by strips (2 rails each =
    // 8 rails) plus 8 rim arcs.
    assert_eq!(modified, 8, "shortened box edges stay Modified");
    assert_eq!(lateral, 16, "strip rails and rims are Lateral");
}

#[test]
fn a_query_after_a_fillet_resolves_and_spells_its_names() {
    let mut m = Model::new();
    let body = box_body(&mut m, 0.75);
    body.fillet(&mut m, edges().straight().direction(Dir::Z), 0.1)
        .expect("fillet");
    let brep = m.to_brep().expect("brep");
    // The record path: resolve → names (D3), no "carry no stable name".
    let top = m
        .resolve_edges_in(&brep, &edges().straight().at_height(0.75))
        .expect("top edges resolve after the fillet");
    assert_eq!(top.names.len(), 4);
    for name in &top.names {
        assert_eq!(name.role.category, RoleCategory::Modified);
        let spelled = m.spell_entity(name).expect("spellable");
        assert!(spelled.starts_with("fillet-1/Modified/"), "{spelled}");
        assert_eq!(&m.parse_entity(&spelled).expect("parse"), name);
    }
    let rims = m
        .resolve_edges_in(&brep, &edges().circular())
        .expect("rim arcs resolve");
    assert_eq!(rims.names.len(), 8);
    assert!(
        rims.names
            .iter()
            .all(|n| n.role.category == RoleCategory::Lateral)
    );
    // The recorded document replays to the same named body.
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("blend_output_names");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("posts.pkn");
    m.save(&path).expect("save");
    let mut again = Model::open(&path).expect("open");
    let report = again.replay().expect("replay");
    assert_eq!(report.failures().count(), 0);
    let (_, unnamed_edges, _, unnamed_vertices) = census(&mut again);
    assert_eq!((unnamed_edges, unnamed_vertices), (0, 0));
}

#[test]
fn a_shell_leaves_every_edge_and_vertex_named() {
    let mut m = Model::new();
    let body = box_body(&mut m, 1.0);
    body.shell(
        &mut m,
        faces().facing(Dir::Z),
        0.1,
        ShellDirectionKind::Inward,
    )
    .expect("shell");
    let (_, unnamed_edges, _, unnamed_vertices) = census(&mut m);
    assert_eq!(unnamed_edges, 0, "unnamed edges after shell");
    assert_eq!(unnamed_vertices, 0, "unnamed vertices after shell");
}
