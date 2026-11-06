use crate::geometry::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn two_disjoint_cubes_append_into_one_valid_body_of_two_regions() {
    let tol = Tolerance::default();
    let mut a = BRep::unit_cube(tol).expect("cube");
    let b = BRep::unit_cube(tol)
        .and_then(|c| c.transform(&Transform3::translation(3.0, 0.0, 0.0)))
        .expect("moved cube");
    let map = a.append_body(&b);
    assert_eq!(a.faces.len(), 12);
    assert_eq!(a.solid_regions().len(), 2);
    assert_eq!(a.back_shell_count(), 2);
    assert_eq!(map.faces.len(), 6);
    assert!(a.validate().is_ok(), "{:?}", a.validate());
    for (&old, &new) in &map.faces {
        assert_eq!(
            a.shell_region(a.face_back_shell(new).expect("shell")),
            map.regions
                .get(
                    &b.shell_region(b.face_back_shell(old).expect("shell"))
                        .expect("region")
                )
                .copied()
        );
    }
}

/// Every ring link is remapped, so a ring of three survives the copy
/// (the copy it replaces kept only a mate).
#[test]
fn a_ring_of_three_fins_survives_an_append() {
    let tol = Tolerance::default();
    let mut src = BRep::unit_cube(tol).expect("cube");
    let edge = src.edges.keys().next().expect("an edge");
    src.thread_test_third_fin(edge);
    let mut into = BRep::new(tol);
    let map = into.append_body(&src);
    assert_eq!(into.edge_fin_count(map.edges[&edge]), 3);
}
