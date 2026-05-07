use super::*;

/// Smoke test: the Open rung must short-circuit cleanly without mutating topology.
#[test]
fn open_rung_is_a_clean_noop() {
    use crate::apex_construction::plan::plan_apex_resolution;
    use crate::apex_construction::plan_e2e_test_helpers::{
        build_fillet_surfaces, force_unsupported_kind, pick_apex_and_incident_edges, test_tolerance,
    };
    use crate::types::{FilletRadius, FilletSpec};

    let cube = BRep::unit_cube(test_tolerance()).expect("cube builds");
    let (apex, edges) = pick_apex_and_incident_edges(&cube);
    let radius = 0.1;
    let mut surfaces = build_fillet_surfaces(&cube, &edges, radius);
    force_unsupported_kind(&mut surfaces, 2);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    let plan = plan_apex_resolution(&cube, apex, &specs, &surfaces, &test_tolerance(), None)
        .expect("Open rung plan");
    assert!(plan.cap_kind.is_open(), "test setup expects Open plan");

    let mut brep = cube.clone();
    let pre_v = brep.vertices.len();
    let pre_e = brep.edges.len();
    let pre_f = brep.faces.len();

    let result = construct_apex_constructively(
        &mut brep,
        &plan,
        &[],
        &test_tolerance(),
        None,
        Default::default(),
    )
    .expect("Open rung construct succeeds");

    assert!(result.created_faces.is_empty());
    assert!(result.apex_vertices.is_empty());
    assert!(result.apex_edges.is_empty());
    assert_eq!(brep.vertices.len(), pre_v);
    assert_eq!(brep.edges.len(), pre_e);
    assert_eq!(brep.faces.len(), pre_f);
}
