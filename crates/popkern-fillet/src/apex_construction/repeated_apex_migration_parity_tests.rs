use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::repeated_apex_side_support::SideSupportRegionDiagnosticSnapshot;
use super::repeated_apex_trim_subcomplex::TrimRewriteSubcomplexBuilder;

#[test]
fn migration_matrix_apex_snapshot_label_contract() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let support_face = faces[0];
    let cylinder_face = faces[1];
    let support_vertices = outer_loop_vertices(&brep, support_face);
    let mut plan = TrimRewriteSubcomplexBuilder::new();
    plan.push_new_generator_trim(
        0,
        cylinder_face,
        support_face,
        support_vertices[0],
        support_vertices[1],
    );
    let plan = plan.finish();

    let snapshot = SideSupportRegionDiagnosticSnapshot::build_read_only(&brep, &plan);
    let label = snapshot.diagnostic_label();

    assert!(snapshot.validate_fresh(&brep));
    assert!(label.contains("knowledge fingerprint"), "{label}");
    assert!(label.contains("trim_endpoint_agreements="), "{label}");
    assert!(label.contains("selected_plan"), "{label}");
    assert!(label.contains("ready_for_executor"), "{label}");
    assert!(label.contains("total_delta=dV=0 dE=1 dF=1 dL=1"), "{label}");
}

fn outer_loop_vertices(brep: &BRep, face: FaceHandle) -> Vec<VertexHandle> {
    let start = brep.loops[brep.faces[face].outer_loop].half_edge;
    let mut vertices = Vec::new();
    let mut cursor = start;
    loop {
        vertices.push(brep.half_edges[cursor].origin);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            break;
        }
    }
    vertices
}
