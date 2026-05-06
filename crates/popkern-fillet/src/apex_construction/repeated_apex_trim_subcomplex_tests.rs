use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::*;

fn cube_handles() -> (FaceHandle, FaceHandle, VertexHandle, VertexHandle) {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let vertices = brep.vertices.keys().collect::<Vec<_>>();
    (faces[0], faces[1], vertices[0], vertices[1])
}

#[test]
fn subcomplex_plan_owns_support_group_derivation() {
    let (cylinder_face, support_face, start, end) = cube_handles();
    let plan = TrimRewriteSubcomplexPlan::from_edge_plans(
        vec![TrimRewriteEdgePlan::new_generator_trim(
            0,
            cylinder_face,
            support_face,
            start,
            end,
        )],
        TopologyCountDelta::zero(),
    );

    assert_eq!(plan.support_groups().len(), 1);
    assert_eq!(plan.support_groups()[0].support_face(), support_face);
    assert_eq!(plan.support_groups()[0].generated_trim_count(), 1);
    assert_eq!(plan.support_groups()[0].cylinder_faces(), &[cylinder_face]);
}

#[test]
fn planned_generator_support_face_requires_unique_owner() {
    let (cylinder_face, support_face, start, end) = cube_handles();
    let edge_plan =
        TrimRewriteEdgePlan::new_generator_trim(0, cylinder_face, support_face, start, end);
    let plan =
        TrimRewriteSubcomplexPlan::from_edge_plans(vec![edge_plan], TopologyCountDelta::zero());

    assert_eq!(
        plan.planned_generator_support_face(cylinder_face, 0, start, end, "generator")
            .expect("unique owner"),
        support_face
    );

    let duplicate = TrimRewriteSubcomplexPlan::from_edge_plans(
        vec![edge_plan, edge_plan],
        TopologyCountDelta::zero(),
    );
    assert!(
        duplicate
            .planned_generator_support_face(cylinder_face, 0, start, end, "generator")
            .is_err()
    );
}

#[test]
fn subcomplex_builder_owns_generator_delta_and_support_groups() {
    let (cylinder_face, support_face, start, end) = cube_handles();
    let mut builder = TrimRewriteSubcomplexBuilder::new();
    builder.push_new_generator_trim(0, cylinder_face, support_face, start, end);
    let plan = builder.finish();

    assert_eq!(plan.edge_plans().len(), 1);
    assert_eq!(plan.support_groups().len(), 1);
    assert!(
        plan.diagnostic_label()
            .contains("edge-only delta dV=0 dE=1 dF=0 dL=0")
    );
}

#[test]
fn subcomplex_plan_projects_seed_faces_and_trim_intents() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let vertices = brep.vertices.keys().collect::<Vec<_>>();
    let cylinder_face = faces[0];
    let cap_face = faces[1];
    let support_face = faces[2];
    let start = vertices[0];
    let end = vertices[1];

    let mut builder = TrimRewriteSubcomplexBuilder::new();
    builder.push_existing_seam(0, cylinder_face, cap_face, start, end);
    builder.push_new_generator_trim(0, cylinder_face, support_face, start, end);
    let plan = builder.finish();

    let seed_faces = plan
        .seed_faces()
        .into_iter()
        .collect::<std::collections::BTreeSet<_>>();
    assert!(seed_faces.contains(&cylinder_face));
    assert!(seed_faces.contains(&cap_face));
    assert!(seed_faces.contains(&support_face));

    let trim_intents = plan.trim_intents();
    assert_eq!(trim_intents.len(), 2);
    let generator_intent = trim_intents
        .iter()
        .find(|intent| intent.is_new_generator_trim)
        .expect("generator intent is projected");
    assert_eq!(generator_intent.cylinder_face, cylinder_face);
    assert_eq!(generator_intent.support_face, Some(support_face));
    assert_eq!(generator_intent.start, start);
    assert_eq!(generator_intent.end, end);
}
