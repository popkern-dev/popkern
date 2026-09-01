//! PK-387 — extrudes of a full-circle profile whose end is UP TO a
//! face on ANOTHER body: the target names its host, the extrude resolves
//! against it, and a Boss merges into it.

mod circle_extrude_support;

use circle_extrude_support::*;
use popkern_core::naming::{EntityName, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_history::*;
use popkern_part::prelude::*;
use popkern_sketch::{AnalyticProfile, ProfileShape};

/// Plate 2.4 x 1.6 x 0.4 with a 0.4 x 0.4 pillar at a corner rising to
/// z = 1.4 (explicit union, as the recorder records tools). Returns the
/// tree and the plate body.
fn plate_with_pillar() -> (FeatureTree, BodyId) {
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Analytic(AnalyticProfile {
                shape: ProfileShape::Rect {
                    width: 2.4.into(),
                    height: 1.6.into(),
                },
                center: [0.0.into(), 0.0.into()],
            }),
        },
        "plate",
    );
    let h = add_param(&mut tree, "plate_h", 0.4);
    let plate_extrude = tree.add_feature(
        Feature::Extrude {
            face: name(1, RoleCategory::Generated, 0),
            direction: None,
            end_condition: ExtrudeEnd::Blind { distance: h },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "plate_extrude",
    );
    let plate_body = tree.feature(plate_extrude).expect("plate feature").body;
    let top = add_param(&mut tree, "plate_top", 0.4);
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::BaseOffset {
                base: BasePlane::XY,
                distance: top,
            },
            sketch_data: SketchPayload::Analytic(AnalyticProfile {
                shape: ProfileShape::Rect {
                    width: 0.4.into(),
                    height: 0.4.into(),
                },
                center: [0.9.into(), 0.5.into()],
            }),
        },
        "pillar_sketch",
    );
    let pillar_h = add_param(&mut tree, "pillar_h", 1.0);
    let pillar = tree.add_feature(
        Feature::Extrude {
            face: name(3, RoleCategory::Generated, 0),
            direction: None,
            end_condition: ExtrudeEnd::Blind { distance: pillar_h },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "pillar",
    );
    let pillar_body = tree.feature(pillar).expect("pillar feature").body;
    tree.add_feature_on(
        Feature::Boolean {
            operation: BoolOpKind::Union,
            tool_body: pillar_body,
        },
        "plate_plus_pillar",
        Some(plate_body),
    );
    tree.replay(Tolerance::default())
        .expect("plate + pillar replay");
    (tree, plate_body)
}

/// The served name of the planar face at height `z` on `body`.
fn planar_face_at(tree: &FeatureTree, body: BodyId, z: f64) -> EntityName {
    let brep = tree.body_brep(body).expect("body");
    brep.faces
        .iter()
        .filter(|(h, _)| {
            let pts = brep.face_boundary_sample_points(*h, 4);
            !pts.is_empty() && pts.iter().all(|p| (p.z - z).abs() < 1e-9)
        })
        .find_map(|(_, f)| f.name.and_then(|id| brep.name_pool.get(id).cloned()))
        .expect("planar face at height is named")
}

/// A Ø0.7 boss sketched at z = 0.4 and extruded UP up to the pillar's
/// top face on the plate body (the target names its host).
fn add_boss_up_to_pillar(tree: &mut FeatureTree, plate_body: BodyId, mode: OperationMode) {
    let pillar_top = planar_face_at(tree, plate_body, 1.4);
    let top = add_param(tree, "boss_z", 0.4);
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::BaseOffset {
                base: BasePlane::XY,
                distance: top,
            },
            sketch_data: analytic_circle(0.0, 0.0, 0.35),
        },
        "boss_sketch",
    );
    tree.add_feature(
        Feature::Extrude {
            face: name(6, RoleCategory::Generated, 0),
            direction: None,
            end_condition: ExtrudeEnd::UpToSurface { face: pillar_top },
            mode,
            draft_angle: None,
            second_direction: None,
        },
        "boss",
    );
}

#[test]
fn analytic_circle_up_to_surface_on_another_body_is_one_periodic_cylinder() {
    // NewBody: the extrude resolves its target on the plate body and is
    // served un-merged — the wire sweep body, exactly Part's.
    let (mut tree, plate_body) = plate_with_pillar();
    add_boss_up_to_pillar(&mut tree, plate_body, OperationMode::NewBody);
    let result = tree.replay(Tolerance::default()).expect("replay");
    assert!(
        result.report.failures().next().is_none(),
        "{:?}",
        result.report
    );
    let boss = tree.feature_by_label("boss").expect("boss feature");
    let boss_body = tree.feature(boss).expect("node").body;
    let body = tree.body_brep(boss_body).expect("boss body").clone();
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    assert_eq!(surface_census(&body), (2, 1, 0), "one wall up to z = 1.4");
    assert_eq!(one_edge_planar_loops(&body), 2);
    let part = Profile::on(Plane::xy_at(0.4))
        .circle(0.35)
        .extrude(1.0)
        .expect("boss");
    assert_same_body(
        "analytic circle up-to-surface (other body)",
        &body,
        part.brep(),
    );
}

#[test]
fn analytic_circle_up_to_surface_boss_merges_into_its_host() {
    // Boss: the same extrude merges into the body its target lives on.
    // (PK-394 restored Part's own three-solid coplanar union, so the
    // merged result is full Part parity.)
    let (mut tree, plate_body) = plate_with_pillar();
    add_boss_up_to_pillar(&mut tree, plate_body, OperationMode::Boss);
    let result = tree.replay(Tolerance::default()).expect("replay");
    assert!(
        result.report.failures().next().is_none(),
        "{:?}",
        result.report
    );
    let body = tree.body_brep(plate_body).expect("plate body").clone();
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    assert_eq!(
        surface_census(&body),
        (12, 1, 0),
        "plate + pillar + one boss wall"
    );
    let boss = tree.feature_by_label("boss").expect("boss feature");
    assert_eq!(
        tree.merge_host(boss),
        Some(plate_body),
        "boss merged into the plate"
    );
    let plate = Profile::on(Plane::XY)
        .rect(2.4, 1.6)
        .extrude(0.4)
        .expect("plate");
    let pillar = Profile::on(Plane::xy_at(0.4))
        .rect(0.4, 0.4)
        .at(0.9, 0.5)
        .extrude(1.0)
        .expect("pillar");
    let boss_tool = Profile::on(Plane::xy_at(0.4))
        .circle(0.35)
        .extrude(1.0)
        .expect("boss");
    let part = plate
        .union(&pillar)
        .expect("union")
        .union(&boss_tool)
        .expect("union");
    assert_same_body("boss up-to merged ≡ Part union", &body, part.brep());
}
