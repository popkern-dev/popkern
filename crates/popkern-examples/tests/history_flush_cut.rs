//! PK-393 — the same flush through-cut built by the feature tree: a
//! circle sketched ON the plate top and cut DOWN (blind or up to the
//! bottom face), so the tool is swept in the opposite direction from
//! the Part matrix. Both sketch dialects.

mod circle_extrude_support;

use circle_extrude_support::*;
use popkern_core::naming::RoleCategory;
use popkern_core::tolerance::Tolerance;
use popkern_history::*;

fn flush_cut(payload: SketchPayload, up_to: bool) -> popkern_core::topology::BRep {
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Analytic(popkern_sketch::AnalyticProfile {
                shape: popkern_sketch::ProfileShape::Rect {
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
    tree.replay(Tolerance::default()).expect("plate replay");
    let bottom = find_face_name(&tree, 2, RoleCategory::Modified);
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::OnFace(find_face_name(&tree, 2, RoleCategory::Cap)),
            sketch_data: payload,
        },
        "hole_sketch",
    );
    let d = add_param(&mut tree, "depth", 0.4);
    let end = if up_to {
        ExtrudeEnd::UpToSurface { face: bottom }
    } else {
        ExtrudeEnd::Blind { distance: d }
    };
    tree.add_feature(
        Feature::Extrude {
            face: name(3, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, -1.0])),
            end_condition: end,
            mode: OperationMode::Cut,
            draft_angle: None,
            second_direction: None,
        },
        "hole",
    );
    let result = tree.replay(Tolerance::default()).expect("replay");
    assert!(
        result.report.failures().next().is_none(),
        "{:?}",
        result.report
    );
    tree.body_brep(plate_body).expect("plate body").clone()
}

#[test]
fn down_swept_flush_cut_serves_the_holed_plate() {
    for (label, payload, up_to) in [
        ("analytic blind", analytic_circle(0.0, 0.0, 0.35), false),
        (
            "analytic up-to-bottom",
            analytic_circle(0.0, 0.0, 0.35),
            true,
        ),
        (
            "constraint blind",
            SketchPayload::Constraint(circle_sketch_json(0.0, 0.0, 0.35)),
            false,
        ),
        (
            "constraint up-to-bottom",
            SketchPayload::Constraint(circle_sketch_json(0.0, 0.0, 0.35)),
            true,
        ),
    ] {
        let body = flush_cut(payload, up_to);
        assert!(
            body.validate_all().is_empty(),
            "{label}: {:?}",
            body.validate_all()
        );
        assert_eq!(surface_census(&body), (6, 1, 0), "{label}: one hole wall");
        assert_eq!(body.faces.len(), 7, "{label}");
    }
}
