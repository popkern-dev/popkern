// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use super::*;
use crate::FeatureTree;
use crate::feature::Feature;
use crate::unit::Unit;
use popkern_core::naming::{EntityRole, OperationId, RoleCategory};

fn name(op: u64, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op),
        role: EntityRole {
            category: RoleCategory::Generated,
            index,
        },
        ancestors: Vec::new(),
    }
}

fn param(tree: &mut FeatureTree, label: &str, value: f64) -> ParameterId {
    tree.parameters_mut()
        .add_global(label, value, Unit::Meters)
        .expect("param")
}

fn full_recipe(end: ParameterId) -> FilletRecipeDesc {
    FilletRecipeDesc {
        radius_law: FilletRadiusLaw::Linear { end },
        continuity: ContinuityKind::G2,
        cap: CapKindDesc::Rounded,
        profile: BlendProfileKind::Conic { rho: 0.4 },
        interference: InterferencePolicyKind::ServeNonInterfering,
    }
}

/// The full recipe and every chamfer kind survive save/load
/// bit-for-bit (PK-361 acceptance: round-trips through `.pkn`).
#[test]
fn recipe_and_chamfer_kind_round_trip_through_pkn() {
    let mut tree = FeatureTree::new();
    let r = param(&mut tree, "r", 0.1);
    let r_end = param(&mut tree, "r_end", 0.2);
    let d = param(&mut tree, "d", 0.05);
    let d2 = param(&mut tree, "d2", 0.08);
    let ang = param(&mut tree, "ang", 0.6);
    let stations = FilletRecipeDesc {
        radius_law: FilletRadiusLaw::ByStations(vec![
            RadiusStation { t: 0.0, radius: r },
            RadiusStation {
                t: 0.5,
                radius: r_end,
            },
            RadiusStation { t: 1.0, radius: r },
        ]),
        ..Default::default()
    };
    tree.add_feature(
        Feature::Fillet {
            edges: vec![name(1, 0)],
            radius: r,
            recipe: full_recipe(r_end),
        },
        "FilletFull",
    );
    tree.add_feature(
        Feature::Fillet {
            edges: vec![name(1, 1)],
            radius: r,
            recipe: stations.clone(),
        },
        "FilletStations",
    );
    tree.add_feature(
        Feature::Chamfer {
            edges: vec![name(1, 2)],
            distance: d,
            kind: ChamferKindDesc::TwoDistance {
                face: name(1, 3),
                distance_2: d2,
            },
        },
        "ChamferTwo",
    );
    tree.add_feature(
        Feature::Chamfer {
            edges: vec![name(1, 4)],
            distance: d,
            kind: ChamferKindDesc::DistanceAngle {
                face: name(1, 5),
                angle: ang,
            },
        },
        "ChamferAngle",
    );

    let bytes = tree.save().expect("save");
    let loaded = FeatureTree::load(&bytes).expect("load");
    let feature =
        |i: usize| -> &Feature { &loaded.feature(loaded.timeline()[i]).expect("node").feature };

    match feature(0) {
        Feature::Fillet { recipe, .. } => assert_eq!(*recipe, full_recipe(r_end)),
        other => panic!("expected fillet, got {other:?}"),
    }
    match feature(1) {
        Feature::Fillet { recipe, .. } => assert_eq!(*recipe, stations),
        other => panic!("expected fillet, got {other:?}"),
    }
    match feature(2) {
        Feature::Chamfer { kind, .. } => assert_eq!(
            *kind,
            ChamferKindDesc::TwoDistance {
                face: name(1, 3),
                distance_2: d2
            }
        ),
        other => panic!("expected chamfer, got {other:?}"),
    }
    match feature(3) {
        Feature::Chamfer { kind, .. } => assert_eq!(
            *kind,
            ChamferKindDesc::DistanceAngle {
                face: name(1, 5),
                angle: ang
            }
        ),
        other => panic!("expected chamfer, got {other:?}"),
    }
}

/// Trees saved before the recipe existed (v3 envelopes with bare
/// `edges`/`radius` and `edges`/`distance` fields) load with the
/// default recipe and kind — the historical replay behaviour.
#[test]
fn pre_recipe_fillet_and_chamfer_json_load_with_defaults() {
    let mut tree = FeatureTree::new();
    let r = param(&mut tree, "r", 0.1);
    let d = param(&mut tree, "d", 0.05);
    // Author the pre-recipe JSON shape: serialize the modern feature,
    // then strip the fields that did not exist in older documents.
    let strip = |feature: &Feature, variant: &str, field: &str| -> Feature {
        let mut value = serde_json::to_value(feature).expect("to json");
        let body = value
            .get_mut(variant)
            .and_then(|v| v.as_object_mut())
            .expect("variant object");
        assert!(body.remove(field).is_some(), "{field} present before strip");
        serde_json::from_value(value).expect("pre-recipe json loads")
    };
    let fillet = strip(
        &Feature::Fillet {
            edges: vec![name(1, 0)],
            radius: r,
            recipe: full_recipe(r),
        },
        "Fillet",
        "recipe",
    );
    match fillet {
        Feature::Fillet { radius, recipe, .. } => {
            assert_eq!(radius, r);
            assert_eq!(recipe, FilletRecipeDesc::default());
        }
        other => panic!("expected fillet, got {other:?}"),
    }
    let chamfer = strip(
        &Feature::Chamfer {
            edges: vec![name(1, 1)],
            distance: d,
            kind: ChamferKindDesc::DistanceAngle {
                face: name(1, 2),
                angle: d,
            },
        },
        "Chamfer",
        "kind",
    );
    match chamfer {
        Feature::Chamfer { distance, kind, .. } => {
            assert_eq!(distance, d);
            assert_eq!(kind, ChamferKindDesc::Equal);
        }
        other => panic!("expected chamfer, got {other:?}"),
    }
    // A partially specified recipe fills the missing fields too.
    let partial: FilletRecipeDesc =
        serde_json::from_str(r#"{"continuity":"G2"}"#).expect("partial recipe");
    assert_eq!(partial.continuity, ContinuityKind::G2);
    assert_eq!(partial.radius_law, FilletRadiusLaw::Constant);
    assert_eq!(partial.cap, CapKindDesc::Open);
}

/// The default recipe maps to exactly `FilletOptions::default()` on
/// every document-visible knob, so replaying an old tree is
/// bit-identical to the pre-recipe hardcoded call.
#[test]
fn default_recipe_matches_kernel_defaults() {
    let opts = FilletRecipeDesc::default().to_kernel_options();
    let kernel = popkern_blend::FilletOptions::default();
    assert_eq!(opts.continuity, kernel.continuity);
    assert_eq!(opts.cap_kind, kernel.cap_kind);
    assert_eq!(opts.profile, kernel.profile);
    assert_eq!(opts.interference_policy, kernel.interference_policy);
    assert_eq!(opts.propagate_tangent, kernel.propagate_tangent);
    assert_eq!(opts.use_rational_gregory, kernel.use_rational_gregory);
}

#[test]
fn recipe_knobs_reach_kernel_options() {
    let mut tree = FeatureTree::new();
    let end = param(&mut tree, "end", 0.2);
    let opts = full_recipe(end).to_kernel_options();
    assert_eq!(opts.continuity, popkern_blend::ContinuityLevel::G2);
    assert_eq!(opts.cap_kind, popkern_blend::CapKind::Rounded);
    assert_eq!(
        opts.profile,
        popkern_blend::BlendProfile::Conic { rho: 0.4 }
    );
    assert_eq!(
        opts.interference_policy,
        popkern_blend::InterferencePolicy::ServeNonInterfering
    );
    assert_eq!(CapKindDesc::Flat.to_kernel(), popkern_blend::CapKind::Flat);
    assert_eq!(
        CapKindDesc::TangentMatch.to_kernel(),
        popkern_blend::CapKind::TangentMatch
    );
    assert_eq!(
        ContinuityKind::G0.to_kernel(),
        popkern_blend::ContinuityLevel::G0
    );
    assert_eq!(
        BlendProfileKind::CurvatureContinuous.to_kernel(),
        popkern_blend::BlendProfile::CurvatureContinuous
    );
}

/// Every parameter a recipe or kind references is visible to the
/// dependency/scope machinery, so a variable-radius end or a second
/// chamfer distance is a tracked driver like the primary value.
#[test]
fn recipe_and_kind_params_are_reported() {
    let mut tree = FeatureTree::new();
    let end = param(&mut tree, "end", 0.2);
    let s0 = param(&mut tree, "s0", 0.1);
    let s1 = param(&mut tree, "s1", 0.3);
    let d2 = param(&mut tree, "d2", 0.04);
    let ang = param(&mut tree, "ang", 0.5);
    assert_eq!(full_recipe(end).param_ids(), vec![("end_radius", end)]);
    let stations = FilletRecipeDesc {
        radius_law: FilletRadiusLaw::ByStations(vec![
            RadiusStation { t: 0.0, radius: s0 },
            RadiusStation { t: 1.0, radius: s1 },
        ]),
        ..Default::default()
    };
    assert_eq!(
        stations.param_ids(),
        vec![("station_radius", s0), ("station_radius", s1)]
    );
    assert!(FilletRecipeDesc::default().param_ids().is_empty());

    let two = ChamferKindDesc::TwoDistance {
        face: name(2, 0),
        distance_2: d2,
    };
    assert_eq!(two.param_ids(), vec![("distance_2", d2)]);
    assert_eq!(two.reference_face(), Some(&name(2, 0)));
    let angle = ChamferKindDesc::DistanceAngle {
        face: name(2, 1),
        angle: ang,
    };
    assert_eq!(angle.param_ids(), vec![("angle", ang)]);
    assert_eq!(angle.reference_face(), Some(&name(2, 1)));
    assert!(ChamferKindDesc::Equal.param_ids().is_empty());
    assert_eq!(ChamferKindDesc::Equal.reference_face(), None);

    // Through the feature-level collector as well.
    let r = param(&mut tree, "r", 0.1);
    let feature = Feature::Fillet {
        edges: vec![name(1, 0)],
        radius: r,
        recipe: full_recipe(end),
    };
    let ids: Vec<ParameterId> = crate::feature_params::feature_param_ids(&feature)
        .into_iter()
        .map(|(_, id)| id)
        .collect();
    assert_eq!(ids, vec![r, end]);
}
