use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::{Constraint, Sketch};

use crate::feature::{ExtrudeEnd, Feature, OperationMode, SketchPlaneDesc, Vec3Desc};
use crate::unit::Unit;
use crate::{CacheStrategy, FeatureTree};

fn entity_name(op_id: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

fn rectangle_sketch_json() -> String {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0);
    let p3 = sk.add_point(0.0, 1.0);
    sk.add_line(p0, p1);
    sk.add_line(p1, p2);
    sk.add_line(p2, p3);
    sk.add_line(p3, p0);
    sk.add_constraint(Constraint::Fixed(p0));
    sk.add_constraint(Constraint::Fixed(p1));
    sk.add_constraint(Constraint::Fixed(p2));
    sk.add_constraint(Constraint::Fixed(p3));
    sk.solve().expect("rectangle should solve");
    serde_json::to_string(&sk).expect("rectangle should serialize")
}

/// Unit-level positive path: hand-authored cached states around a real
/// fillet (bypassing replay — the replay-built extrude body ships
/// INVALID lateral surfaces today, PK-316, so no fillet feature has
/// ever executed through replay). The route must refresh the cached
/// post-fillet state to the new radius and match a fresh fillet.
#[test]
fn fillet_radius_tweak_refreshes_cached_state() {
    let tol = Tolerance::default();
    let mut state_before = BRep::unit_cube(tol).expect("cube");
    let edge = state_before.edges.keys().next().expect("cube edge");
    let edge_name = entity_name(7, RoleCategory::Generated, 0);
    let name_id = state_before.name_pool.intern(edge_name.clone(), None);
    state_before.edges[edge].name = Some(name_id);

    let fillet_at = |radius: f64| -> BRep {
        popkern_blend::fillet_edges(
            &state_before,
            &[popkern_blend::FilletSpec {
                edge,
                radius: popkern_blend::FilletRadius::Constant(radius),
            }],
            &popkern_blend::FilletOptions::default(),
            tol,
            None,
        )
        .expect("fillet serves")
        .brep
    };
    let state_at = fillet_at(0.1);

    // Tree shell: dummy ancestor + the fillet feature, cache primed
    // with the two states the route reads.
    let mut tree = FeatureTree::new();
    tree.set_cache_strategy(CacheStrategy::EveryStep);
    let anchor = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "anchor",
    );
    let r_old = tree
        .parameters_mut()
        .add_global("r", 0.1, Unit::Meters)
        .expect("param");
    let r_new = tree
        .parameters_mut()
        .add_global("r2", 0.2, Unit::Meters)
        .expect("param");
    let old_feature = Feature::Fillet {
        edges: vec![edge_name.clone()],
        radius: r_old,
        recipe: Default::default(),
    };
    let fillet_id = tree.add_feature(old_feature.clone(), "fillet");
    tree.cache
        .store_explicit(None, anchor, state_before.clone());
    tree.cache.store_explicit(None, fillet_id, state_at.clone());

    let new_feature = Feature::Fillet {
        edges: vec![edge_name],
        radius: r_new,
        recipe: Default::default(),
    };
    let refreshed = tree
        .try_parameter_tweak_in_place(fillet_id, &old_feature, &new_feature)
        .expect("radius tweak takes the in-place route");
    refreshed.validate().expect("refreshed state validates");

    let fresh = fillet_at(0.2);
    assert_eq!(refreshed.faces.len(), fresh.faces.len());
    assert_eq!(refreshed.vertices.len(), fresh.vertices.len());
    let corners = |b: &BRep| -> Vec<[f64; 3]> {
        let mut pts: Vec<[f64; 3]> = b
            .vertices
            .values()
            .map(|v| [v.position.x, v.position.y, v.position.z])
            .collect();
        pts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        pts
    };
    for (g, w) in corners(&refreshed).iter().zip(corners(&fresh).iter()) {
        for k in 0..3 {
            assert!(
                (g[k] - w[k]).abs() < 1.0e-9,
                "route-refreshed vertex {g:?} deviates from fresh-fillet vertex {w:?}"
            );
        }
    }
}

/// An edit that changes the EDGE SET is not a parameter tweak — the
/// route declines and `modify_feature` falls back to full
/// invalidation (cache entry at the feature evicted).
#[test]
fn edge_set_change_falls_back_to_full_replay() {
    let mut tree = FeatureTree::new();
    tree.set_cache_strategy(CacheStrategy::EveryStep);
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "sketch",
    );
    let height = tree
        .parameters_mut()
        .add_global("height", 1.0, Unit::Meters)
        .expect("height param");
    tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: height },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "extrude",
    );
    let r = tree
        .parameters_mut()
        .add_global("r", 0.1, Unit::Meters)
        .expect("radius param");
    let edge_name = entity_name(2, RoleCategory::Lateral, 0);
    let fillet = tree.add_feature(
        Feature::Fillet {
            edges: vec![edge_name.clone()],
            radius: r,
            recipe: Default::default(),
        },
        "fillet",
    );
    tree.replay(Tolerance::default()).expect("replay");

    // Same radius, EMPTY edge set — not a parameter tweak.
    tree.modify_feature(
        fillet,
        Feature::Fillet {
            edges: vec![],
            radius: r,
            recipe: Default::default(),
        },
    );
    assert!(
        tree.cache.get(None, fillet).is_none(),
        "a non-tweak edit must fall back to full invalidation"
    );
}

/// PK-316 acceptance: the sketch→extrude replay body validates, a
/// fillet feature EXECUTES through replay (impossible before — every
/// lateral face carried the base plane's surface), and a radius tweak
/// then takes the in-place route end-to-end through modify_feature.
#[test]
fn replay_extrude_validates_and_fillet_routes_in_place() {
    let mut tree = FeatureTree::new();
    tree.set_cache_strategy(CacheStrategy::EveryStep);
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "sketch",
    );
    let height = tree
        .parameters_mut()
        .add_global("height", 1.0, Unit::Meters)
        .expect("height param");
    let extrude = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: height },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "extrude",
    );
    tree.replay(Tolerance::default()).expect("base replay");
    {
        let timeline = tree.timeline().to_vec();
        let sketch_state = tree.cache.get(None, timeline[0]).expect("sketch cached");
        sketch_state
            .validate()
            .expect("replay-built sketch sheet validates (PK-317)");
    }
    let extruded = tree.cache.get(None, extrude).expect("extrude cached");
    extruded
        .validate()
        .expect("replay-built extrude body validates (PK-316)");

    // A named lateral edge with a genuine dihedral.
    let edge_name = {
        let state = extruded;
        state
            .name_pool
            .names
            .iter()
            .find(|n| {
                crate::resolve::resolve_edge(n, state).is_ok_and(|eh| {
                    let ed = &state.edges[eh];
                    if ed.curve.is_none() {
                        return false;
                    }
                    let hed = &state.half_edges[ed.half_edge];
                    let twin = &state.half_edges[hed.twin];
                    let nf = |lh| {
                        let f = state.loops.get(lh).map(|l| l.face)?;
                        match &state.faces.get(f)?.surface.as_ref()?.geometry {
                            popkern_core::geometry::surface::SurfaceGeometry::Plane {
                                normal,
                                ..
                            } => Some(*normal),
                            _ => None,
                        }
                    };
                    match (nf(hed.loop_), nf(twin.loop_)) {
                        (Some(a), Some(b)) => a.dot(&b).abs() < 0.9,
                        _ => false,
                    }
                })
            })
            .cloned()
            .expect("extruded body carries a filletable named edge")
    };

    let r_old = tree
        .parameters_mut()
        .add_global("r", 0.1, Unit::Meters)
        .expect("param");
    let fillet = tree.add_feature(
        Feature::Fillet {
            edges: vec![edge_name.clone()],
            radius: r_old,
            recipe: Default::default(),
        },
        "fillet",
    );
    let rr = tree.replay(Tolerance::default()).expect("fillet replay");
    assert!(
        matches!(
            rr.feature_results.get(&fillet),
            Some(crate::FeatureResult::Success { .. })
        ),
        "fillet feature executes through replay (PK-316), got {:?}",
        rr.feature_results.get(&fillet)
    );

    // The end-to-end route: modify_feature with a new radius must
    // refresh the cache entry at the fillet in place.
    let r_new = tree
        .parameters_mut()
        .add_global("r2", 0.2, Unit::Meters)
        .expect("param");
    tree.modify_feature(
        fillet,
        Feature::Fillet {
            edges: vec![edge_name],
            radius: r_new,
            recipe: Default::default(),
        },
    );
    assert!(
        tree.cache.get(None, fillet).is_some(),
        "radius tweak takes the in-place route through modify_feature"
    );
    tree.replay(Tolerance::default())
        .expect("post-tweak replay");
    let final_brep = tree.current_brep.clone().expect("final body");
    final_brep.validate().expect("tweaked body validates");
}

/// A recipe change (here G1 → G2) is not a radius tweak: the in-place
/// radius-modify route must decline so the tree falls back to a full
/// replay with the new options honoured (PK-361).
#[test]
fn fillet_recipe_change_declines_in_place_route() {
    let tol = Tolerance::default();
    let mut state_before = BRep::unit_cube(tol).expect("cube");
    let edge = state_before.edges.keys().next().expect("cube edge");
    let edge_name = entity_name(7, RoleCategory::Generated, 0);
    let name_id = state_before.name_pool.intern(edge_name.clone(), None);
    state_before.edges[edge].name = Some(name_id);
    let state_at = popkern_blend::fillet_edges(
        &state_before,
        &[popkern_blend::FilletSpec {
            edge,
            radius: popkern_blend::FilletRadius::Constant(0.1),
        }],
        &popkern_blend::FilletOptions::default(),
        tol,
        None,
    )
    .expect("fillet serves")
    .brep;

    let mut tree = FeatureTree::new();
    tree.set_cache_strategy(CacheStrategy::EveryStep);
    let anchor = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "anchor",
    );
    let r_old = tree
        .parameters_mut()
        .add_global("r", 0.1, Unit::Meters)
        .expect("param");
    let r_new = tree
        .parameters_mut()
        .add_global("r2", 0.2, Unit::Meters)
        .expect("param");
    let old_feature = Feature::Fillet {
        edges: vec![edge_name.clone()],
        radius: r_old,
        recipe: Default::default(),
    };
    let fillet_id = tree.add_feature(old_feature.clone(), "fillet");
    tree.cache
        .store_explicit(None, anchor, state_before.clone());
    tree.cache.store_explicit(None, fillet_id, state_at.clone());

    let g2 = Feature::Fillet {
        edges: vec![edge_name.clone()],
        radius: r_new,
        recipe: crate::feature_blend::FilletRecipeDesc {
            continuity: crate::feature_blend::ContinuityKind::G2,
            ..Default::default()
        },
    };
    assert!(
        tree.try_parameter_tweak_in_place(fillet_id, &old_feature, &g2)
            .is_none(),
        "continuity change must not take the radius-only in-place route"
    );

    // Same recipe on both sides, but a variable law: also not a
    // constant-radius tweak.
    let variable = crate::feature_blend::FilletRecipeDesc {
        radius_law: crate::feature_blend::FilletRadiusLaw::Linear { end: r_new },
        ..Default::default()
    };
    let old_var = Feature::Fillet {
        edges: vec![edge_name.clone()],
        radius: r_old,
        recipe: variable.clone(),
    };
    let new_var = Feature::Fillet {
        edges: vec![edge_name],
        radius: r_new,
        recipe: variable,
    };
    assert!(
        tree.try_parameter_tweak_in_place(fillet_id, &old_var, &new_var)
            .is_none(),
        "variable-law fillets must not take the constant-radius in-place route"
    );
}
