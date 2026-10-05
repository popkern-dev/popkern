// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
//
// PK-677: references rebind by provenance — a verified exact name is
// kept, a missing or repurposed name is rewritten to the entity its
// snapshot describes, and an ambiguous snapshot declines.

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, NameId, OperationId, RoleCategory};
use popkern_core::naming_provenance::{EntityKind, NameValue};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{ReferenceProvenance, rebind};
use crate::FeatureTree;
use crate::error::ReplayError;
use crate::feature::{ExtrudeEnd, Feature, OperationMode, SketchPlaneDesc, Vec3Desc};
use crate::replay_report::RecordSource;
use crate::unit::Unit;

fn name(op: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

/// A 1×1×1 box: the body after sketch + extrude.
fn box_brep() -> BRep {
    let mut sk = Sketch::new();
    let p: Vec<_> = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        .iter()
        .map(|&(x, y)| sk.add_point(x, y))
        .collect();
    for i in 0..4 {
        sk.add_line(p[i], p[(i + 1) % 4]);
        sk.add_constraint(Constraint::Fixed(p[i]));
    }
    sk.solve().expect("solve");
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: crate::SketchPayload::Constraint(
                serde_json::to_string(&sk).expect("json"),
            ),
        },
        "Sketch 1",
    );
    let dist = tree
        .parameters_mut()
        .add_global("dist", 1.0, Unit::Meters)
        .expect("param");
    tree.add_feature(
        Feature::Extrude {
            face: name(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude 1",
    );
    tree.brep_at(2, Tolerance::default()).expect("box")
}

/// A feature whose one reference is `reference` (an extrude's face
/// field stands in for any reference slot).
fn referencing(reference: EntityName) -> Feature {
    Feature::Extrude {
        face: reference,
        direction: None,
        end_condition: ExtrudeEnd::ThroughAll,
        mode: OperationMode::NewBody,
        draft_angle: None,
        second_direction: None,
    }
}

fn face_of(feature: &Feature) -> &EntityName {
    match feature {
        Feature::Extrude { face, .. } => face,
        _ => unreachable!("the test feature is an extrude"),
    }
}

/// Two named edges of the box and their name ids.
fn two_edges(brep: &BRep) -> [(EntityName, NameId); 2] {
    let mut named = brep
        .edges
        .values()
        .filter_map(|e| Some((brep.name_pool.get(e.name?)?.clone(), e.name?)));
    let a = named.next().expect("edge");
    let b = named.find(|(n, _)| *n != a.0).expect("another edge");
    [a, b]
}

fn snapshot(name: &EntityName, value: NameValue) -> ReferenceProvenance {
    ReferenceProvenance {
        name: name.clone(),
        kind: EntityKind::Edge,
        value,
    }
}

#[test]
fn a_name_that_still_names_its_entity_is_kept() {
    let brep = box_brep();
    let [(a, id), _] = two_edges(&brep);
    let value = brep.name_pool.value_of(id).expect("value");
    let (feature, fallbacks) =
        rebind(&referencing(a.clone()), &[snapshot(&a, value)], &brep).expect("rebind");
    assert_eq!(face_of(&feature), &a);
    assert!(fallbacks.is_empty());
}

#[test]
fn a_missing_name_rebinds_to_the_entity_its_snapshot_describes() {
    let brep = box_brep();
    let [(a, id), _] = two_edges(&brep);
    let gone = name(999, RoleCategory::Modified, 7);
    let value = brep.name_pool.value_of(id).expect("value");
    let (feature, fallbacks) =
        rebind(&referencing(gone.clone()), &[snapshot(&gone, value)], &brep).expect("rebind");
    assert_eq!(face_of(&feature), &a);
    assert_eq!(fallbacks.len(), 1);
    assert_eq!(fallbacks[0].source, RecordSource::Naming);
}

/// The hazard: after an upstream edit, the recorded numbers name a
/// different entity. The snapshot disagrees, so the reference moves
/// to the entity it was recorded on instead of silently binding.
#[test]
fn a_name_that_now_denotes_another_entity_follows_its_snapshot() {
    let brep = box_brep();
    let [(a, _), (b, id_b)] = two_edges(&brep);
    let recorded_on_b = brep.name_pool.value_of(id_b).expect("value");
    let (feature, fallbacks) = rebind(
        &referencing(a.clone()),
        &[snapshot(&a, recorded_on_b)],
        &brep,
    )
    .expect("rebind");
    assert_eq!(face_of(&feature), &b);
    assert_eq!(fallbacks.len(), 1);
}

#[test]
fn a_snapshot_several_entities_match_declines() {
    let brep = box_brep();
    let [(_, id), _] = two_edges(&brep);
    let mut value = brep.name_pool.value_of(id).expect("value");
    // Every index renumbered: only operations and categories remain,
    // which the box's sibling edges share.
    fn renumber(v: &mut NameValue) {
        v.role.index += 1000;
        v.ancestors.iter_mut().for_each(renumber);
    }
    renumber(&mut value);
    let gone = name(999, RoleCategory::Modified, 7);
    match rebind(&referencing(gone.clone()), &[snapshot(&gone, value)], &brep) {
        Err(ReplayError::AmbiguousReference { candidates, .. }) => assert!(candidates > 1),
        other => panic!("expected an ambiguity decline, got {other:?}"),
    }
}
