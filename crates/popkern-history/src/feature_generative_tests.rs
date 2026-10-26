use super::*;
use crate::feature::Feature;
use popkern_core::naming::{EntityRole, OperationId, RoleCategory};

fn name(op: u64, i: u32) -> EntityName {
    EntityName {
        operation: OperationId(op),
        role: EntityRole {
            category: RoleCategory::Generated,
            index: i,
        },
        ancestors: Vec::new(),
    }
}

/// Four parameter ids (volume fraction, filter, cell, force).
fn ids() -> [ParameterId; 4] {
    let mut keys = slotmap::SlotMap::<ParameterId, ()>::with_key();
    [(); 4].map(|_| keys.insert(()))
}

fn desc() -> GenerativeDesc {
    let [vf, filter, cell, force] = ids();
    GenerativeDesc {
        material: "PLA (FDM)".into(),
        supports: vec![GenerativeSupportDesc {
            faces: vec![name(1, 0)],
            kind: GenerativeSupportKindDesc::Fixed,
        }],
        loads: vec![GenerativeLoadDesc {
            faces: vec![name(1, 3)],
            force,
            direction: [0.0, 0.0, -1.0],
        }],
        keep: vec![vec![name(1, 0)], vec![name(1, 3)]],
        volume_fraction: vf,
        filter,
        keep_depth: None,
        cell: Some(cell),
        build_direction: Some(BuildDirectionDesc {
            axis: 1,
            negative: false,
        }),
        max_iterations: 300,
        robust: false,
        max_cells: None,
        mm_per_unit: 1.0,
    }
}

#[test]
fn the_recipe_round_trips_through_serde() {
    let feature = Feature::Generative {
        spec: Box::new(desc()),
    };
    let text = serde_json::to_string(&feature).unwrap_or_else(|e| panic!("{e}"));
    let back: Feature = serde_json::from_str(&text).unwrap_or_else(|e| panic!("{e}"));
    match back {
        Feature::Generative { spec } => assert_eq!(*spec, desc()),
        other => panic!("{other:?}"),
    }
}

#[test]
fn parameters_and_names_are_listed_in_a_fixed_order() {
    let d = desc();
    let [vf, filter, cell, force] = ids();
    assert_eq!(
        d.param_ids(),
        vec![
            ("volume_fraction", vf),
            ("filter", filter),
            ("cell", cell),
            ("force", force),
        ]
    );
    // Supports, loads, then kept faces; a name used twice is listed
    // twice (each reference rebinds on its own).
    assert_eq!(
        d.entity_names(),
        vec![&name(1, 0), &name(1, 3), &name(1, 0), &name(1, 3)]
    );
    let mut m = d.clone();
    for n in m.entity_names_mut() {
        n.operation = OperationId(9);
    }
    assert!(
        m.entity_names()
            .iter()
            .all(|n| n.operation == OperationId(9))
    );
}

#[test]
fn build_directions_spell_as_the_skin_takes_them() {
    let s = |axis, negative| BuildDirectionDesc { axis, negative }.spelling();
    assert_eq!([s(0, false), s(1, true), s(2, false)], ["+x", "-y", "+z"]);
}
