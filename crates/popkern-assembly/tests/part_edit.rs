//! PK-678: a part definition edited in place — every instance follows,
//! mates re-solve on the edited geometry, and a refused edit (no body,
//! mates that no longer resolve, a family variant, a sub-assembly)
//! leaves the assembly exactly as it was.

use std::collections::BTreeMap;

use popkern_assembly::{
    Assembly, AssemblyError, FamilyVariant, Mate, MateEntity, MateKind, MateRef, Parent, Placement,
};
use popkern_script::Plane;
use popkern_script::prelude::*;

fn block(height: f64) -> Model {
    let mut m = Model::new();
    m.profile(Plane::XY)
        .rect(1.0, 0.5)
        .sketch("block")
        .expect("sketch")
        .extrude(&mut m, End::Blind(height.into()))
        .expect("extrude");
    m
}

fn face(instance: &str, name: &str) -> MateRef {
    MateRef {
        instance: instance.into(),
        entity: MateEntity::Face { name: name.into() },
    }
}

/// `moving` hangs under `base`: its top cap flush on base's bottom.
fn hanging() -> Assembly {
    let mut asm = Assembly::new("two");
    let d = asm.add_part("block", block(0.2)).expect("part");
    asm.place(Parent::Root, d, Some("base"), Placement::identity())
        .expect("base");
    asm.place(
        Parent::Root,
        d,
        Some("moving"),
        Placement::translation(3.0, -1.0, 0.7),
    )
    .expect("moving");
    asm.add_mate(Mate {
        kind: MateKind::Coincident { aligned: false },
        a: face("moving", "extrude-1/Cap/0"),
        b: face("base", "extrude-1/Modified/0@0"),
        drive: None,
    })
    .expect("flush");
    asm
}

fn moving_z(asm: &Assembly) -> f64 {
    asm.leaves()
        .into_iter()
        .find(|r| r.path.to_string() == "moving")
        .expect("moving")
        .world
        .origin()
        .z
}

#[test]
fn a_taller_block_moves_what_hangs_under_it() {
    let mut asm = hanging();
    assert!((moving_z(&asm) + 0.2).abs() < 1e-12);
    let id = asm.definition_named("block").expect("block");
    asm.replace_part_model(id, block(0.5)).expect("edit");
    assert!(
        (moving_z(&asm) + 0.5).abs() < 1e-9,
        "re-solved on the taller part: {}",
        moving_z(&asm)
    );
    assert_eq!(asm.instances(Parent::Root).expect("root").len(), 2);
    assert_eq!(asm.check_history(), Ok(()));
}

#[test]
fn a_refused_edit_leaves_the_assembly_as_it_was() {
    let mut asm = hanging();
    let id = asm.definition_named("block").expect("block");
    let before = asm.tree();
    // No body.
    assert!(matches!(
        asm.replace_part_model(id, Model::new()),
        Err(AssemblyError::EmptyDefinition { .. })
    ));
    // The mated face no longer has its name: the mates cannot solve.
    let mut renamed = block(0.5);
    let extrude = renamed.features().into_iter().nth(1).expect("extrude");
    renamed.rename(&extrude, "slab").expect("rename");
    assert!(asm.replace_part_model(id, renamed).is_err());
    assert_eq!(asm.tree(), before, "nothing moved");
    let model = asm.part_model_mut(id).expect("block");
    assert_eq!(
        model.features()[1].label,
        "extrude-1",
        "the original model is back"
    );
}

#[test]
fn variants_and_sub_assemblies_are_not_edited_in_place() {
    let mut asm = Assembly::new("kinds");
    let variant = asm
        .add_family_variant(
            "block[h=0.2]",
            block(0.2),
            FamilyVariant {
                family: "block".into(),
                params: BTreeMap::new(),
            },
        )
        .expect("variant");
    let stage = asm.add_subassembly("stage").expect("stage");
    for id in [variant, stage] {
        assert!(matches!(
            asm.replace_part_model(id, block(0.5)),
            Err(AssemblyError::WrongDefinitionKind { .. })
        ));
    }
}
