//! PK-666 assembly construction history: one step per authored place,
//! move, mate or solve; the replayed placements end at the document's
//! own; mates show the move they caused; opened files say their log is
//! not an authoring record.

use popkern_assembly::history::AssemblyStepKind;
use popkern_assembly::{
    Assembly, Datum, Document, InstanceRef, Mate, MateEntity, MateKind, MateRef, Parent, Placement,
};
use popkern_script::Plane;
use popkern_script::prelude::*;

fn block() -> Model {
    let mut m = Model::new();
    m.profile(Plane::XY)
        .rect(1.0, 0.5)
        .sketch("block")
        .expect("sketch")
        .extrude(&mut m, End::Blind(0.2.into()))
        .expect("extrude");
    m
}

fn datum(instance: &str, datum: Datum) -> MateRef {
    MateRef {
        instance: instance.into(),
        entity: MateEntity::Datum { datum },
    }
}

fn root(index: usize) -> InstanceRef {
    InstanceRef {
        parent: Parent::Root,
        index,
    }
}

fn at(asm: &Assembly, count: usize, index: usize) -> Option<Placement> {
    asm.placements_after(count)
        .into_iter()
        .find(|(r, _)| *r == root(index))
        .map(|(_, p)| p)
}

/// base, then moving placed and translated (one step), then a flush
/// mate that pulls moving down onto base's plane.
fn authored() -> Assembly {
    let mut asm = Assembly::new("two");
    let d = asm.add_part("block", block()).expect("part");
    asm.place(Parent::Root, d, Some("base"), Placement::identity())
        .expect("base");
    let moving = asm
        .place(Parent::Root, d, Some("moving"), Placement::identity())
        .expect("moving");
    asm.set_placement(moving, Placement::translation(3.0, -1.0, 0.7))
        .expect("translate");
    asm.add_mate(Mate {
        kind: MateKind::Coincident { aligned: true },
        a: datum("moving", Datum::Xy),
        b: datum("base", Datum::Xy),
    })
    .expect("flush");
    asm
}

#[test]
fn a_place_and_its_translate_are_one_step_and_a_mate_shows_its_move() {
    let asm = authored();
    let steps = &asm.history().steps;
    let kinds: Vec<&AssemblyStepKind> = steps.iter().map(|s| &s.kind).collect();
    assert_eq!(
        kinds,
        [
            &AssemblyStepKind::Place,
            &AssemblyStepKind::Place,
            &AssemblyStepKind::Mate(0)
        ]
    );
    assert_eq!(steps[1].label, "place moving (block)");
    assert!(steps[2].label.starts_with("mate: "), "{}", steps[2].label);
    assert!(asm.history().recorded);

    // Before the mate: moving sits where it was authored, at z = 0.7.
    assert_eq!(at(&asm, 0, 0), None, "nothing placed before step 1");
    assert_eq!(at(&asm, 1, 1), None, "moving not placed yet");
    let authored = at(&asm, 2, 1).expect("moving after step 2");
    assert_eq!(authored, Placement::translation(3.0, -1.0, 0.7));
    // After it: moved, and only moving is in the mate step's record.
    let solved = at(&asm, 3, 1).expect("moving after the mate");
    assert_ne!(solved, authored, "the mate moved it");
    assert_eq!(steps[2].placements.len(), 1, "base is the ground");
    assert_eq!(asm.check_history(), Ok(()));
    // The viewer's tree at a step: the occurrences placed by then.
    assert_eq!(asm.tree_after(1).len(), 1);
    let before_mate = asm.tree_after(2);
    assert_eq!(before_mate[1].world, authored);
    assert_eq!(
        asm.tree_after(3),
        asm.tree(),
        "the last step is the document"
    );
}

#[test]
fn a_refused_mate_and_a_no_op_solve_record_nothing() {
    let mut asm = authored();
    let before = asm.history().steps.len();
    let refused = asm.add_mate(Mate {
        kind: MateKind::Distance {
            distance: 0.5,
            aligned: true,
        },
        a: datum("base", Datum::Xy),
        b: datum("moving", Datum::Xy),
    });
    assert!(refused.is_err());
    asm.solve_mates().expect("solve");
    assert_eq!(asm.history().steps.len(), before);
    assert_eq!(asm.check_history(), Ok(()));
}

#[test]
fn consecutive_re_solves_merge_into_one_step() {
    let mut asm = authored();
    let moving = root(1);
    // Knock moving off its mate twice, re-solving each time: the moves
    // are steps, each solve is one step (the second has nothing new to
    // merge into because a move separates them).
    asm.set_placement(moving, Placement::translation(0.0, 0.0, 5.0))
        .expect("knock");
    asm.solve_mates().expect("solve");
    let n = asm.history().steps.len();
    asm.solve_mates().expect("solve again");
    assert_eq!(asm.history().steps.len(), n, "nothing moved: no step");
    assert_eq!(
        asm.history().steps[n - 1].kind,
        AssemblyStepKind::Solve,
        "the knock was solved back"
    );
    assert_eq!(asm.check_history(), Ok(()));
}

#[test]
fn an_opened_assembly_lists_its_instances_and_says_it_was_not_recorded() {
    let asm = authored();
    let bytes = asm.to_pkn_bytes().expect("save");
    let Document::Assembly(opened) = Document::from_pkn_bytes(&bytes).expect("open") else {
        panic!("not an assembly");
    };
    assert!(!opened.history().recorded);
    assert_eq!(opened.history().steps.len(), 2, "one step per instance");
    assert_eq!(opened.check_history(), Ok(()));
    assert_eq!(
        at(&opened, 2, 1),
        at(&asm, 3, 1),
        "at the stored (solved) placement"
    );
}

#[test]
fn sub_assembly_placements_are_steps_and_replay_exactly() {
    let mut asm = Assembly::new("nested");
    let d = asm.add_part("block", block()).expect("part");
    let stage = asm.add_subassembly("stage").expect("stage");
    asm.place(
        Parent::Definition(stage),
        d,
        Some("a"),
        Placement::translation(1.0, 0.0, 0.0),
    )
    .expect("a");
    asm.place(Parent::Root, stage, Some("s"), Placement::identity())
        .expect("s");
    let labels: Vec<&str> = asm
        .history()
        .steps
        .iter()
        .map(|s| s.label.as_str())
        .collect();
    assert_eq!(labels, ["place stage/a (block)", "place s (stage)"]);
    assert_eq!(asm.check_history(), Ok(()));
    assert!(
        asm.tree_after(1).is_empty(),
        "stage holds a, but s is not placed yet"
    );
    let rows = asm.tree_after(2);
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[1].path.to_string(), "s/a");
}
