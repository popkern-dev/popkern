//! PK-631 mate solver: full placement, free motions, redundancy,
//! conflict (transactional), signed angles, face mates.

use popkern_assembly::{
    Assembly, AssemblyError, Datum, Mate, MateEntity, MateKind, MateRef, Parent, Placement,
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

fn two() -> Assembly {
    let mut asm = Assembly::new("two");
    let d = asm.add_part("block", block()).expect("part");
    asm.place(Parent::Root, d, Some("base"), Placement::identity())
        .expect("base");
    asm.place(
        Parent::Root,
        d,
        Some("moving"),
        Placement::translation(3.0, -1.0, 0.7),
    )
    .expect("moving");
    asm
}

fn datum(instance: &str, datum: Datum) -> MateRef {
    MateRef {
        instance: instance.into(),
        entity: MateEntity::Datum { datum },
    }
}

fn mate(kind: MateKind, a: MateRef, b: MateRef) -> Mate {
    Mate { kind, a, b }
}

fn moving(asm: &Assembly) -> Placement {
    asm.leaves()
        .into_iter()
        .find(|r| r.path.to_string() == "moving")
        .expect("moving")
        .world
}

#[test]
fn three_mates_place_an_instance_fully_and_report_free_motions_on_the_way() {
    let mut asm = two();
    let r = asm
        .add_mate(mate(
            MateKind::Concentric,
            datum("moving", Datum::Z),
            datum("base", Datum::Z),
        ))
        .expect("concentric");
    assert_eq!(
        r.free_dof["moving"], 2,
        "slide along and turn about the axis remain"
    );
    let r = asm
        .add_mate(mate(
            MateKind::Coincident { aligned: true },
            datum("moving", Datum::Xy),
            datum("base", Datum::Xy),
        ))
        .expect("flush");
    assert_eq!(r.free_dof["moving"], 1, "the turn remains");
    assert_eq!(
        r.redundant,
        Vec::<usize>::new()
            .into_iter()
            .chain([1])
            .collect::<Vec<_>>(),
        "the flush mate's normal rows repeat the axis rows"
    );
    let about = datum("base", Datum::Z);
    let r = asm
        .add_mate(mate(
            MateKind::Angle {
                degrees: 30.0,
                about: Some(about),
            },
            datum("base", Datum::Zx),
            datum("moving", Datum::Zx),
        ))
        .expect("angle");
    assert_eq!(r.free_dof["moving"], 0);
    assert!(r.converged && r.residual <= 1e-10);
    let want = Placement::rotation(
        popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
        popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        30.0,
    )
    .expect("turn");
    assert!(
        moving(&asm).transform().approx_eq(want.transform(), 1e-12),
        "{:?}",
        moving(&asm)
    );
    assert_eq!(asm.mates().len(), 3);
}

#[test]
fn a_conflicting_mate_is_refused_and_nothing_moves() {
    let mut asm = two();
    asm.add_mate(mate(
        MateKind::Coincident { aligned: true },
        datum("moving", Datum::Xy),
        datum("base", Datum::Xy),
    ))
    .expect("flush");
    let before = moving(&asm);
    match asm.add_mate(mate(
        MateKind::Distance {
            distance: 0.5,
            aligned: true,
        },
        datum("base", Datum::Xy),
        datum("moving", Datum::Xy),
    )) {
        Err(AssemblyError::MateConflict { report }) => {
            assert!(!report.converged);
            assert_eq!(
                report.conflicting.len(),
                2,
                "both mates on the one height: {report:?}"
            );
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(moving(&asm), before, "transactional");
    assert_eq!(asm.mates().len(), 1, "the conflicting mate was removed");
}

#[test]
fn a_face_mate_uses_the_part_face_and_a_wrong_pairing_is_refused() {
    let mut asm = two();
    // The top cap of "moving" (extrude-1/Cap/0) flush against the
    // bottom face of "base" (the sketch face the extrude modified):
    // the moving block sits under it (opposed normals, distance 0).
    let top = MateRef {
        instance: "moving".into(),
        entity: MateEntity::Face {
            name: "extrude-1/Cap/0".into(),
        },
    };
    let bottom = MateRef {
        instance: "base".into(),
        entity: MateEntity::Face {
            name: "extrude-1/Modified/0@0".into(),
        },
    };
    let r = asm.add_mate(mate(
        MateKind::Coincident { aligned: false },
        top.clone(),
        bottom,
    ));
    match r {
        Ok(report) => {
            assert!(report.converged);
            let z = moving(&asm).origin().z;
            assert!(
                (z + 0.2).abs() < 1e-12,
                "the moving block's top at z = 0: origin z {z}"
            );
        }
        Err(e) => panic!("{e}"),
    }
    // A plane cannot be concentric with anything.
    assert!(matches!(
        asm.add_mate(mate(MateKind::Concentric, top, datum("base", Datum::Z))),
        Err(AssemblyError::Mate { .. })
    ));
}

#[test]
fn mates_survive_save_and_open() {
    let mut asm = two();
    asm.add_mate(mate(
        MateKind::Concentric,
        datum("moving", Datum::Z),
        datum("base", Datum::Z),
    ))
    .expect("concentric");
    let bytes = asm.to_pkn_bytes().expect("save");
    let popkern_assembly::Document::Assembly(mut back) =
        popkern_assembly::Document::from_pkn_bytes(&bytes).expect("open")
    else {
        panic!("assembly")
    };
    assert_eq!(back.mates(), asm.mates());
    assert_eq!(moving(&back), moving(&asm));
    back.solve_mates().expect("re-solve holds");
    assert!(
        moving(&back)
            .transform()
            .approx_eq(moving(&asm).transform(), 1e-12)
    );
}
