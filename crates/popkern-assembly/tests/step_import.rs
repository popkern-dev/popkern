//! PK-631: STEP assembly import — a third-party writer's product
//! structure (OpenCASCADE XCAF, `tools/occ-oracle/occ_write_assembly.py`)
//! and popkern's own exports come back as an assembly with the same
//! tree, names and world placements.

use popkern_assembly::{Assembly, Parent, Placement};
use popkern_core::topology::BRep;
use popkern_script::Plane;
use popkern_script::prelude::*;

fn bounds(brep: &BRep) -> [[f64; 3]; 2] {
    let mut b = [[f64::MAX; 3], [f64::MIN; 3]];
    for v in brep.vertices.values() {
        for (k, x) in [v.position.x, v.position.y, v.position.z]
            .into_iter()
            .enumerate()
        {
            b[0][k] = b[0][k].min(x);
            b[1][k] = b[1][k].max(x);
        }
    }
    b
}

#[test]
fn an_opencascade_assembly_imports_with_its_tree_names_and_placements() {
    let bytes = include_bytes!("fixtures/occ_bench_assembly.step");
    let (mut asm, report) = Assembly::from_step_bytes(bytes, "occ_bench_assembly").expect("import");
    assert!(
        report.refused_solids.is_empty() && report.refused_instances.is_empty(),
        "{report:?}"
    );
    assert_eq!(asm.name(), "bench");
    let rows: Vec<(String, String)> = asm
        .tree()
        .into_iter()
        .map(|r| {
            (
                r.path.to_string(),
                asm.definition(r.definition).expect("def").name.clone(),
            )
        })
        .collect();
    assert_eq!(
        rows,
        [
            ("top".to_string(), "plate".to_string()),
            ("leg_front".into(), "leg".into()),
            ("leg_back".into(), "leg".into()),
            ("pair".into(), "leg_pair".into()),
            ("pair/left".into(), "leg".into()),
            ("pair/right".into(), "leg".into()),
        ]
    );
    assert_eq!(report.parts, 2, "the leg is ONE part, placed four times");
    // World positions as OpenCASCADE reads them (occ_assembly.py), to
    // its bounding-box precision.
    let want = [
        ("top", [-100.0, -40.0, 60.0]),
        ("leg_front", [77.0, 22.0, 0.0]),
        ("leg_back", [-93.0, 22.0, 0.0]),
        ("pair/left", [77.0, -38.0, 0.0]),
        ("pair/right", [-93.0, -38.0, 0.0]),
    ];
    let flat = asm.flatten().expect("flatten");
    for (path, lo) in want {
        let body = flat
            .iter()
            .find(|b| b.path.to_string() == path)
            .expect(path);
        let b = bounds(&body.brep);
        for k in 0..3 {
            // Vertex bounds of a cylinder sit inside its true box by
            // at most the radius on the seam-less sides.
            assert!(
                b[0][k] >= lo[k] - 1e-6 && b[0][k] <= lo[k] + 16.0 + 1e-6,
                "{path} {k}: {} vs {}",
                b[0][k],
                lo[k]
            );
        }
    }
    let top = flat
        .iter()
        .find(|b| b.path.to_string() == "top")
        .expect("top");
    assert_eq!(
        bounds(&top.brep),
        [[-100.0, -40.0, 60.0], [100.0, 40.0, 70.0]]
    );
}

#[test]
fn a_popkern_export_imports_back_to_the_same_tree() {
    let mut block = Model::new();
    block
        .profile(Plane::XY)
        .rect(1.0, 0.5)
        .sketch("block")
        .expect("sketch")
        .extrude(&mut block, End::Blind(0.2.into()))
        .expect("extrude");
    let mut asm = Assembly::new("stack");
    let part = asm.add_part("block", block).expect("part");
    let pair = asm.add_subassembly("pair").expect("pair");
    asm.place(
        Parent::Definition(pair),
        part,
        Some("low"),
        Placement::identity(),
    )
    .expect("low");
    asm.place(
        Parent::Definition(pair),
        part,
        Some("high"),
        Placement::translation(0.0, 0.0, 0.3),
    )
    .expect("high");
    asm.place(Parent::Root, part, Some("base"), Placement::identity())
        .expect("base");
    let turned = Placement::rotation(
        popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
        popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        45.0,
    )
    .expect("turn")
    .then(&Placement::translation(2.0, 0.0, 0.0));
    asm.place(Parent::Root, pair, Some("tower"), turned)
        .expect("tower");
    let bytes = asm.to_step_bytes("stack").expect("export");
    let (back, report) = Assembly::from_step_bytes(&bytes, "stack").expect("import");
    assert!(report.refused_instances.is_empty(), "{report:?}");
    let (a, b) = (asm.tree(), back.tree());
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(&b) {
        assert_eq!(x.path, y.path);
        assert!(
            x.world.transform().approx_eq(y.world.transform(), 1e-9),
            "{}",
            x.path
        );
    }
}
