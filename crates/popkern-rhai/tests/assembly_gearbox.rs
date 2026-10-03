//! The gallery's planetary gearbox as an assembly (PK-630, PK-631,
//! PK-638): measured clearance between the meshing gears, a negative
//! interference control, placement by mates, and one assembly
//! parameter driving every gear's face width.

use popkern_assembly::{Assembly, Document};
use popkern_core::topology::BRep;
use popkern_rhai::ScriptHost;
use std::collections::BTreeMap;

fn root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).to_path_buf()
}

fn rewritten(stem: &str) -> Assembly {
    let source = std::fs::read_to_string(root().join(format!("examples/{stem}.rhai")))
        .expect("gallery script");
    let mut outcome = ScriptHost::new().run(&source, stem).expect("assembly runs");
    match outcome.take_document() {
        Some(Document::Assembly(assembly)) => *assembly,
        _ => panic!("{stem}: the gallery script authors an assembly"),
    }
}

#[test]
fn the_gearbox_meshes_with_measured_clearance_and_no_interference() {
    let mut gearbox = rewritten("planetary_gearbox");
    let report = gearbox
        .interference(&popkern_assembly::InterferenceOptions {
            clearance: 0.05,
            ..Default::default()
        })
        .expect("report");
    let mut sun_planet = Vec::new();
    for pair in &report.pairs {
        match &pair.status {
            popkern_assembly::PairStatus::Clear {
                distance,
                lower_bound,
                ..
            } => {
                assert!(*lower_bound > 0.0, "{} / {}: {lower_bound}", pair.a, pair.b);
                if pair.a.starts_with("sun") && pair.b.starts_with("planet") {
                    sun_planet.push(*distance);
                }
            }
            other => panic!("{} / {}: {other:?}", pair.a, pair.b),
        }
    }
    eprintln!("[pk630] gearbox pairs within 0.05: {:#?}", report.pairs);
    assert_eq!(
        sun_planet.len(),
        3,
        "each planet meshes the sun within 0.05"
    );
    // Backlash: every flank pair is apart, and the three planets see the
    // same clearance (they are one part placed by symmetry).
    let (lo, hi) = sun_planet
        .iter()
        .fold((f64::INFINITY, 0.0f64), |(l, h), d| (l.min(*d), h.max(*d)));
    assert!(lo > 0.0 && hi - lo < 1e-9, "{sun_planet:?}");
}

#[test]
fn a_planet_turned_into_the_sun_interferes() {
    let source =
        std::fs::read_to_string(root().join("examples/planetary_gearbox.rhai")).expect("script");
    // Half a tooth off its mesh phase: tooth on tooth.
    let clash = source.replace(
        ".rotate(z_axis, psi * 180.0 / pi)",
        ".rotate(z_axis, psi * 180.0 / pi + 20.0)",
    );
    assert_ne!(clash, source);
    let mut outcome = ScriptHost::new().run(&clash, "clash").expect("runs");
    let Some(Document::Assembly(mut gearbox)) = outcome.take_document() else {
        panic!("assembly")
    };
    let report = gearbox
        .interference(&popkern_assembly::InterferenceOptions::default())
        .expect("report");
    let clashes: Vec<_> = report
        .pairs
        .iter()
        .filter(|p| matches!(p.status, popkern_assembly::PairStatus::Interfering { .. }))
        .map(|p| (p.a.clone(), p.b.clone(), p.status.clone()))
        .collect();
    eprintln!("[pk630] clashes: {clashes:#?}");
    assert!(
        clashes
            .iter()
            .any(|(a, b, _)| a == "sun_1" && b == "planet_0"),
        "{clashes:?}"
    );
    for (_, _, status) in &clashes {
        if let popkern_assembly::PairStatus::Interfering {
            volume: Some(v), ..
        } = status
        {
            assert!(*v > 0.0);
        }
    }
}

#[test]
fn the_gearbox_placed_by_mates_is_the_gearbox_placed_by_transforms() {
    let source = std::fs::read_to_string(root().join("tests/fixtures/gearbox_by_mates.rhai"))
        .expect("fixture");
    let mut outcome = ScriptHost::new()
        .run(&source, "gearbox_by_mates")
        .expect("mates solve");
    assert_eq!(
        outcome.printed,
        [
            "planet_0 free motions: 0",
            "planet_1 free motions: 0",
            "planet_2 free motions: 0"
        ],
        "each planet fully placed"
    );
    let Some(Document::Assembly(mated)) = outcome.take_document() else {
        panic!("assembly")
    };
    let placed = rewritten("planetary_gearbox");
    let by_path = |asm: &Assembly| -> BTreeMap<String, popkern_assembly::Placement> {
        asm.leaves()
            .into_iter()
            .map(|r| (r.path.to_string(), r.world))
            .collect()
    };
    let (want, got) = (by_path(&placed), by_path(&mated));
    assert_eq!(
        want.keys().collect::<Vec<_>>(),
        got.keys().collect::<Vec<_>>()
    );
    let mut worst: f64 = 0.0;
    for (path, w) in &want {
        let g = &got[path];
        for r in 0..3 {
            for c in 0..4 {
                worst = worst.max((w.transform().matrix[r][c] - g.transform().matrix[r][c]).abs());
            }
        }
    }
    eprintln!("[pk631] mates vs transforms: worst placement entry difference {worst:e}");
    assert!(worst < 1e-9, "{worst:e}");
    assert_eq!(mated.mates().len(), 3 + 3 + 3 * 3);
}

#[test]
fn one_assembly_parameter_drives_every_gear_width() {
    // PK-638: setting the assembly's `width` equals setting each gear's.
    let mut driven = rewritten("planetary_gearbox");
    let mut by_hand = rewritten("planetary_gearbox");
    let mut doc = Document::Assembly(Box::new(std::mem::replace(&mut driven, Assembly::new(""))));
    doc.set_param("width", 0.3).expect("assembly parameter");
    let Document::Assembly(mut driven) = doc else {
        unreachable!()
    };
    for part in ["sun", "planet", "ring"] {
        let id = by_hand.definition_named(part).expect(part);
        let model = by_hand.part_model_mut(id).expect("model");
        model.set_param("width", 0.3).expect("set");
        model.replay().expect("replay");
    }
    let (a, b) = (driven.flatten().expect("a"), by_hand.flatten().expect("b"));
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(&b) {
        let zmax = |brep: &BRep| {
            brep.vertices
                .values()
                .map(|v| v.position.z)
                .fold(f64::MIN, f64::max)
        };
        if ["sun_1", "ring_1"].contains(&x.path.to_string().as_str())
            || x.path.to_string().starts_with("planet")
        {
            assert!(
                (zmax(&x.brep) - 0.3).abs() < 1e-12,
                "{}: {}",
                x.path,
                zmax(&x.brep)
            );
        }
        let mut va: Vec<[f64; 3]> = x
            .brep
            .vertices
            .values()
            .map(|v| [v.position.x, v.position.y, v.position.z])
            .collect();
        let mut vb: Vec<[f64; 3]> = y
            .brep
            .vertices
            .values()
            .map(|v| [v.position.x, v.position.y, v.position.z])
            .collect();
        va.sort_by(|p, q| p.partial_cmp(q).expect("finite"));
        vb.sort_by(|p, q| p.partial_cmp(q).expect("finite"));
        assert_eq!(va, vb, "{}", x.path);
    }
    // A driven part parameter refuses a direct edit.
    let mut doc = Document::Assembly(driven);
    assert!(doc.set_param("planet:width", 0.2).is_err());
    // The parameter and its drives survive .pkn.
    let Document::Assembly(driven) = doc else {
        unreachable!()
    };
    let bytes = driven.to_pkn_bytes().expect("save");
    let Document::Assembly(back) = Document::from_pkn_bytes(&bytes).expect("open") else {
        panic!()
    };
    assert_eq!(back.params(), driven.params());
    assert_eq!(back.drives().len(), 3);
}
