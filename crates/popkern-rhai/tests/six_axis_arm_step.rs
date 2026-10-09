//! The six-axis arm's STEP product structure as an independent reader
//! sees it (design §9.7, PK-725): OpenCASCADE's XCAF reader finds every
//! part occurrence at its path, naming its part, at its placed bounding
//! box, in its colour, with the part's volume.
//!
//! Runs when the OCP venv (`~/venvs/ocp`) is installed; otherwise it
//! says so and checks nothing.

mod arm_support;

use std::path::{Path, PathBuf};

use arm_support::*;
use serde_json::Value;

/// The OCC assembly oracle's report on `step`, when the venv exists.
fn occ_assembly(step: &Path) -> Option<Value> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let python =
        PathBuf::from(std::env::var("HOME").unwrap_or_default()).join("venvs/ocp/bin/python");
    let script = root.join("tools/occ-oracle/occ_assembly.py");
    if !python.exists() || !script.exists() {
        eprintln!("[pk725] OCC oracle absent; the arm's STEP is not cross-checked");
        return None;
    }
    let output = std::process::Command::new(python)
        .arg(script)
        .arg(step)
        .output()
        .expect("oracle");
    Some(serde_json::from_slice(&output.stdout).expect("oracle json"))
}

#[test]
fn an_independent_reader_finds_every_occurrence_where_and_as_we_place_it() {
    let mut asm = build();
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("pk725");
    std::fs::create_dir_all(&dir).expect("dir");
    let step = dir.join("six_axis_arm.step");
    asm.to_step(&step).expect("STEP export");
    let Some(report) = occ_assembly(&step) else {
        return;
    };
    assert_eq!(
        report["check_messages"],
        serde_json::json!([]),
        "reader messages"
    );
    let instances = report["instances"].as_array().expect("instances");
    let rows: Vec<_> = asm.leaves().into_iter().filter(|r| r.is_part).collect();
    assert_eq!(instances.len(), rows.len(), "occurrence count");

    // Each definition body once: its mesh (for placed bounds) and volume.
    let parts = asm.definition_parts().expect("parts");
    let mut ours = std::collections::HashMap::new();
    for (id, bodies) in &parts {
        let mut meshes = Vec::new();
        let mut volume = 0.0;
        for part in bodies {
            let brep = part.brep();
            meshes.push(
                popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
                    .expect("mesh"),
            );
            volume += popkern_measure::body_volume(brep, &brep.tolerance).expect("volume");
        }
        ours.insert(*id, (meshes, volume));
    }

    let mut worst_volume: (f64, String) = (0.0, String::new());
    let mut occ_split = std::collections::BTreeSet::new();
    for (occ, row) in instances.iter().zip(&rows) {
        let path = row.path.to_string();
        assert_eq!(occ["path"], path.as_str(), "occurrence order");
        let definition = asm.definition(row.definition).expect("definition");
        assert_eq!(occ["part"], definition.name.as_str(), "{path}: part name");
        // Placed bounds: our mesh vertices lie on the surfaces, so their
        // box sits inside the true box by at most the chord sagitta.
        let (meshes, volume) = &ours[&row.definition];
        let m = row.world.transform().matrix;
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for mesh in meshes {
            for v in &mesh.vertices {
                let p = v.position;
                for i in 0..3 {
                    let c = m[i][0] * p.x + m[i][1] * p.y + m[i][2] * p.z + m[i][3];
                    lo[i] = lo[i].min(c);
                    hi[i] = hi[i].max(c);
                }
            }
        }
        for i in 0..3 {
            for (k, ours) in [(0, lo[i]), (1, hi[i])] {
                let theirs = occ["bbox"][k][i].as_f64().expect("bbox");
                assert!(
                    (theirs - ours).abs() < 0.1,
                    "{path}: bbox[{k}][{i}] OCC {theirs} vs ours {ours}"
                );
            }
        }
        // Colour: the definition's appearance, as written.
        let colour = definition
            .model()
            .and_then(|model| model.appearance())
            .map(|a| a.rgb_f32());
        match colour {
            Some(rgb) => {
                for (c, theirs) in rgb.iter().zip(occ["color"].as_array().expect("colour")) {
                    let theirs = theirs.as_f64().expect("rgb");
                    assert!(
                        (theirs - *c as f64).abs() < 1e-3,
                        "{path}: colour {theirs} vs {c}"
                    );
                }
            }
            None => assert!(
                occ["color"].is_null(),
                "{path}: an uncoloured part read a colour"
            ),
        }
        // Volume against OCC's own tessellation of the STEP solid, and
        // against its exact integration wherever OCC agrees with itself:
        // BRepGProp diverges on multi-turn NURBS threads (the screws read
        // 0.5% low adaptively, the cycloid discs 0.34% low by default).
        let meshed = occ["mesh_volume"].as_f64().expect("mesh volume");
        let exact = occ["volume"].as_f64().expect("volume");
        let mut against = vec![meshed];
        if ((exact - meshed) / meshed).abs() < 1e-3 {
            against.push(exact);
        } else {
            occ_split.insert(definition.name.clone());
        }
        for theirs in against {
            let rel = ((volume - theirs) / theirs).abs();
            if rel > worst_volume.0 {
                worst_volume = (
                    rel,
                    format!("{} ours {volume} OCC {theirs}", definition.name),
                );
            }
        }
    }
    eprintln!(
        "worst volume disagreement: {:.2e} ({})",
        worst_volume.0, worst_volume.1
    );
    eprintln!("parts whose OCC integration and tessellation disagree: {occ_split:?}");
    assert!(worst_volume.0 < 1e-3, "volume: {}", worst_volume.1);
}
