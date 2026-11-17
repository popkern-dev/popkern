//! The OCC oracle (design §6.2 oracle 9): a served result and operand A
//! exported to STEP, checked by OpenCASCADE (reader, BRepCheck, free
//! bounds) and measured by its adaptive volume. Ratios to A's volume are
//! compared, so the units of the export do not enter.

use std::path::{Path, PathBuf};

use popkern_core::topology::BRep;

/// The OCP interpreter, when installed.
pub fn python() -> Option<PathBuf> {
    let p = Path::new(&std::env::var("HOME").ok()?).join("venvs/ocp/bin/python");
    p.exists().then_some(p)
}

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn export(brep: &BRep, path: &Path) -> Result<(), String> {
    let name = path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "part".into());
    let doc = popkern_io::StepDocument {
        products: vec![popkern_io::StepProductNode {
            name: name.clone(),
            id: name.to_uppercase(),
            children: vec![],
            shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                kind: popkern_io::ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: popkern_io::RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: popkern_io::ProductMetadata::default(),
        }],
        roots: vec![0],
        units: popkern_io::UnitContext::default(),
        opaque_bag: popkern_io::OpaqueBag::default(),
    };
    popkern_io::export_step_file(&doc, path, &popkern_io::StepExportOptions::default())
        .map_err(|e| format!("export: {e:?}"))
}

/// OCC's verdict on a body: `Ok(volume)` when it reads and checks clean.
pub fn measure(brep: &BRep, tag: &str) -> Result<f64, String> {
    let py = python().ok_or("OCC oracle absent")?;
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("corpus-occ");
    std::fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
    let path = dir.join(format!("{tag}.step"));
    export(brep, &path)?;
    let check = std::process::Command::new(&py)
        .arg(root().join("tools/occ-oracle/occ_check.py"))
        .arg(&path)
        .output()
        .map_err(|e| e.to_string())?;
    let report = String::from_utf8_lossy(&check.stdout);
    if !report.contains("\"passed\": true") {
        return Err(format!(
            "OCC check rejected the STEP: {}",
            report.chars().take(300).collect::<String>()
        ));
    }
    let solids = std::process::Command::new(&py)
        .arg(root().join("tools/occ-oracle/occ_solids.py"))
        .arg(&path)
        .arg("--adaptive")
        .output()
        .map_err(|e| e.to_string())?;
    let v: serde_json::Value =
        serde_json::from_slice(&solids.stdout).map_err(|e| format!("occ_solids output: {e}"))?;
    let total: f64 = v["solids"]
        .as_array()
        .ok_or("no solids")?
        .iter()
        .filter_map(|s| s["volume"].as_f64())
        .sum();
    let _ = std::fs::remove_file(&path);
    Ok(total)
}
