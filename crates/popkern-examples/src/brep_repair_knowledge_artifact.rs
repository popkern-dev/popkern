//! Lightweight BRep repair knowledge JSON artifact example.

use std::path::{Path, PathBuf};

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BrepRepairDiagnosticExport};
use popkern_examples::repair_artifacts::write_whole_brep_repair_knowledge_artifact;

fn main() {
    let output_dir = output_dir();
    std::fs::create_dir_all(&output_dir).expect("create repair knowledge output dir");
    let brep = BRep::unit_cube(Tolerance::default()).expect("unit cube");
    let path = output_dir.join("unit_cube_repair_knowledge.json");
    write_whole_brep_repair_knowledge_artifact(
        &path,
        "unit_cube_non_filleting",
        &brep,
        5,
        BrepRepairDiagnosticExport::CompactJson,
    )
    .expect("write repair knowledge artifact");
    println!("repair knowledge {}", path.display());
}

fn output_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("examples")
        .join("output")
        .join("brep_repair_knowledge")
}
