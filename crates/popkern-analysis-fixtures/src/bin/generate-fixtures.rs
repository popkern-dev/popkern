//! Writes every fixture's STEP file into `fixtures/`. Run only on a
//! popkern pin bump, and record the diff in `PINNED.md`.

use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/analysis");
    std::fs::create_dir_all(&dir)?;
    for fixture in &popkern_analysis_fixtures::FIXTURES {
        let bytes = popkern_analysis_fixtures::step_bytes(fixture)?;
        let path = dir.join(format!("{}.step", fixture.name));
        std::fs::write(&path, &bytes)?;
        eprintln!("wrote {} ({} bytes)", path.display(), bytes.len());
    }
    Ok(())
}
