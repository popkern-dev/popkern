//! STEP assembly import corpus measurement (PK-631, design §8.3).
//!
//!     cargo run --release -p popkern-assembly --example step_corpus -- <dir> [report.md]
//!
//! For every `.stp`/`.step` file in `<dir>`: import it as an assembly
//! and record what loaded — products, occurrences placed against the
//! file's `NEXT_ASSEMBLY_USAGE_OCCURRENCE` count, solids recorded
//! against the file's `MANIFOLD_SOLID_BREP` count, refused solids and
//! occurrences, part occurrences against OpenCASCADE's (XCAF, when the
//! OCP oracle is installed), and how many part bodies tessellate
//! watertight. Prints a markdown table (and writes it to the report
//! path when given).

use std::path::{Path, PathBuf};
use std::process::Command;

use popkern_assembly::Assembly;

struct Row {
    file: String,
    ok: Result<(), String>,
    nauo: usize,
    msb: usize,
    products: usize,
    parts: usize,
    occurrences: usize,
    occ_occurrences: Option<usize>,
    solids: usize,
    refused_solids: usize,
    refused_instances: usize,
    bodies: usize,
    watertight: usize,
    ms: u128,
}

fn count(text: &str, entity: &str) -> usize {
    text.matches(&format!("{entity}(")).count()
}

fn occ_occurrences(path: &Path) -> Option<usize> {
    let python = PathBuf::from(std::env::var("HOME").ok()?).join("venvs/ocp/bin/python");
    let script =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle/occ_assembly.py");
    if !python.exists() {
        return None;
    }
    let out = Command::new(python).arg(script).arg(path).output().ok()?;
    let report: serde_json::Value = serde_json::from_slice(&out.stdout).ok()?;
    report["instances"].as_array().map(Vec::len)
}

fn measure(path: &Path) -> Row {
    let file = path
        .file_name()
        .map(|f| f.to_string_lossy().into_owned())
        .unwrap_or_default();
    let bytes = std::fs::read(path).unwrap_or_default();
    let text = String::from_utf8_lossy(&bytes);
    let mut row = Row {
        file: file.clone(),
        ok: Ok(()),
        nauo: count(&text, "NEXT_ASSEMBLY_USAGE_OCCURRENCE"),
        msb: count(&text, "MANIFOLD_SOLID_BREP"),
        products: 0,
        parts: 0,
        occurrences: 0,
        occ_occurrences: occ_occurrences(path),
        solids: 0,
        refused_solids: 0,
        refused_instances: 0,
        bodies: 0,
        watertight: 0,
        ms: 0,
    };
    let started = popkern_core::time::Stopwatch::start();
    let stem = path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    match Assembly::from_step_bytes(&bytes, &stem) {
        Err(e) => row.ok = Err(e.to_string()),
        Ok((mut assembly, report)) => {
            row.products = report.products;
            row.parts = report.parts;
            row.occurrences = assembly.leaves().len();
            row.solids = report.solids;
            row.refused_solids = report.refused_solids.len();
            row.refused_instances = report.refused_instances.len();
            for line in report
                .refused_solids
                .iter()
                .chain(&report.refused_instances)
                .take(4)
            {
                eprintln!("  {file}: {line}");
            }
            // Per part: a part whose bodies will not validate as parts
            // counts as not watertight (one bad part must not hide the rest).
            for (_, model) in assembly.part_models_mut() {
                match model.to_parts_all() {
                    Ok(parts) => {
                        for part in parts {
                            row.bodies += 1;
                            if let Ok(mesh) = popkern_mesh::tessellate(
                                part.brep(),
                                &popkern_mesh::TessellationOptions::default(),
                            ) && popkern_mesh::welded_open_edges(&mesh) == 0
                            {
                                row.watertight += 1;
                            }
                        }
                    }
                    Err(_) => row.bodies += model.alive_bodies().len().max(1),
                }
            }
        }
    }
    row.ms = started.elapsed().map_or(0, |d| d.as_millis());
    row
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let Some(dir) = args.get(1) else {
        eprintln!("usage: step_corpus <dir> [report.md]");
        std::process::exit(2);
    };
    let mut files: Vec<PathBuf> = std::fs::read_dir(dir)
        .map(|d| d.filter_map(|e| e.ok().map(|e| e.path())).collect())
        .unwrap_or_default();
    files.retain(|p| {
        p.extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("stp") || e.eq_ignore_ascii_case("step"))
    });
    files.sort();
    let mut table = String::from(
        "| file | loaded | products | NAUO in file | part occurrences (ours / OCC) | solids recorded / MSB in file | refused solids | refused occurrences | watertight bodies | ms |\n|---|---|---|---|---|---|---|---|---|---|\n",
    );
    let (mut loaded, mut complete) = (0, 0);
    for path in &files {
        let r = measure(path);
        let occ = r
            .occ_occurrences
            .map(|n| n.to_string())
            .unwrap_or_else(|| "–".into());
        let ok = match &r.ok {
            Ok(()) => {
                loaded += 1;
                if r.refused_solids == 0
                    && r.refused_instances == 0
                    && r.occ_occurrences.is_none_or(|n| n == r.occurrences)
                {
                    complete += 1;
                }
                "yes".to_string()
            }
            Err(e) => format!("no: {}", e.chars().take(80).collect::<String>()),
        };
        table.push_str(&format!(
            "| {} | {} | {} | {} | {} / {} | {} / {} | {} | {} | {}/{} | {} |\n",
            r.file,
            ok,
            r.products,
            r.nauo,
            r.occurrences,
            occ,
            r.solids,
            r.msb,
            r.refused_solids,
            r.refused_instances,
            r.watertight,
            r.bodies,
            r.ms
        ));
    }
    let summary = format!(
        "\n{loaded}/{} files loaded; {complete}/{} complete (every solid recorded, every occurrence placed, occurrence count equal to OpenCASCADE's).\n",
        files.len(),
        files.len()
    );
    table.push_str(&summary);
    print!("{table}");
    if let Some(out) = args.get(2) {
        let _ = std::fs::write(out, &table);
    }
}
