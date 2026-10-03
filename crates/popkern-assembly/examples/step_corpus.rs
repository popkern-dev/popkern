//! STEP assembly import corpus measurement (PK-631, design §8.3).
//!
//!     cargo run --release -p popkern-assembly --example step_corpus -- <dir> [report.md]
//!
//! For every `.stp`/`.step` file in `<dir>`: import it as an assembly
//! and record what loaded — products, occurrences placed against the
//! file's `NEXT_ASSEMBLY_USAGE_OCCURRENCE` count, solids recorded
//! against the file's `MANIFOLD_SOLID_BREP` count, refused solids and
//! occurrences, part occurrences against OpenCASCADE's (XCAF, when the
//! OCP oracle is installed), how many part bodies tessellate
//! watertight, and how many of those match a solid OpenCASCADE reads
//! from the file (same face count, volume within 1.5 %; PK-649 — a
//! watertight mesh of the wrong region or winding is not correct).
//! Prints a markdown table (and writes it to the report path when
//! given).

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
    /// Watertight bodies matching an OCC solid; `None` without the oracle.
    volume_matched: Option<usize>,
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

/// Every solid OpenCASCADE reads from the file as (faces, volume);
/// `None` when the OCP oracle is not installed.
fn occ_solids(path: &Path) -> Option<Vec<(usize, f64)>> {
    let python = PathBuf::from(std::env::var("HOME").ok()?).join("venvs/ocp/bin/python");
    let script = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle/occ_solids.py");
    if !python.exists() {
        return None;
    }
    let out = Command::new(python).arg(script).arg(path).output().ok()?;
    let report: serde_json::Value = serde_json::from_slice(&out.stdout).ok()?;
    let rows = report["solids"].as_array()?;
    Some(
        rows.iter()
            .filter_map(|r| Some((r["faces"].as_u64()? as usize, r["volume"].as_f64()?)))
            .collect(),
    )
}

/// The mesh's enclosed volume (outward winding positive).
fn mesh_volume(mesh: &popkern_mesh::TriMesh) -> f64 {
    mesh.triangles
        .iter()
        .map(|t| {
            let (a, b, c) = (
                mesh.vertices[t[0]].position,
                mesh.vertices[t[1]].position,
                mesh.vertices[t[2]].position,
            );
            (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0
        })
        .sum()
}

/// Whether a watertight body matches one of OCC's solids: the same
/// face count and a mesh volume within 1.5 %. Small parts whose chord
/// sag alone exceeds that at the default chord height are re-measured
/// at a tenth of it (4pinplug's pins: 8 % low at 0.005, matched at
/// 0.0005).
fn matches_occ(brep: &popkern_core::topology::BRep, volume: f64, occ: &[(usize, f64)]) -> bool {
    let faces = brep.faces.len();
    let near = |v: f64| {
        occ.iter()
            .any(|&(f, o)| f == faces && (v - o).abs() <= 0.015 * o.abs().max(1e-12))
    };
    if near(volume) {
        return true;
    }
    let fine = popkern_mesh::TessellationOptions {
        chord_height: 0.1 * popkern_mesh::TessellationOptions::default().chord_height,
        ..Default::default()
    };
    popkern_mesh::tessellate(brep, &fine).is_ok_and(|m| near(mesh_volume(&m)))
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
        volume_matched: None,
        ms: 0,
    };
    let occ = occ_solids(path);
    let mut matched = 0;
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
                                if let Some(occ) = &occ {
                                    let volume = mesh_volume(&mesh);
                                    if matches_occ(part.brep(), volume, occ) {
                                        matched += 1;
                                    } else {
                                        eprintln!(
                                            "  {file}: watertight body with {} faces, volume {volume:.6}, matches no OpenCASCADE solid",
                                            part.brep().faces.len()
                                        );
                                    }
                                }
                            }
                        }
                    }
                    Err(_) => row.bodies += model.alive_bodies().len().max(1),
                }
            }
        }
    }
    row.volume_matched = occ.map(|_| matched);
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
        "| file | loaded | products | NAUO in file | part occurrences (ours / OCC) | solids recorded / MSB in file | refused solids | refused occurrences | watertight bodies | of which OCC-volume-matched | ms |\n|---|---|---|---|---|---|---|---|---|---|---|\n",
    );
    let (mut loaded, mut complete) = (0, 0);
    let (mut watertight, mut matched, mut bodies) = (0, 0, 0);
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
        watertight += r.watertight;
        bodies += r.bodies;
        matched += r.volume_matched.unwrap_or(0);
        table.push_str(&format!(
            "| {} | {} | {} | {} | {} / {} | {} / {} | {} | {} | {}/{} | {} | {} |\n",
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
            r.volume_matched
                .map(|m| format!("{m}/{}", r.watertight))
                .unwrap_or_else(|| "–".into()),
            r.ms
        ));
    }
    let summary = format!(
        "\n{loaded}/{} files loaded; {complete}/{} complete (every solid recorded, every occurrence placed, occurrence count equal to OpenCASCADE's); {watertight}/{bodies} bodies watertight, {matched} of them matching an OpenCASCADE solid's volume.\n",
        files.len(),
        files.len()
    );
    table.push_str(&summary);
    print!("{table}");
    if let Some(out) = args.get(2) {
        let _ = std::fs::write(out, &table);
    }
}
