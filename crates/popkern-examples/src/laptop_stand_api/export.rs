//! Piece acceptance gates and output artifacts for the laptop stand:
//! per-piece watertight + STEP round-trip gating with STEP/STL/OBJ
//! export, and the merged assembly-preview OBJ writer.

use popkern_part::prelude::*;

/// Write a merged multi-piece OBJ (meshed copies, each yawed about
/// z then translated — a visual stack preview, not a boolean
/// assembly).
pub(crate) fn write_assembly_preview(
    pieces: &[(&Part, [f64; 3], f64)],
    path: &std::path::Path,
) -> std::io::Result<()> {
    use std::io::Write as _;
    let mut out = String::new();
    let mut base_index = 1usize;
    let options = popkern_mesh::TessellationOptions::default();
    for (part, offset, yaw) in pieces {
        let mesh = part.mesh(&options).expect("preview mesh");
        let (c, s) = (yaw.cos(), yaw.sin());
        for v in &mesh.vertices {
            out.push_str(&format!(
                "v {} {} {}\n",
                v.position.x * c - v.position.y * s + offset[0],
                v.position.x * s + v.position.y * c + offset[1],
                v.position.z + offset[2]
            ));
        }
        for t in &mesh.triangles {
            out.push_str(&format!(
                "f {} {} {}\n",
                t[0] + base_index,
                t[1] + base_index,
                t[2] + base_index
            ));
        }
        base_index += mesh.vertices.len();
    }
    let mut f = std::fs::File::create(path)?;
    f.write_all(out.as_bytes())
}

/// Gate a finished piece (non-empty, watertight, STEP round-trip)
/// and export STEP + STL + OBJ. `PK_OPEN_DEBUG=1` dumps the first
/// open mesh edges when the watertight gate is about to fail.
pub(crate) fn gate_and_export(part: Part, name: &str, out_dir: &std::path::Path) -> Part {
    if std::env::var_os("PK_OPEN_DEBUG").is_some()
        && let Ok(n) = part.open_edge_count()
        && n > 0
    {
        let mesh = part
            .mesh(&popkern_mesh::TessellationOptions::default())
            .expect("mesh");
        use std::collections::HashMap;
        let mut count: HashMap<(usize, usize), i32> = HashMap::new();
        for t in &mesh.triangles {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                *count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
        let mut shown = 0;
        for ((a, b), c) in &count {
            if *c != 2 && shown < 10 {
                let pa = &mesh.vertices[*a].position;
                let pb = &mesh.vertices[*b].position;
                eprintln!(
                    "  [{name}-open] uses={c} ({:.2},{:.2},{:.2})->({:.2},{:.2},{:.2})",
                    pa.x, pa.y, pa.z, pb.x, pb.y, pb.z
                );
                shown += 1;
            }
        }
    }
    // An EMPTY body meshes trivially watertight — a silently swallowed
    // boolean must not pass the gate (measured: an arm tip whose run
    // had wedged out to nothing exported as a valid 0-face part).
    assert!(part.face_count() > 0, "{name}: empty body");
    let part = part.assert_watertight().expect("watertight gate");
    let step = out_dir.join(format!("{name}.step"));
    part.export_step(&step).expect("STEP export");
    part.export_stl(out_dir.join(format!("{name}.stl")))
        .expect("STL export");
    part.export_obj(out_dir.join(format!("{name}.obj")))
        .expect("OBJ export");
    let part = part.assert_step_roundtrip(&step).expect("round-trip gate");
    println!(
        "  {name}: watertight ✓ round-trip ✓ ({} faces) → {}",
        part.face_count(),
        step.display()
    );
    part
}
