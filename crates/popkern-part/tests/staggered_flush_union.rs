//! PK-379 — unions of identical boxes staggered by a fraction of their
//! width: the operands' coplanar faces overlap partially on several
//! parallel plane pairs at once (the "brick" class), so the
//! coincidence engine sees crossing overlaps instead of the contained
//! ones PK-184 closed.
//!
//! Before this task both flush cases failed at the split stage ("trim
//! curve end point does not match any vertex on the face boundary"):
//! a trim lying ALONG a fragment's boundary was routed to the other
//! fragment of the same original face. The x-only stagger is a
//! 1.5 × 1 × 1 box (6 faces); the (0.5, 0.3, 0) stagger is the prism
//! of the two squares' union — an 8-gon — so 8 walls + top + bottom;
//! the generic (0.5, 0.3, 0.2) overlap is the control.
use popkern_part::prelude::*;

fn cube() -> Part {
    Profile::on(Plane::xy_at(0.0))
        .rect(1.0, 1.0)
        .extrude(1.0)
        .expect("unit cube")
}

fn signed_mesh_volume(part: &Part) -> f64 {
    let mesh = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    v6 / 6.0
}

fn assert_union(offset: (f64, f64, f64), faces: usize, volume: f64) -> Part {
    let (dx, dy, dz) = offset;
    let what = format!("union offset ({dx},{dy},{dz})");
    let out = cube()
        .union(&cube().translate(dx, dy, dz).expect("translate"))
        .unwrap_or_else(|e| panic!("{what}: {e}"));
    let issues = out.brep().validate_all();
    assert!(issues.is_empty(), "{what}: validate_all: {issues:?}");
    assert_eq!(
        out.open_edge_count().expect("mesh"),
        0,
        "{what}: watertight"
    );
    let b = out.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "{what}: genus 0");
    let got = signed_mesh_volume(&out);
    assert!(
        (got - volume).abs() < 1e-9,
        "{what}: volume {got:.12}, expected {volume:.12}"
    );
    assert_eq!(out.face_count(), faces, "{what}: coplanar skins merged");
    out
}

#[test]
fn half_width_stagger_flush_on_four_planes_is_one_box() {
    assert_union((0.5, 0.0, 0.0), 6, 1.5);
}

#[test]
fn two_axis_stagger_flush_on_two_planes_is_the_union_prism() {
    // 2-D union of two unit squares offset (0.5, 0.3): area 1 + 1 − 0.35.
    assert_union((0.5, 0.3, 0.0), 10, 1.65);
}

#[test]
fn generic_overlap_control() {
    // Overlap box 0.5 × 0.7 × 0.8.
    assert_union((0.5, 0.3, 0.2), 12, 2.0 - 0.28);
}

/// OCC oracle leg (skips without the OCP venv).
#[test]
fn staggered_unions_pass_the_occ_oracle() {
    let oracle_python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let oracle_script = std::path::Path::new("../../tools/occ-oracle/occ_check.py");
    if !oracle_python.exists() || !oracle_script.exists() {
        eprintln!("[pk379-occ] OCP venv absent — leg skipped");
        return;
    }
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("staggered_flush");
    std::fs::create_dir_all(&dir).expect("dir");
    for (name, offset, faces, volume) in [
        ("stagger_x", (0.5, 0.0, 0.0), 6, 1.5),
        ("stagger_xy", (0.5, 0.3, 0.0), 10, 1.65),
    ] {
        let out = assert_union(offset, faces, volume);
        let path = dir.join(format!("{name}.step"));
        out.export_step(&path).expect("export");
        let back = out.assert_step_roundtrip(&path).expect("round-trip");
        assert_eq!(back.face_count(), faces);
        let output = std::process::Command::new(&oracle_python)
            .arg(oracle_script)
            .arg(&path)
            .output()
            .expect("oracle runs");
        let report = String::from_utf8_lossy(&output.stdout);
        assert!(
            report.contains("\"passed\": true"),
            "{name}: OCC oracle: {}",
            report.chars().take(400).collect::<String>()
        );
    }
}
