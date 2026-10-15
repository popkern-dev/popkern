//! PK-378 — planar chamfers meeting three at a corner build the exact
//! mitre apex.
//!
//! The equal-distance chamfer of every edge of a cube is the
//! half-space intersection: each bevel is a hexagon ending in a point
//! where the three bevel planes concur (`(d/2, d/2, d/2)` from the
//! corner), so the body has 6 + 12 = 18 faces and no corner triangle.
//! Volume by inclusion–exclusion over the twelve edge prisms
//! (`d²/2` each), the 24 adjacent-pair overlaps (`d³/3`) and the eight
//! triple overlaps (`d³/4`): `V = 1 − 6d² + 6d³` for `d < 1/2`. One
//! corner alone: `V = 1 − 3d²/2 + 3d³/4`, 9 faces.
//!
//! Before this task the 12-edge chamfer produced an invalid body (an
//! edge `d/√2` off its bevel plane) for every `d`, and a single corner
//! served 12 faces with the corner sub-cube left as material.
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

fn assert_chamfered_body(part: &Part, faces: usize, expected_volume: f64, what: &str) {
    assert_eq!(part.face_count(), faces, "{what}: face count");
    let issues = part.brep().validate_all();
    assert!(issues.is_empty(), "{what}: validate_all: {issues:?}");
    let opens = part.open_edge_count().expect("mesh");
    assert_eq!(opens, 0, "{what}: must be watertight");
    let b = part.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "{what}: genus 0");
    let got = signed_mesh_volume(part);
    assert!(
        (got - expected_volume).abs() < 1.0e-9,
        "{what}: volume {got:.12}, expected {expected_volume:.12} (delta {:.2e})",
        got - expected_volume
    );
}

#[test]
fn every_cube_edge_chamfered_is_the_exact_half_space_intersection() {
    for d in [0.05, 0.1, 0.2, 0.3] {
        let out = cube()
            .chamfer(d, |e: EdgeQuery| e)
            .unwrap_or_else(|e| panic!("12-edge chamfer d={d}: {e}"));
        let expected = 1.0 - 6.0 * d * d + 6.0 * d * d * d;
        assert_chamfered_body(&out, 18, expected, &format!("12-edge chamfer d={d}"));
        // Every bevel is a hexagon: two rails, two mitre edges per end.
        let b = out.brep();
        let mut hexagons = 0;
        for (fh, _) in b.faces.iter() {
            let n = b.face_outer_loop_positions(fh).len();
            match n {
                6 => hexagons += 1,
                4 => {}
                other => panic!("d={d}: face {fh:?} has {other} boundary vertices"),
            }
        }
        assert_eq!(hexagons, 12, "d={d}: twelve hexagonal bevels");
        // The apex sits at (d/2, d/2, d/2) from each corner.
        let apexes = b
            .vertices
            .iter()
            .filter(|(_, v)| {
                // Profile rect is centred in x/y; the extrude spans z ∈ [0, 1].
                let p = v.position;
                let at = |c: f64, lo: f64, hi: f64| {
                    (c - (lo + d / 2.0)).abs() < 1e-12 || (c - (hi - d / 2.0)).abs() < 1e-12
                };
                at(p.x, -0.5, 0.5) && at(p.y, -0.5, 0.5) && at(p.z, 0.0, 1.0)
            })
            .count();
        assert_eq!(apexes, 8, "d={d}: one mitre apex per corner");
    }
}

#[test]
fn single_corner_three_chamfers_meet_at_the_mitre_apex() {
    let d = 0.2;
    // The three edges at the (-0.5,-0.5,0) corner.
    let out = cube()
        .chamfer(d, |e: EdgeQuery| {
            e.filter(|info| {
                // Edge midpoints: vertical (-0.5,-0.5,0.5); bottom
                // (-0.5,0,0) and (0,-0.5,0).
                let m = info.midpoint;
                let near = |a: f64, b: f64| (a - b).abs() < 1e-9;
                (near(m.x, -0.5) && near(m.y, -0.5))
                    || (near(m.x, -0.5) && near(m.z, 0.0))
                    || (near(m.y, -0.5) && near(m.z, 0.0))
            })
        })
        .expect("three chamfers at one corner");
    let expected = 1.0 - 1.5 * d * d + 0.75 * d * d * d;
    assert_chamfered_body(&out, 9, expected, "single trihedral corner");
}

/// OCC oracle leg (skips without the OCP venv): the 12-edge chamfer
/// round-trips through OpenCASCADE's BRepCheck / free-bounds check.
#[test]
fn every_cube_edge_chamfered_passes_the_occ_oracle() {
    let oracle_python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let oracle_script = std::path::Path::new("../../tools/occ-oracle/occ_check.py");
    if !oracle_python.exists() || !oracle_script.exists() {
        eprintln!("[pk378-occ] OCP venv absent — leg skipped");
        return;
    }
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("trihedral_mitre");
    std::fs::create_dir_all(&dir).expect("dir");
    for d in [0.1, 0.3] {
        let out = cube()
            .chamfer(d, |e: EdgeQuery| e)
            .expect("12-edge chamfer");
        let path = dir.join(format!("cube_chamfer_{d}.step"));
        out.export_step(&path).expect("export");
        let back = out.assert_step_roundtrip(&path).expect("round-trip");
        assert_eq!(back.face_count(), 18);
        let output = std::process::Command::new(&oracle_python)
            .arg(oracle_script)
            .arg(&path)
            .output()
            .expect("oracle runs");
        let report = String::from_utf8_lossy(&output.stdout);
        assert!(
            report.contains("\"passed\": true"),
            "d={d}: OCC oracle: {}",
            report.chars().take(400).collect::<String>()
        );
    }
}
