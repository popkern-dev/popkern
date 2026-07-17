//! PK-135: chamfer lead-ins around POLYGON end loops — a closed
//! chain of straight-edge chamfers whose strips meet corner-to-corner
//! at every polygon vertex (convex 60° star points and reflex 240°
//! inner corners on the hexagram).
//!
//! Acceptance is measured, not visual: watertight mesh, genus-0
//! Euler characteristic, STEP round-trip, and an exact prismatoid
//! volume oracle for the removed chamfer wedge ring.

use popkern_part::{Dir, EdgeQuery, Part, Plane, Profile, SurfaceKind};

fn hexagram(r_outer: f64) -> Vec<(f64, f64)> {
    (0..12)
        .map(|k| {
            let a = std::f64::consts::TAU * k as f64 / 12.0;
            let r = if k % 2 == 0 {
                r_outer
            } else {
                r_outer / 3.0_f64.sqrt()
            };
            (r * a.cos(), r * a.sin())
        })
        .collect()
}

/// Shoelace area of a polygon.
fn area(poly: &[(f64, f64)]) -> f64 {
    let n = poly.len();
    0.5 * (0..n)
        .map(|i| {
            let (x0, y0) = poly[i];
            let (x1, y1) = poly[(i + 1) % n];
            x0 * y1 - x1 * y0
        })
        .sum::<f64>()
        .abs()
}

/// Inset each edge line by `d` and re-intersect (same construction
/// as the laptop-stand hexagram socket clearance).
fn inset(poly: &[(f64, f64)], d: f64) -> Vec<(f64, f64)> {
    let n = poly.len();
    let lines: Vec<((f64, f64), (f64, f64))> = (0..n)
        .map(|i| {
            let (x0, y0) = poly[i];
            let (x1, y1) = poly[(i + 1) % n];
            let (dx, dy) = (x1 - x0, y1 - y0);
            let len = (dx * dx + dy * dy).sqrt();
            let (nx, ny) = (dy / len, -dx / len);
            ((x0 - nx * d, y0 - ny * d), (dx, dy))
        })
        .collect();
    (0..n)
        .map(|i| {
            let ((px1, py1), (dx1, dy1)) = lines[(i + n - 1) % n];
            let ((px2, py2), (dx2, dy2)) = lines[i];
            let det = dx1 * dy2 - dy1 * dx2;
            let t = ((px2 - px1) * dy2 - (py2 - py1) * dx2) / det;
            (px1 + t * dx1, py1 + t * dy1)
        })
        .collect()
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

fn assert_chamfered_prism(part: &Part, poly: &[(f64, f64)], h: f64, d: f64, ends: usize) {
    let opens = part.open_edge_count().expect("mesh");
    assert_eq!(opens, 0, "chamfered prism must be watertight");

    let b = part.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "chamfered prism must stay genus-0");

    // Exact volume: full prism minus, per chamfered end, the ring
    // between the outline and its inset — a prismatoid of height d
    // (exact for planar-faced solids): V = d/6 · (A_bot + 4·A_mid + A_top).
    let a0 = area(poly);
    let a1 = area(&inset(poly, d));
    let am = area(&inset(poly, d / 2.0));
    let frustum = d / 6.0 * (a0 + 4.0 * am + a1);
    let expected = a0 * h - (a0 * d - frustum) * ends as f64;
    let got = signed_mesh_volume(part);
    assert!(
        (got - expected).abs() < 1.0e-6 * expected.abs().max(1.0),
        "volume oracle: got {got:.6}, expected {expected:.6} (delta {:.2e})",
        got - expected
    );
}

/// One end loop chamfered: the 12 bevel strips must close into a
/// ring with mitre edges at all 12 corners (6 convex star points at
/// 60°, 6 reflex inner corners at 240°).
#[test]
fn hexagram_rim_chamfer_chain_top_is_watertight_and_exact() {
    let poly = hexagram(12.0);
    let h = 30.0;
    let d = 1.5;
    let prism = Profile::on(Plane::XY)
        .polygon(&poly)
        .extrude(h)
        .expect("prism");
    let out = prism
        .chamfer(d, |e: EdgeQuery| {
            e.straight()
                .at_height(h)
                .between(SurfaceKind::Plane, SurfaceKind::Plane)
        })
        .expect("rim chamfer chain");
    assert_chamfered_prism(&out, &poly, h, d, 1);
    let dir = std::env::temp_dir().join("pk135_rim_chamfer");
    std::fs::create_dir_all(&dir).expect("tmp dir");
    let path = dir.join("hexagram_rim.step");
    out.export_step(&path).expect("step export");
    let re = out.assert_step_roundtrip(&path).expect("roundtrip gate");
    assert_eq!(re.open_edge_count().expect("mesh"), 0, "roundtrip mesh");
}

/// Both end loops chamfered — the peg's actual feature chain.
#[test]
fn hexagram_rim_chamfer_chain_both_ends_is_watertight_and_exact() {
    let poly = hexagram(12.0);
    let h = 30.0;
    let d = 1.5;
    let prism = Profile::on(Plane::XY)
        .polygon(&poly)
        .extrude(h)
        .expect("prism");
    let out = prism
        .chamfer(d, |e: EdgeQuery| {
            e.straight()
                .at_height(h)
                .between(SurfaceKind::Plane, SurfaceKind::Plane)
        })
        .expect("top rim")
        .chamfer(d, |e: EdgeQuery| {
            e.straight()
                .at_height(0.0)
                .between(SurfaceKind::Plane, SurfaceKind::Plane)
        })
        .expect("bottom rim");
    assert_chamfered_prism(&out, &poly, h, d, 2);
}

/// Convex-only control: a square rim chain (four 90° corners) — the
/// simplest closed bevel-bevel chain. If this fails the defect is in
/// chain corners generally, not the star's 60°/240° extremes.
#[test]
fn square_rim_chamfer_chain_is_watertight_and_exact() {
    let poly: Vec<(f64, f64)> = vec![(-10.0, -10.0), (10.0, -10.0), (10.0, 10.0), (-10.0, 10.0)];
    let h = 20.0;
    let d = 2.0;
    let prism = Profile::on(Plane::XY)
        .polygon(&poly)
        .extrude(h)
        .expect("prism");
    let out = prism
        .chamfer(d, |e: EdgeQuery| {
            e.straight()
                .at_height(h)
                .between(SurfaceKind::Plane, SurfaceKind::Plane)
        })
        .expect("rim chamfer chain");
    assert_chamfered_prism(&out, &poly, h, d, 1);
    let _ = Dir::Z;
}

/// Diagnostic probe (run with --ignored --nocapture): dump the
/// square-chain outcome — open edges, faces, per-face surfaces.
#[test]
#[ignore = "diagnostic probe for PK-135 corner failure analysis"]
fn square_rim_chamfer_chain_probe() {
    let poly: Vec<(f64, f64)> = vec![(-10.0, -10.0), (10.0, -10.0), (10.0, 10.0), (-10.0, 10.0)];
    let prism = Profile::on(Plane::XY)
        .polygon(&poly)
        .extrude(20.0)
        .expect("prism");
    let out = prism
        .chamfer(2.0, |e: EdgeQuery| {
            e.straight()
                .at_height(20.0)
                .between(SurfaceKind::Plane, SurfaceKind::Plane)
        })
        .expect("rim chamfer chain");
    eprintln!("open mesh edges: {:?}", out.open_edge_count());
    let b = out.brep();
    eprintln!(
        "V={} E={} F={} L={}",
        b.vertices.len(),
        b.edges.len(),
        b.faces.len(),
        b.loops.len()
    );
    for (fh, fd) in &b.faces {
        let kind = match fd.surface.as_ref().map(|s| &s.geometry) {
            Some(popkern_core::geometry::SurfaceGeometry::Plane { origin, normal, .. }) => format!(
                "Plane o=({:.1},{:.1},{:.1}) n=({:.2},{:.2},{:.2})",
                origin.x, origin.y, origin.z, normal.x, normal.y, normal.z
            ),
            Some(_) => "other".into(),
            None => "none".into(),
        };
        let pts = b.face_outer_loop_positions(fh);
        eprintln!("  {fh:?} verts={} {kind}", pts.len());
    }
    for err in b.validate_all() {
        eprintln!("  validate: {err}");
    }
}
