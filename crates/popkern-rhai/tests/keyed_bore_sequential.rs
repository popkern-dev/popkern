//! PK-428: a keyway slot subtracted SEQUENTIALLY after its bore (not
//! pre-unioned into one tool) serves watertight.
//!
//! The slot's four wall records chain into a closed rectangle that
//! STRADDLES the bore hole on each cap face — one side wholly inside
//! the bore, two sides crossing the rim mid-record. The fragment walk
//! owns the fix: record pieces lying inside a face's own holes are
//! phantoms for that face (they bound the other operand's material)
//! and are dropped from its walk list, and a ring partially twinned
//! with an already-walked hole ring is the region piece the hole
//! bounds — the winding rule elects it a separate fragment. Before:
//! the boolean served a corrupt body — 74 open mesh edges, cellular
//! declined on a non-manifold edge, six degenerate junction stubs.
//!
//! Gate: mesh/volume oracle (validator-clean is NOT correct — the
//! PK-454 lesson): watertight and volume within 0.5% of analytic.

use popkern_rhai::ScriptHost;

const BLOCK_H: f64 = 0.3;
const BORE_R: f64 = 0.15;

/// Block, bore subtract, then a keyway slot subtract whose slot opens
/// into the bore. `rect` = (width, height), `at` = slot centre.
fn keyed_bore_script(rect: (f64, f64), at: (f64, f64)) -> String {
    let (w, h) = rect;
    let (cx, cy) = at;
    format!(
        r#"
let m = model();
let b = m.profile(xy()).rect(1.0, 1.0).sketch("block").extrude(blind({BLOCK_H}));
b = b.subtract(m.tool(|t| t.profile(xy_at(-0.1)).circle({BORE_R}).sketch("bore").extrude(blind(0.5))), #{{ label: "bore_cut" }});
b = b.subtract(m.tool(|t| t.profile(xy_at(-0.14)).rect({w}, {h}).at({cx}, {cy}).sketch("key").extrude(blind(0.58))), #{{ label: "key_cut" }});
m.assert_watertight();
m
"#
    )
}

/// Watertight census + signed volume of the served body's mesh; open
/// segments come back as positions so a failure names WHERE it tore.
fn mesh_volume_and_open_edges(script: &str) -> (f64, Vec<(f64, f64, f64)>) {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(script, "keyed_bore.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
    let mut model = outcome.model.expect("model");
    let part = model.to_part().expect("part");
    let brep = part.brep();
    brep.validate().expect("validation clean");
    let options = popkern_mesh::TessellationOptions {
        chord_height: 5.0e-4,
        min_subdivisions: 4,
        max_subdivisions: 512,
        ..Default::default()
    };
    let mesh = popkern_mesh::tessellate(brep, &options).expect("mesh");
    let mut vol6 = 0.0;
    for t in &mesh.triangles {
        let a = mesh.vertices[t[0]].position;
        let b = mesh.vertices[t[1]].position;
        let c = mesh.vertices[t[2]].position;
        vol6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    type QuantizedPoint = (i64, i64, i64);
    let mut edge_uses: std::collections::HashMap<(QuantizedPoint, QuantizedPoint), i64> =
        std::collections::HashMap::new();
    let key = |p: popkern_core::geometry::Point3| -> (i64, i64, i64) {
        let q = |x: f64| (x * 1.0e7).round() as i64;
        (q(p.x), q(p.y), q(p.z))
    };
    for t in &mesh.triangles {
        for (i, j) in [(0usize, 1usize), (1, 2), (2, 0)] {
            let a = key(mesh.vertices[t[i]].position);
            let b = key(mesh.vertices[t[j]].position);
            *edge_uses.entry((a.min(b), a.max(b))).or_insert(0) += if a < b { 1 } else { -1 };
        }
    }
    let open: Vec<(f64, f64, f64)> = edge_uses
        .iter()
        .filter(|&(_, &c)| c != 0)
        .map(|((a, _), _)| (a.0 as f64 / 1e7, a.1 as f64 / 1e7, a.2 as f64 / 1e7))
        .collect();
    (vol6.abs() / 6.0, open)
}

/// Analytic: block minus (bore ∪ slot) footprint times height. The
/// circle/rect overlap integrates numerically at 2e5 stations —
/// deterministic and converged far below the mesh-volume gate.
fn expected_volume(rect: (f64, f64), at: (f64, f64)) -> f64 {
    let (w, h) = rect;
    let (cx, cy) = at;
    // Clamp the slot footprint to the block (an edge keyway's tool
    // overshoots the stock).
    let (x0, x1) = ((cx - w * 0.5).max(-0.5), (cx + w * 0.5).min(0.5));
    let (y0, y1) = ((cy - h * 0.5).max(-0.5), (cy + h * 0.5).min(0.5));
    let (w, h) = (x1 - x0, y1 - y0);
    let mut overlap = 0.0;
    let n = 200_000;
    for i in 0..n {
        let y = y0 + (y1 - y0) * (i as f64 + 0.5) / n as f64;
        if y.abs() >= BORE_R {
            continue;
        }
        let half = (BORE_R * BORE_R - y * y).sqrt();
        let lo = x0.max(-half);
        let hi = x1.min(half);
        if hi > lo {
            overlap += hi - lo;
        }
    }
    overlap *= (y1 - y0) / n as f64;
    let union = std::f64::consts::PI * BORE_R * BORE_R + w * h - overlap;
    (1.0 - union) * BLOCK_H
}

fn assert_keyed_bore_serves(rect: (f64, f64), at: (f64, f64), label: &str) {
    let (volume, open) = mesh_volume_and_open_edges(&keyed_bore_script(rect, at));
    assert!(
        open.is_empty(),
        "[{label}] torn mesh: {} open edges at {:?}",
        open.len(),
        &open[..open.len().min(6)]
    );
    let expected = expected_volume(rect, at);
    let err = (volume - expected).abs() / expected;
    assert!(
        err < 0.005,
        "[{label}] volume {volume:.6} vs analytic {expected:.6} ({:.2}%)",
        err * 100.0
    );
}

#[test]
fn keyway_slot_into_bore_serves_watertight() {
    // The PK-428 minimal repro: slot opens into the bore from +y.
    assert_keyed_bore_serves((0.08, 0.06), (0.0, 0.17), "slot +y");
}

#[test]
fn keyway_slot_into_bore_rotated_serves_watertight() {
    // Same class rotated a quarter turn — the fix must not be axis-lucky.
    assert_keyed_bore_serves((0.06, 0.08), (0.17, 0.0), "slot +x");
}

#[test]
fn keyway_slot_from_block_edge_into_bore_serves_watertight() {
    // The mixed-boundary sibling: the slot runs from the block's
    // OUTER edge into the bore, so its surviving records arrive as
    // open chains (the outer clip truncates them at the block edge)
    // and the open-chain machinery owns the cut. Pinned alongside the
    // closed-chain class so both keyway shapes stay covered.
    assert_keyed_bore_serves((0.08, 0.44), (0.0, 0.35), "edge keyway");
}

#[test]
fn wide_slot_crossing_the_bore_serves_watertight() {
    // The slot spans the whole bore: the chain crosses the rim FOUR
    // times and the region clip must emit TWO rim-to-rim runs.
    assert_keyed_bore_serves((0.44, 0.06), (0.0, 0.17), "slot across");
}

/// The CURVED-wall sibling: a radial bore through a tube wall, then a
/// radial keyway slot overlapping that bore. The slot records straddle
/// the hole on the periodic wall chart — the walk's hole-interior
/// record filter must branch-normalize its containment there (before:
/// validator-clean serve with 89 open mesh edges; the volume was even
/// right, the wall mesh just tore around the phantom rings).
#[test]
fn radial_keyway_into_tube_bore_serves_watertight() {
    const SCRIPT: &str = r#"
let m = model();
let b = m.profile(xy()).circle(0.3).sketch("tube").extrude(blind(0.6));
b = b.subtract(m.tool(|t| t.profile(plane_at([0.4, 0.0, 0.3], [-1.0, 0.0, 0.0])).circle(0.05).sketch("bore").extrude(blind(0.8))), #{ label: "bore_cut" });
b = b.subtract(m.tool(|t| t.profile(plane_at([0.4, 0.0, 0.36], [-1.0, 0.0, 0.0])).rect(0.05, 0.04).sketch("key").extrude(blind(0.8))), #{ label: "key_cut" });
m.assert_watertight();
m
"#;
    let (volume, open) = mesh_volume_and_open_edges(SCRIPT);
    assert!(
        open.is_empty(),
        "[radial keyway] torn mesh: {} open edges at {:?}",
        open.len(),
        &open[..open.len().min(6)]
    );
    // Grid CSG volume: tube minus the x-invariant bore and slot
    // profiles (both cutters run parallel to x, so each removed (y,z)
    // cell removes the tube's whole x-chord).
    let n = 400usize;
    let step = 0.62 / n as f64;
    let mut expected = 0.0;
    for iy in 0..n {
        let y = -0.31 + (iy as f64 + 0.5) * step;
        for iz in 0..n {
            let z = -0.01 + (iz as f64 + 0.5) * step;
            if !(0.0..=0.6).contains(&z) {
                continue;
            }
            let r2 = 0.09 - y * y;
            if r2 <= 0.0 {
                continue;
            }
            let in_bore = y * y + (z - 0.3) * (z - 0.3) < 0.0025;
            let in_slot = y.abs() < 0.025 && (0.34..0.38).contains(&z);
            if !(in_bore || in_slot) {
                expected += 2.0 * r2.sqrt() * step * step;
            }
        }
    }
    let err = (volume - expected).abs() / expected;
    assert!(
        err < 0.01,
        "[radial keyway] volume {volume:.6} vs grid CSG {expected:.6} ({:.2}%)",
        err * 100.0
    );
}
