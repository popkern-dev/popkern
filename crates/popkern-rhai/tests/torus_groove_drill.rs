//! PK-459 (partial): drills through a v-partial torus GROOVE face
//! serve watertight.
//!
//! A torus-ring tool subtracted from a cylinder leaves a groove whose
//! wall is a v-bounded band of a doubly-periodic surface (the bearing
//! race construction). Drilling through that band used to DECLINE
//! ("assembly stranded half-claimed record sides"): the chart walk
//! refused every analytic-torus face whose chart cell could not
//! resolve, and the 3D-Newell fallback elected a drill rim as the
//! band's outer. The fix routes cell-less but v-BOUNDED faces to the
//! tuned u-periodic chart machinery — for the walk they are exactly
//! u-closed walls.
//!
//! The groove-only volume is pinned too: authoring the torus tool as
//! split half-tube spans was measured to glue a tool protrusion into
//! the result (volume 0.185 vs 0.1658) — the tool-side machinery
//! relies on the doubly-wrapped single-face form (see the PK-459
//! design-fork ledger record).
//!
//! Census quantum is 1e-9: the sew leaves legitimate distinct
//! boundary samples ~6e-8 apart, and a 1e-7 census collapses them
//! into false self-edge "cracks" (the PK-325 lesson).

use popkern_rhai::ScriptHost;

fn groove_script(drill_z: Option<f64>) -> String {
    let drill = match drill_z {
        Some(z) => format!(
            r#"b = b.subtract(m.tool(|t| t.profile(plane_at([0.4, 0.0, {z}], [-1.0, 0.0, 0.0])).circle(0.02).sketch("bore").extrude(blind(0.8))), #{{ label: "drill" }});"#
        ),
        None => String::new(),
    };
    format!(
        r#"
let m = model();
let b = m.profile(xy()).circle(0.3).sketch("cyl").extrude(blind(0.6));
let groove = m.tool(|t| t.profile(xz()).circle(0.05).at(0.32, 0.3).sketch("ring")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle()));
b = b.subtract(groove, #{{ label: "groove_cut" }});
{drill}
m
"#
    )
}

/// Mesh volume + open census at 1e-9 quantum.
fn census(script: &str) -> (f64, usize) {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(script, "groove.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
    let mut model = outcome.model.expect("model");
    let part = model.to_part().expect("part");
    let brep = part.brep();
    let options = popkern_mesh::TessellationOptions {
        chord_height: 5.0e-4,
        min_subdivisions: 4,
        max_subdivisions: 512,
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
    let mut uses: std::collections::HashMap<((i64, i64, i64), (i64, i64, i64)), i64> =
        Default::default();
    let key = |p: popkern_core::geometry::Point3| {
        (
            (p.x * 1e9).round() as i64,
            (p.y * 1e9).round() as i64,
            (p.z * 1e9).round() as i64,
        )
    };
    for t in &mesh.triangles {
        for (i, j) in [(0usize, 1usize), (1, 2), (2, 0)] {
            let a = key(mesh.vertices[t[i]].position);
            let b = key(mesh.vertices[t[j]].position);
            *uses.entry((a.min(b), a.max(b))).or_insert(0) += if a < b { 1 } else { -1 };
        }
    }
    let open = uses.values().filter(|&&c| c != 0).count();
    (vol6.abs() / 6.0, open)
}

/// Grid CSG: cylinder r 0.3 h 0.6 minus torus (R 0.32, r 0.05, plane
/// z = 0.3) minus an optional x-axis drill (r 0.02 at y = 0, z).
fn expected(drill_z: Option<f64>) -> f64 {
    let n = 260usize;
    let step = 0.62 / n as f64;
    let mut vol = 0.0;
    for iy in 0..n {
        let y = -0.31 + (iy as f64 + 0.5) * step;
        for iz in 0..n {
            let z = -0.01 + (iz as f64 + 0.5) * step;
            if !(0.0..=0.6).contains(&z) {
                continue;
            }
            // x-extent inside the cylinder at this (y, z).
            let r2 = 0.09 - y * y;
            if r2 <= 0.0 {
                continue;
            }
            let half = r2.sqrt();
            // Removed x-intervals along this scan line.
            let mut kept = 2.0 * half;
            // Torus: (sqrt(x^2+y^2) - 0.32)^2 + (z-0.3)^2 < 0.05^2 —
            // solve for |x| interval: sqrt(x^2+y^2) in [0.32-w, 0.32+w].
            let w2 = 0.0025 - (z - 0.3) * (z - 0.3);
            if w2 > 0.0 {
                let w = w2.sqrt();
                let (rho0, rho1) = (0.32 - w, 0.32 + w);
                let x_of = |rho: f64| (rho * rho - y * y).max(0.0).sqrt();
                if rho1 * rho1 > y * y {
                    let x1 = x_of(rho1).min(half);
                    let x0 = if rho0 * rho0 > y * y { x_of(rho0) } else { 0.0 };
                    let x0 = x0.min(half);
                    if x1 > x0 {
                        kept -= 2.0 * (x1 - x0);
                    }
                }
            }
            if let Some(dz) = drill_z {
                if y * y + (z - dz) * (z - dz) < 0.0004 {
                    kept = 0.0_f64.max(kept - 2.0 * half.min(half));
                }
            }
            vol += kept.max(0.0) * step * step;
        }
    }
    vol
}

fn assert_serves(drill_z: Option<f64>, label: &str) {
    let (volume, open) = census(&groove_script(drill_z));
    assert_eq!(open, 0, "[{label}] torn mesh: {open} open segments");
    let want = expected(drill_z);
    let err = (volume - want).abs() / want;
    assert!(
        err < 0.01,
        "[{label}] volume {volume:.6} vs grid CSG {want:.6} ({:.2}%)",
        err * 100.0
    );
}

#[test]
fn torus_groove_alone_serves_with_true_volume() {
    assert_serves(None, "groove only");
}

#[test]
fn drill_through_torus_groove_serves_watertight() {
    assert_serves(Some(0.3), "drill mid-groove");
}

#[test]
fn asymmetric_drill_through_torus_groove_serves_watertight() {
    assert_serves(Some(0.31), "drill off-centre");
}
