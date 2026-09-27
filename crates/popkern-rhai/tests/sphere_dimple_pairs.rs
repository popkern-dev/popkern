//! PK-455: paired sphere-dimple subtracts on one face, per revolve axis.
//!
//! Subtracting TWO disjoint spheres revolved about the Y (or X) axis
//! into the same block face used to serve validator-clean with a torn
//! mesh (109 open edges): the revolved NURBS sphere stores its angle
//! in v, the chart walk had no v-closure machinery, and the SECOND
//! boolean's re-walk of the surviving cap face closed the rim ring
//! onto itself, stranded the seam pair as a zero-area slit loop, and
//! the roles coin-flipped on zero shoelaces — the trim clip then
//! shredded the cap's grid mesh and the boundary sew snapped 2168 of
//! its vertices onto the rim polyline (the whole dimple sits inside
//! the sew band). Z-revolved spheres always worked: they promote to
//! the exact analytic `Sphere`.
//!
//! Root: geometric v-closure now declines the chart walk to the
//! recorded 3D-Newell fallback exactly like an analytic v-period
//! (`SubdivisionIssue::VPeriodicChart`) — the walk's own doc already
//! held that geometric closure counts like an analytic period.
//!
//! The gate is the mesh/volume oracle (validator-clean is NOT
//! correct — the PK-454 lesson): watertight, and volume within 1% of
//! block minus two spherical caps.

use popkern_rhai::ScriptHost;

const BLOCK: f64 = 0.55;
const RP: f64 = 0.072; // pip sphere radius
const SINK: f64 = 0.052; // sphere centre beyond the face
const G: f64 = 0.15; // dimple offset on the face

/// Two dimples in the +axis face of a block, spheres revolved about
/// `axis`. `fillet` reproduces the 28-open-edge variant of the class.
fn dimple_script(axis: &str, fillet: bool) -> String {
    let h = BLOCK / 2.0;
    let c = h + SINK;
    // Semicircle profile containing the revolve axis, per axis.
    let (profile, spin, center) = match axis {
        "x" => (
            "t.profile(xy()).path(-0.072, 0.0).arc_to([0.0, 0.072], [0.072, 0.0]).line_to(-0.072, 0.0)",
            "[1.0, 0.0, 0.0]",
            format!("{c}, off[0], off[1]"),
        ),
        "y" => (
            "t.profile(xy()).path(0.0, -0.072).arc_to([0.072, 0.0], [0.0, 0.072]).line_to(0.0, -0.072)",
            "[0.0, 1.0, 0.0]",
            format!("off[0], {c}, off[1]"),
        ),
        _ => (
            "t.profile(xz()).path(0.0, -0.072).arc_to([0.072, 0.0], [0.0, 0.072]).line_to(0.0, -0.072)",
            "[0.0, 0.0, 1.0]",
            format!("off[0], off[1], {c}"),
        ),
    };
    let fillet_line = if fillet {
        "b = b.fillet(edges().straight(), 0.05, #{ label: \"rounds\" });"
    } else {
        ""
    };
    format!(
        r#"
let m = model();
let g = {G};
let b = m.profile(xy_at(-{h})).rect({BLOCK}, {BLOCK}).sketch("block").extrude(blind({BLOCK}));
{fillet_line}
let k = 0;
for off in [[g, -1.0 * g], [-1.0 * g, g]] {{
    let ball = m.tool(|t| {{
        let s = {profile}.close().sketch("ball_" + k)
            .revolve(axis_at([0.0, 0.0, 0.0], {spin}), full_circle());
        s.translate({center})
    }});
    b = b.subtract(ball, #{{ label: "cut_" + k }});
    k += 1;
}}
m
"#
    )
}

/// Watertight census + signed volume of the served body's mesh; open
/// segments come back as positions so a failure names WHERE it tore.
fn mesh_volume_and_open_edges(script: &str) -> (f64, Vec<(f64, f64, f64)>) {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(script, "dimples.rhai")
        .unwrap_or_else(|e| panic!("{e}"));
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

/// Analytic: block minus two spherical caps of depth rp - sink.
fn expected_volume() -> f64 {
    let h = RP - SINK;
    let cap = std::f64::consts::PI * h * h * (3.0 * RP - h) / 3.0;
    BLOCK * BLOCK * BLOCK - 2.0 * cap
}

fn assert_dimples_serve(axis: &str, fillet: bool) {
    let (volume, open) = mesh_volume_and_open_edges(&dimple_script(axis, fillet));
    assert!(
        open.is_empty(),
        "[{axis} fillet={fillet}] torn mesh: {} open edges at {:?}",
        open.len(),
        &open[..open.len().min(6)]
    );
    if !fillet {
        let expected = expected_volume();
        let err = (volume - expected).abs() / expected;
        assert!(
            err < 0.01,
            "[{axis}] volume {volume:.6} vs analytic {expected:.6} ({:.2}%)",
            err * 100.0
        );
    }
}

#[test]
fn paired_dimples_from_x_revolved_spheres_serve_watertight() {
    assert_dimples_serve("x", false);
}

#[test]
fn paired_dimples_from_y_revolved_spheres_serve_watertight() {
    assert_dimples_serve("y", false);
}

#[test]
fn paired_dimples_from_y_revolved_spheres_survive_a_prior_fillet() {
    assert_dimples_serve("y", true);
}

#[test]
fn paired_dimples_from_z_revolved_spheres_serve_watertight() {
    assert_dimples_serve("z", false);
}
