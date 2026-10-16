//! A sphere face keeps a valid chart under rotation.
//!
//! A sphere's chart (poles along its axis, u from the first vector of
//! `orthonormal_basis(axis)`) must turn with the material, as cylinder,
//! cone and torus charts do. While `Sphere` stored only centre and
//! radius its chart stayed world-fixed: a dome turned on its side kept
//! pcurves for its upright pose (1.41 off), and once re-derived, its
//! base circle ran through both chart poles, which no continuous chart
//! curve can follow (0.097 off, 8 open edges in the display mesh); a
//! tilted dome put a pole inside its face (8 open edges). The six-axis
//! arm's j5 housing showed it: its world copy meshed through the
//! uniform-grid fallback in 44 s with 512 open edges. The base disc
//! also caught a second defect: `attach_planar_pcurves` took a closed
//! edge's sense from its (tied) endpoints, so the turned disc wound
//! against its orientation flag.

use std::collections::HashMap;

use popkern_rhai::ScriptHost;

fn domed(turn: &str) -> popkern_core::topology::BRep {
    let source = format!(
        r#"let m = model();
let body = m.profile(xz())
    .path(0.0, 0.0).line_to(1.0, 0.0).line_to(1.0, 2.0)
    .arc_to([0.70710678118654752, 2.70710678118654752], [0.0, 3.0])
    .close()
    .sketch("dome")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
{turn}
m
"#
    );
    let mut run = ScriptHost::new()
        .run(&source, "dome.rhai")
        .expect("dome builds");
    let mut model = run.model.take().expect("model");
    let parts = model.to_parts_all().expect("parts");
    parts[0].brep().clone()
}

/// Open and non-manifold edges of the mesh welded at bit-identical positions.
fn welded_open(mesh: &popkern_mesh::TriMesh) -> usize {
    let key = |p: &popkern_core::geometry::Point3| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
    let mut weld: HashMap<[u64; 3], usize> = HashMap::new();
    let ids: Vec<usize> = mesh
        .vertices
        .iter()
        .map(|v| {
            let n = weld.len();
            *weld.entry(key(&v.position)).or_insert(n)
        })
        .collect();
    let mut edges: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (ids[t[k]], ids[t[(k + 1) % 3]]);
            *edges.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    edges.values().filter(|&&n| n != 2).count()
}

#[test]
fn a_dome_turned_on_its_side_validates_and_meshes_watertight() {
    for (label, turn) in [
        ("upright", ""),
        (
            "spun",
            "body = body.rotate(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), 50.0);",
        ),
        (
            "on its side",
            "body = body.rotate(axis_at([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), 90.0);",
        ),
        (
            "tilted",
            "body = body.rotate(axis_at([0.3, -0.2, 0.5], [0.6, 0.0, 0.8]), 37.0);",
        ),
    ] {
        let brep = domed(turn);
        let issues = brep.validate_all();
        assert!(issues.is_empty(), "{label}: validate_all: {issues:?}");
        let watch = std::time::Instant::now();
        let mesh = popkern_mesh::tessellate(&brep, &popkern_mesh::TessellationOptions::display())
            .expect("mesh");
        let elapsed = watch.elapsed();
        assert_eq!(
            welded_open(&mesh),
            0,
            "{label}: watertight in the bit-identical weld"
        );
        // The upright dome meshes in ~10 ms; the pole-crossing chart
        // took ~80-160 ms and the arm's fallback 44 s.
        assert!(
            elapsed.as_secs_f64() < 5.0,
            "{label}: a dome meshes in well under 5 s ({elapsed:?})"
        );
    }
}
