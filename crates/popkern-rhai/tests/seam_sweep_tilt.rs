//! Display facets stay within the tilt bound wherever a periodic
//! chart's seam falls.
//!
//! The tilt finish (bisection and vertex relaxation) planned in raw
//! chart coordinates, but a seam vertex is stored once, at one end of
//! the face's window: a facet across the seam mixed corners a period
//! apart, its candidate midpoints landed on the far side of the tube,
//! and every repair was refused. The six-axis arm's M4 screws showed
//! it only in their world copies (one facet at 15.0 deg): placing a
//! screw re-anchors its cylinders' u, which moved the seam. Turning a
//! screw about its own axis sweeps the seam through every position.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, SurfaceEval, Vec3};
use popkern_rhai::ScriptHost;

/// Facets whose normal leans past the 12 deg bound from the surface
/// normal at one of their vertices (chart-singular vertices excepted),
/// and the worst lean.
fn leaning(brep: &popkern_core::topology::BRep, mesh: &popkern_mesh::TriMesh) -> (usize, f64) {
    let (mut over, mut worst) = (0usize, 0.0f64);
    for fm in &mesh.face_meshes {
        let Some(surface) = brep.faces[fm.face].surface.as_ref().map(|s| &s.geometry) else {
            continue;
        };
        for t in &fm.triangles {
            let p: [Point3; 3] = t.map(|i| fm.vertices[i].position);
            let n: Vec3 = (p[1] - p[0]).cross(&(p[2] - p[0]));
            if n.length() < 1e-300 {
                continue;
            }
            let mut tilt = 0.0f64;
            let mut singular = false;
            for &i in t {
                let v = &fm.vertices[i];
                let (su, sv) = (
                    surface.partial_u(v.uv.0, v.uv.1),
                    surface.partial_v(v.uv.0, v.uv.1),
                );
                if su.cross(&sv).length() < 1e-9 * (su.length() * sv.length()).max(1e-300) {
                    singular = true;
                    break;
                }
                let c = n.dot(&v.normal) / (n.length() * v.normal.length());
                tilt = tilt.max(c.clamp(-1.0, 1.0).acos().to_degrees());
            }
            if singular {
                continue;
            }
            over += usize::from(tilt > 12.0);
            worst = worst.max(tilt);
        }
    }
    (over, worst)
}

#[test]
fn a_screw_turned_about_its_axis_keeps_every_facet_within_the_bound() {
    let source = r#"let m = model();
m.bolt("M4", 12.0, #{ head: "socket", label: "bolt" });
m
"#;
    let mut run = ScriptHost::new()
        .run(source, "screw.rhai")
        .expect("screw builds");
    let mut model = run.model.take().expect("model");
    let parts = model.to_parts_all().expect("parts");
    let screw = parts[0].brep().clone();
    let mut failures = Vec::new();
    for k in 0..12 {
        let turn = (k as f64 + 0.37) * std::f64::consts::TAU / 12.0;
        let turned = screw
            .transform(&Transform3::rotation_z(turn))
            .expect("turn");
        let mesh = popkern_mesh::tessellate(&turned, &popkern_mesh::TessellationOptions::display())
            .expect("mesh");
        let (over, worst) = leaning(&turned, &mesh);
        if over > 0 {
            failures.push(format!(
                "turn {turn:.3}: {over} facets lean, worst {worst:.1} deg"
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
