//! PK-569: default (chord-only) meshes of the coil spring and propeller
//! honour the chord contract — every facet's centroid lies within the
//! chord height of the surface — and carry no near-antipodal vertex
//! normals (the task's 174-176 deg triangles).

use popkern_core::geometry::{Point3, SurfaceEval};

fn check(script: &str, body: usize) -> Vec<String> {
    let src = std::fs::read_to_string(format!(
        "{}/examples/{script}.rhai",
        env!("CARGO_MANIFEST_DIR")
    ))
    .expect("script");
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let mut model = host.run(&src, script).expect("runs").model.expect("model");
    let parts = model.to_parts_all().expect("parts");
    let brep = parts[body].brep();
    let options = popkern_mesh::TessellationOptions::default();
    let mesh = popkern_mesh::tessellate(brep, &options).expect("mesh");
    let angle = |a: popkern_core::geometry::Vec3, b: popkern_core::geometry::Vec3| {
        (a.dot(&b) / (a.length() * b.length()))
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    };
    let mut failures = Vec::new();
    for fm in &mesh.face_meshes {
        let surface = &brep.faces[fm.face]
            .surface
            .as_ref()
            .expect("surface")
            .geometry;
        for t in &fm.triangles {
            let n = t.map(|i| fm.vertices[i].normal);
            let spread = angle(n[0], n[1])
                .max(angle(n[1], n[2]))
                .max(angle(n[0], n[2]));
            if spread > 90.0 {
                failures.push(format!(
                    "{script}#{body}: vertex normals spread {spread:.1} deg"
                ));
            }
            let p = t.map(|i| fm.vertices[i].position);
            let centroid = Point3::new(
                (p[0].x + p[1].x + p[2].x) / 3.0,
                (p[0].y + p[1].y + p[2].y) / 3.0,
                (p[0].z + p[1].z + p[2].z) / 3.0,
            );
            let deviation = popkern_core::geometry::uv_projection::project_point_to_uv_multisheet(
                surface, &centroid,
            )
            .map(|(u, v)| surface.point_at(u, v).distance_to(&centroid))
            .unwrap_or(f64::INFINITY);
            if deviation > options.chord_height * (1.0 + 1e-9) {
                failures.push(format!(
                    "{script}#{body}: facet centroid {deviation:.3e} off the surface"
                ));
            }
        }
    }
    failures
}

#[test]
fn coil_and_propeller_default_meshes_keep_the_chord_contract() {
    let failures: Vec<String> = [
        ("coil_spring", 0),
        ("coil_spring", 1),
        ("coil_spring", 2),
        ("propeller", 0),
    ]
    .into_iter()
    .flat_map(|(s, b)| check(s, b))
    .collect();
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
