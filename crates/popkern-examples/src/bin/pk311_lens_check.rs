// scratch: does the micro-lens serve now?
fn main() {
    let tol = popkern_core::tolerance::Tolerance::default();
    for depth in [0.0005_f64, 0.002, 0.008, 0.03] {
        let d = 2.0 * (1.0 + 0.35) - depth;
        let ta = popkern_core::topology::BRep::torus(1.0, 0.35, tol).expect("a");
        let tb = popkern_core::topology::BRep::torus(1.0, 0.35, tol)
            .expect("b")
            .transform(&popkern_core::geometry::Transform3::translation(
                d, 0.0, 0.0,
            ))
            .expect("t");
        match popkern_bool::boolean(popkern_bool::BoolOperation::Union, ta, tb, tol, None) {
            Ok(b) => {
                let vol =
                    popkern_mesh::tessellate(&b, &popkern_mesh::TessellationOptions::default())
                        .ok()
                        .map(|mesh| {
                            let mut vol6 = 0.0_f64;
                            for tri in &mesh.triangles {
                                let a = mesh.vertices[tri[0]].position;
                                let bb = mesh.vertices[tri[1]].position;
                                let c = mesh.vertices[tri[2]].position;
                                vol6 += a.x * (bb.y * c.z - bb.z * c.y)
                                    - a.y * (bb.x * c.z - bb.z * c.x)
                                    + a.z * (bb.x * c.y - bb.y * c.x);
                            }
                            vol6.abs() / 6.0
                        });
                println!(
                    "depth {depth}: SERVES {} faces validate={:?} validate_all={} vol {vol:?}",
                    b.faces.len(),
                    b.validate().err().map(|e| e.to_string()),
                    b.validate_all().len()
                );
            }
            Err(e) => println!("depth {depth}: DECLINES {e}"),
        }
    }
}
