//! Scratch: flagship STEP round-trip face census.
#[test]
#[ignore = "scratch"]
fn flagship_step_census() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples/output/flagship_boss_plate.step");
    let doc = popkern_io::import_step_file(&path, &popkern_io::StepImportOptions::default())
        .expect("import")
        .document;
    for product in &doc.products {
        for shape in &product.shapes {
            for solid in &shape.solids {
                let mut kinds = std::collections::BTreeMap::new();
                for (_, f) in solid.faces.iter() {
                    let k = match f.surface.as_ref().map(|s| &s.geometry) {
                        Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) => "Plane",
                        Some(popkern_core::geometry::SurfaceGeometry::Cylinder { .. }) => "Cyl",
                        Some(popkern_core::geometry::SurfaceGeometry::Torus { .. }) => "Torus",
                        Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) => "Nurbs",
                        _ => "?",
                    };
                    *kinds.entry(k).or_insert(0usize) += 1;
                }
                eprintln!(
                    "imported solid: faces={} edges={} kinds={kinds:?} validate={:?}",
                    solid.faces.len(),
                    solid.edges.len(),
                    solid.validate().err()
                );
            }
        }
    }
}
