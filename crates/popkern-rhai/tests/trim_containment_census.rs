//! PK-792: `BRep::face_uv_trim_containment` agrees with sound
//! tessellation. Every display triangle's centroid, projected into its
//! face's chart, classifies Interior or Boundary on gallery bodies whose
//! faces were measured reading Outside: a rounded box's plane faces
//! (six chords per corner arc left the sliver between chord and arc
//! Outside) and a die's spherical pip, a polar cap whose seam pcurves
//! both sit on one branch (the cap's polygon collapsed to zero area).

use popkern_core::validation::FaceUvContainment;

fn outside_counts(script: &str) -> Vec<(String, usize, usize)> {
    let path = format!("{}/examples/{script}.rhai", env!("CARGO_MANIFEST_DIR"));
    let source = std::fs::read_to_string(&path).expect("script");
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let mut model = host
        .run(&source, script)
        .expect("runs")
        .model
        .expect("model");
    let parts = model.to_parts_all().expect("parts");
    let mut rows = Vec::new();
    for (pi, part) in parts.iter().enumerate() {
        let brep = part.brep();
        let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::display())
            .expect("mesh");
        for fm in &mesh.face_meshes {
            let surface = &brep.faces[fm.face]
                .surface
                .as_ref()
                .expect("surface")
                .geometry;
            let mut outside = 0;
            for t in &fm.triangles {
                let [a, b, c] = t.map(|i| fm.vertices[i].position);
                let centroid = popkern_core::geometry::Point3::new(
                    (a.x + b.x + c.x) / 3.0,
                    (a.y + b.y + c.y) / 3.0,
                    (a.z + b.z + c.z) / 3.0,
                );
                let Some((u, v)) =
                    popkern_core::geometry::uv_projection::project_point_to_uv(surface, &centroid)
                else {
                    continue;
                };
                if brep.face_uv_trim_containment(fm.face, surface, u, v)
                    == Some(FaceUvContainment::Outside)
                {
                    outside += 1;
                }
            }
            if outside > 0 {
                rows.push((
                    format!("{script}#{pi} {:?} {}", fm.face, surface.kind_name()),
                    outside,
                    fm.triangles.len(),
                ));
            }
        }
    }
    rows
}

#[test]
fn rounded_box_faces_contain_their_own_triangles() {
    let rows = outside_counts("hinged_box");
    assert!(rows.is_empty(), "triangles classified Outside: {rows:?}");
}

#[test]
fn a_die_pip_cap_contains_its_own_triangles() {
    let rows: Vec<_> = outside_counts("polyhedral_dice")
        .into_iter()
        .filter(|(name, _, _)| name.contains("Sphere"))
        .collect();
    assert!(
        rows.is_empty(),
        "sphere triangles classified Outside: {rows:?}"
    );
}
