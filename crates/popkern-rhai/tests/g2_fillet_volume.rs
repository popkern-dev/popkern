//! PK-1193: a body with G2 rounds measures its volume. Its corner caps
//! are Gregory patches, and a separate volume path for such bodies
//! signed each face's divergence term by the face's `orientation` flag
//! (a winding contract, not a material side): the G2 boss plate read
//! 1.0043 for a mesh volume of 1.7593. Every face is now integrated on
//! its exact surface, oriented by the measuring mesh.

use popkern_core::tolerance::Tolerance;

#[test]
fn the_g2_boss_plate_measures_its_mesh_volume() {
    let source = r#"let m = model();
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate_sketch").extrude(blind(0.4));
let boss = m.tool(|m| m.profile(xy_at(0.2)).circle(0.35).sketch("boss_sketch").extrude(blind(0.8)));
plate = plate.union(boss);
plate.fillet(
    edges().between("cylinder", "plane").at_height(0.4)
        .or(edges().straight().at_height(0.4))
        .or(edges().straight().direction("z")),
    0.08, #{ continuity: "g2", label: "rounds" });
m
"#;
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let b = host
        .run(source, "g2_boss_plate")
        .unwrap_or_else(|e| panic!("{e}"));
    let b = b.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone();
    let has_gregory = b.faces.values().any(|f| {
        f.surface.as_ref().is_some_and(|s| {
            matches!(
                s.geometry,
                popkern_core::geometry::surface::SurfaceGeometry::Gregory(_)
            )
        })
    });
    assert!(has_gregory, "the G2 corners are Gregory caps");
    // The fine mesh's divergence sum: short of the exact volume by its
    // facets' sagitta only.
    let options = popkern_mesh::TessellationOptions {
        chord_height: 1e-4,
        ..Default::default()
    };
    let mesh = popkern_mesh::tessellate(&b, &options).unwrap();
    let mesh_volume: f64 = mesh
        .face_meshes
        .iter()
        .flat_map(|fm| {
            fm.triangles.iter().map(move |t| {
                let [a, c, d] = t.map(|i| fm.vertices[i].position);
                (a.x * (c.y * d.z - c.z * d.y) - a.y * (c.x * d.z - c.z * d.x)
                    + a.z * (c.x * d.y - c.y * d.x))
                    / 6.0
            })
        })
        .sum::<f64>()
        .abs();
    let volume = popkern_measure::body_volume(&b, &Tolerance::default()).unwrap();
    assert!(
        volume >= mesh_volume && (volume - mesh_volume) / volume < 1e-4,
        "body_volume {volume} against the mesh's {mesh_volume}"
    );
}
