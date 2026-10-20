//! A tapped hole's thread faces mesh within their chord height at
//! every certification level (PK-900). A tap flank is a NURBS chart
//! whose iso-lines wind six turns round the hole: probed at fixed
//! fractions, one of an M6's three flanks returned to its own chord by
//! phase. The chart took the ruled Delaunay metric and its interior
//! lattice stayed at nine points, so its facets sagged 0.20 at a 0.011
//! chord and 2.8e-3 at 1.75e-4, the arm's gripper screws certified only
//! on the finest whole-flange meshes (23 s a screw), and a 0.05 gap
//! read as crossing meshes two levels deep.

use popkern_core::geometry::SurfaceGeometry;
use popkern_core::topology::BRep;

fn tapped_block() -> BRep {
    let src = r#"
let m = model();
let b = m.profile(xy()).rect(40.0, 40.0).sketch("block").extrude(blind(15.0));
b = b.tapped_hole(axis_at([0.0, 0.0, 15.0], [0.0, 0.0, -1.0]), "M6", blind(10.0),
    #{ clearance: 0.1, label: "tap" });
m
"#;
    popkern_rhai::ScriptHost::new()
        .run(src, "tapped_block.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()))
        .model
        .expect("model")
        .to_brep()
        .expect("brep")
}

#[test]
fn every_tap_flank_meshes_within_its_chord_at_each_certification_level() {
    let body = tapped_block();
    let base = popkern_mesh::TessellationOptions::default();
    // The chords and subdivision caps of the interference certifier's
    // three refinement levels on the arm's tool flange.
    for (level, chord) in [(1u32, 1.12e-2), (2, 1.4e-3), (3, 1.75e-4)] {
        let options = popkern_mesh::TessellationOptions {
            chord_height: chord,
            max_subdivisions: base.max_subdivisions << (3 * level),
            ..base
        };
        let mesh = popkern_mesh::tessellate(&body, &options).expect("mesh");
        let mut flanks = 0;
        for fm in &mesh.face_meshes {
            let nurbs = matches!(
                body.faces[fm.face].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Nurbs(_))
            );
            if !nurbs {
                continue;
            }
            flanks += 1;
            // The certifier's own facet bound: centroid and edge
            // midpoints against the surface.
            let worst = fm
                .triangles
                .iter()
                .map(|t| {
                    let corners = t.map(|i| fm.vertices[i].position);
                    let uvs = t.map(|i| fm.vertices[i].uv);
                    popkern_assembly::mesh_bvh::deviation(&body, fm.face, &corners, &uvs)
                })
                .fold(0.0_f64, f64::max);
            assert!(
                worst <= 1.5 * chord,
                "level {level}: flank {:?} sags {worst:.3e} at chord {chord:.2e} ({} triangles)",
                fm.face,
                fm.triangles.len()
            );
        }
        assert_eq!(flanks, 3, "level {level}: the tap's three flank faces");
    }
}
