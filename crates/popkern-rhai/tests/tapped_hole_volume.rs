//! A tapped hole's measured volume is smooth in its fit clearance
//! (PK-901). body_volume integrates the fine tessellation's triangles on
//! the exact surface; on an M6 thread flank at a 0.08-0.1 fit, slivers
//! chording 27° of the helix tipped their facet normals across the
//! chart's, 50477 triangles fell back to facets whose errors no longer
//! cancelled, and a tapped block read 0.7% over its true volume (an arm
//! tool flange 73064 against OpenCASCADE's 56743).

use popkern_core::topology::BRep;

fn tapped_block(clearance: f64) -> BRep {
    let src = format!(
        r#"
let m = model();
let b = m.profile(xy()).rect(40.0, 40.0).sketch("block").extrude(blind(15.0));
b = b.tapped_hole(axis_at([0.0, 0.0, 15.0], [0.0, 0.0, -1.0]), "M6", blind(10.0),
    #{{ clearance: {clearance}, label: "tap" }});
m
"#
    );
    popkern_rhai::ScriptHost::new()
        .run(&src, "tapped_block.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()))
        .model
        .expect("model")
        .to_brep()
        .expect("brep")
}

fn mesh_volume(brep: &BRep) -> f64 {
    let options = popkern_mesh::TessellationOptions {
        chord_height: 0.002,
        ..Default::default()
    };
    let mesh = popkern_mesh::tessellate(brep, &options).expect("mesh");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.vertices[i].position);
        v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    v6 / 6.0
}

#[test]
fn a_tapped_holes_volume_is_smooth_in_its_fit_clearance() {
    let mut volumes = Vec::new();
    for clearance in [0.05, 0.08, 0.1, 0.12, 0.15] {
        let body = tapped_block(clearance);
        let measured = popkern_measure::body_volume(&body, &body.tolerance).expect("volume");
        let meshed = mesh_volume(&body);
        // The quadrature is the exact surface; the mesh sits a chord's
        // sagitta off it.
        assert!(
            ((measured - meshed) / meshed).abs() < 1e-4,
            "clearance {clearance}: body_volume {measured} vs mesh {meshed}"
        );
        volumes.push(measured);
    }
    // A wider fit removes more material, monotonically.
    assert!(
        volumes.windows(2).all(|w| w[1] < w[0]),
        "not monotone in clearance: {volumes:?}"
    );
}
