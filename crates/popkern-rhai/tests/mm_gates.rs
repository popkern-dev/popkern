//! The watertight gate meshes a millimetre part at the part's own scale.
//! Its default chord height is an absolute 5 mm: a 16 mm die's four
//! 2.1 mm spherical pip dimples on one face meshed at the minimum
//! subdivisions and the welded mesh opened edges, though the same die in metres (and the
//! page's display mesh, with its facet-angle bound) was closed.

use popkern_rhai::ScriptHost;

fn four_pip_face(unit: &str) -> Result<(), String> {
    let source = format!(
        "let m = model(#{{ unit: \"{unit}\" }});\n\
         let die = m.profile(xy_at(-8.0)).rect(16.0, 16.0).sketch(\"block\").extrude(blind(16.0));\n\
         // The four pips of the die's 4 face (x = -8), as the dice set lays them.\n\
         let k = 0;\n\
         for yz in [[-4.4, -4.4], [-4.4, 4.4], [4.4, -4.4], [4.4, 4.4]] {{\n\
             let pip = m.tool(|t| t.profile(xz()).path(0.0, -2.1).arc_to([2.1, 0.0], [0.0, 2.1]).line_to(0.0, -2.1)\n\
                 .close().sketch(\"ball_\" + k).revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())\n\
                 .translate(-9.5, yz[0], yz[1]));\n\
             die = die.subtract(pip, #{{ label: \"dimple_\" + k }});\n\
             k += 1;\n\
         }}\n\
         m.assert_watertight();\nm"
    );
    ScriptHost::new()
        .run(&source, "dimple.rhai")
        .map(|_| ())
        .map_err(|e| e.render())
}

#[test]
fn a_millimetre_pip_face_passes_the_watertight_gate() {
    for unit in ["m", "mm"] {
        four_pip_face(unit).unwrap_or_else(|e| panic!("{unit}: {e}"));
    }
}

/// The four-pip face meshes watertight at any chord (PK-839). Each
/// dimple's rim is two semicircular arc edges; at a chord coarser than
/// the 1.47 rim radius the edge sampler kept only their ends, the
/// planar face's holes collapsed to two-point rings and were dropped
/// (the face meshed as if solid), and 64 welded edges opened against
/// the dimples. A span turning past ~110 deg now splits at any chord.
#[test]
fn the_pip_face_meshes_watertight_at_a_coarse_chord() {
    use std::collections::HashMap;
    let source = "let m = model();\n\
        let die = m.profile(xy_at(-8.0)).rect(16.0, 16.0).sketch(\"block\").extrude(blind(16.0));\n\
        let k = 0;\n\
        for yz in [[-4.4, -4.4], [-4.4, 4.4], [4.4, -4.4], [4.4, 4.4]] {\n\
            let pip = m.tool(|t| t.profile(xz()).path(0.0, -2.1).arc_to([2.1, 0.0], [0.0, 2.1]).line_to(0.0, -2.1)\n\
                .close().sketch(\"ball_\" + k).revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())\n\
                .translate(-9.5, yz[0], yz[1]));\n\
            die = die.subtract(pip, #{ label: \"dimple_\" + k });\n\
            k += 1;\n\
        }\nm";
    let mut run = ScriptHost::new().run(source, "pips.rhai").expect("builds");
    let brep = run.model.take().expect("model").to_brep().expect("brep");
    for chord in [5.0, 50.0] {
        let options = popkern_mesh::TessellationOptions {
            chord_height: chord,
            ..Default::default()
        };
        let mesh = popkern_mesh::tessellate(&brep, &options).expect("mesh");
        let key =
            |p: &popkern_core::geometry::Point3| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
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
        let open = edges.values().filter(|&&n| n != 2).count();
        assert_eq!(open, 0, "chord {chord}: welded open or non-manifold edges");
    }
}
