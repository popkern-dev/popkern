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
