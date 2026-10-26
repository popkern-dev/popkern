// A wall bracket carved by generative design (PK-956). The design space
// is a solid 50 x 16 x 36 mm block against the wall at x = 0: two M4
// screws pass through its 6 mm wall flange, each with room behind it
// for its head and driver, and an eye at the tip carries the load. The
// solver keeps a quarter of the material where it makes the stiffest
// bracket under 200 N down on the eye, printable flat on its side
// (+y) without supports; the kept faces (the wall, the bores, the eye)
// stay the design space's.
//
// The result is a body like any other: its faces are planar facets,
// the wall and eye keep their names, and the Simulate tab checks it.
// Measured there (wall held, 200 N down on the eye, PLA): the bracket
// moves 0.63 mm with a 21 MPa peak (inter-layer safety factor 1.3)
// where the solid block moves 0.074 mm with 5.1 MPa (5.6), at a
// quarter of the block's material.
// In the browser a design takes minutes to solve, so an edit that
// changes it asks before re-solving; this page's design was solved
// when the gallery was built.
let m = model();
m.param("reach", 50.0);
m.param("load", 200.0);

let space = m.profile(yz()).rect(16.0, 36.0).sketch("space").extrude(blind(p("reach")));
// The eye, across the bracket near its tip (at 0.88 of the reach).
let eye = m.tool(|m| m.profile(xz_at(-9.0)).circle(3.0).at(expr("reach * 0.88"), 0.0)
    .sketch("eye").extrude(blind(18.0)));
space = space.subtract(eye);
m.label_face("eye", faces().kind("cylinder"));
// Two screws through the flange; behind each, a pocket for the head.
for z in [12.0, -12.0] {
    let bore = m.tool(|m| m.profile(yz_at(-1.0)).circle(2.2).at(0.0, z).sketch("bore").extrude(blind(8.0)));
    space = space.subtract(bore);
    m.label_face(if z > 0.0 { "bolt_top" } else { "bolt_bottom" }, faces().kind("cylinder").at_height(z));
    let pocket = m.tool(|m| m.profile(yz_at(6.0)).circle(4.5).at(0.0, z).sketch("pocket")
        .extrude(blind(expr("reach * 0.92"))));
    space = space.subtract(pocket);
}
m.label_face("wall", faces().facing_neg("x"));

let bracket = space.generative(#{
    supports: [#{ faces: faces().labeled("wall"), kind: "fixed" }],
    loads: [#{ faces: faces().labeled("eye"), force: p("load"), direction: [0.0, 0.0, -1.0] }],
    keep: [faces().labeled("wall"), faces().labeled("eye"), faces().labeled("bolt_top"), faces().labeled("bolt_bottom")],
    volume_fraction: 0.25,
    filter: 3.75,
    keep_depth: 2.5,
    cell: 1.25,
    build_direction: "+y",
    max_iterations: 300,
    label: "carve",
});
m.assert_watertight();
m
