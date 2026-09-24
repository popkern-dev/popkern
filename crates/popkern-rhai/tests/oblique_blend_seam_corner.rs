//! Oblique pair corners between equal-radius fillets serve with exact
//! seams (PK-522).
//!
//! The guide's triangular prism fillets its three vertical edges and
//! the front top edge at r = 0.08. The two top corners of the front
//! edge join a vertical strip and the front-top strip at 58° (the
//! prism's base angle), with the third edge sharp. The strips' terminal circles
//! on the top face do not coincide, so the near anchor lies on one
//! strip only, which is the stub topology. The bisector-ellipse path
//! declined these corners `BranchUnresolved`, and the unresolved-corner
//! gate (PK-516) then declined the whole fillet. Dispatch now routes
//! by closure topology, and the stub-closing march builds the seam.

use popkern_rhai::ScriptHost;

const GUIDE_PRISM: &str = r#"
let m = model();
let prism = m.profile(xy()).polygon([[0.0, 0.0], [1.0, 0.0], [0.5, 0.8]]).sketch("tri").extrude(blind(1.0));
m.label_edge("front_top", edges().straight().at_height(1.0).direction("x").chord_length_between(0.9, 1.1));
prism.fillet(edges().straight().direction("z").or(edges().labeled("front_top")), 0.08, #{ label: "rounds" });
m
"#;

#[test]
fn guide_prism_oblique_corners_serve_check9_clean() {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(GUIDE_PRISM, "prism.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
    let mut model = outcome.model.expect("model");
    let part = model.to_part().expect("part");
    let brep = part.brep();
    let errors = brep.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
    // 5 planes + 4 strips; the two front-top corners resolved as
    // strip×strip seams (no phantom flap faces left on the top).
    assert_eq!(brep.faces.len(), 9, "face count");
}
