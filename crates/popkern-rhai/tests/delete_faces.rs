//! PK-1199: `b.delete_faces(faces)` takes blend faces out of a solid and
//! closes each gap with its neighbours run on until they meet. A round
//! or bevel along a straight edge or round a circle, convex or concave,
//! a countersink and a corner's rounds each give back the solid before
//! the blend; a face that is no blend declines.

use popkern_core::tolerance::Tolerance;
use std::f64::consts::PI;

fn volume(source: &str) -> f64 {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(source, "delete_faces")
        .unwrap_or_else(|e| panic!("{e}"));
    let b = run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone();
    let errors = b.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
    popkern_measure::body_volume(&b, &Tolerance::default()).unwrap()
}

fn close(v: f64, want: f64) {
    assert!((v - want).abs() < 1e-7 * want, "{v} want {want}");
}

#[test]
fn a_straight_round_deletes_back_to_the_block() {
    close(
        volume(
            "let m = model();\nlet b = m.box(20.0, 10.0, 10.0);\nb = b.fillet(edges().straight().direction(\"y\").at_vertex([20.0, 0.0, 10.0]), 2.0);\nb = b.delete_faces(faces().kind(\"cylinder\"));\nm\n",
        ),
        2000.0,
    );
}

#[test]
fn a_concave_root_round_deletes_back_to_the_sharp_boss() {
    close(
        volume(
            "let m = model();\nlet b = m.cylinder(8.0, 4.0, #{ at: [0.0, 0.0, -4.0] });\nb = b.union(m.tool(|m| m.cylinder(3.0, 8.0)));\nb = b.fillet(edges().between(\"cylinder\", \"plane\").at_height(0.0), 1.0);\nb = b.delete_faces(faces().kind(\"torus\"));\nm\n",
        ),
        PI * (64.0 * 4.0 + 9.0 * 8.0),
    );
}

#[test]
fn a_rim_round_and_a_countersink_delete() {
    close(
        volume(
            "let m = model();\nlet b = m.cylinder(6.0, 10.0);\nb = b.fillet(edges().at_height(10.0), 2.0);\nb = b.delete_faces(faces().kind(\"torus\"));\nm\n",
        ),
        PI * 36.0 * 10.0,
    );
    close(
        volume(
            "let m = model();\nlet b = m.box(20.0, 10.0, 10.0);\nb = b.subtract(m.tool(|m| m.cylinder(2.0, 12.0, #{ at: [10.0, 5.0, -1.0] })));\nb = b.subtract(m.tool(|m| m.cone(2.0, 4.0, 2.0, #{ at: [10.0, 5.0, 8.0] })));\nb = b.delete_faces(faces().kind(\"cone\"));\nm\n",
        ),
        2000.0 - PI * 4.0 * 10.0,
    );
}

#[test]
fn a_corner_of_rounds_deletes_back_to_the_cube() {
    close(
        volume(
            "let m = model();\nlet b = m.box(10.0, 10.0, 10.0);\nb = b.fillet(edges().at_vertex([10.0, 10.0, 10.0]), 2.0);\nb = b.delete_faces(faces().kind(\"cylinder\").or(faces().kind(\"sphere\")));\nm\n",
        ),
        1000.0,
    );
}

#[test]
fn a_face_that_is_no_blend_declines() {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let Err(e) = host.run(
            "let m = model();\nlet b = m.sphere(5.0);\nb = b.delete_faces(faces().kind(\"sphere\"));\nm\n",
            "delete_faces",
    ) else {
        panic!("a whole ball is no corner");
    };
    assert!(e.to_string().contains("corner ball"), "{e}");
}

#[test]
fn a_flat_face_moves_and_a_hole_widens() {
    // The top of a drilled block moved up 6: the sides and the hole run on.
    close(
        volume(
            "let m = model();\nlet b = m.box(40.0, 30.0, 10.0);\nb = b.subtract(m.tool(|m| m.cylinder(5.0, 12.0, #{ at: [20.0, 15.0, -1.0] })));\nb = b.offset_faces(faces().facing(\"z\").at_height(10.0), 6.0);\nm\n",
        ),
        (1200.0 - PI * 25.0) * 16.0,
    );
    // A hole's wall offset 2 into the plate: radius 5 to 7.
    close(
        volume(
            "let m = model();\nlet b = m.box(40.0, 30.0, 8.0);\nb = b.subtract(m.tool(|m| m.cylinder(5.0, 10.0, #{ at: [20.0, 15.0, -1.0] })));\nb = b.offset_faces(faces().kind(\"cylinder\"), -2.0);\nm\n",
        ),
        (1200.0 - PI * 49.0) * 8.0,
    );
}

#[test]
fn a_face_whose_neighbours_lean_onto_it_does_not_move() {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let Err(e) = host.run(
            "let m = model();\nlet b = m.box(10.0, 10.0, 10.0);\nb = b.chamfer(edges().at_height(10.0), 2.0);\nb = b.offset_faces(faces().facing(\"z\").at_height(10.0), 1.0);\nm\n",
            "offset_faces",
    ) else {
        panic!("the chamfers lean onto the top");
    };
    assert!(e.to_string().contains("leans onto the face"), "{e}");
}

/// PK-1220: a box rounded on every edge, every round and corner ball
/// deleted, is the box again: the kept faces' regions intersected.
#[test]
fn a_box_with_every_round_deleted_is_the_box() {
    close(
        volume(
            "let m = model();\nlet b = m.box(20.0, 16.0, 12.0);\nb = b.fillet(edges(), 2.0);\nb = b.delete_faces(faces().kind(\"cylinder\").or(faces().kind(\"sphere\")));\nm\n",
        ),
        3840.0,
    );
}

/// PK-1221: a block's top moved up 6 carries the round that blends it
/// into the back: deleted, the top moved, the new edge rounded again at
/// radius 3. Volume 40·30·16 − (1 − π/4)·3²·40 = 18840 + 90π.
#[test]
fn a_moved_top_carries_its_round() {
    close(
        volume(
            "let m = model();\nlet b = m.box(40.0, 30.0, 10.0);\nb = b.fillet(edges().straight().direction(\"x\").at_vertex([0.0, 30.0, 10.0]), 3.0);\nb = b.offset_faces(faces().facing(\"z\").at_height(10.0), 6.0);\nm\n",
        ),
        18840.0 + 90.0 * PI,
    );
}
