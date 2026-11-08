//! Sheet primitives (PK-1083): one planar face with free edges (the
//! canonical sheet form), the area its arguments name, replaying from
//! its emitted script; a fill bounds the face of a closed planar wire.

use popkern_core::topology::BodyKind;
use popkern_rhai::{ScriptHost, to_rhai};
use std::f64::consts::PI;

fn run(source: &str) -> Result<popkern_script::Model, String> {
    let mut host = ScriptHost::with_options(Default::default());
    host.run(source, "sheets")
        .map_err(|e| e.to_string())?
        .model
        .ok_or_else(|| "no model".into())
}

/// The total area of the model's sheets, each checked to be one face
/// whose every edge is free.
fn sheet_area(model: &mut popkern_script::Model) -> f64 {
    let mut area = 0.0;
    for body in model.alive_bodies() {
        let brep = model.brep_of(body).expect("body");
        assert_eq!(brep.body_kind(), BodyKind::Sheet);
        assert_eq!(brep.faces.len(), 1);
        assert_eq!(brep.free_edges().len(), brep.edges.len());
        area += popkern_measure::body_surface_area(&brep, &brep.tolerance).expect("area");
    }
    area
}

fn assert_area(body: &str, area: f64) {
    let source = format!("let m = model();\n{body}\nm\n");
    let mut model = run(&source).unwrap_or_else(|e| panic!("{e}\n{source}"));
    let a = sheet_area(&mut model);
    assert!(
        (a - area).abs() <= 1e-9 * area,
        "{body}: area {a} vs {area}"
    );
    let emitted = to_rhai(&mut model).expect("emit");
    let a2 = sheet_area(&mut run(&emitted).unwrap_or_else(|e| panic!("{e}\n{emitted}")));
    assert!((a - a2).abs() <= 1e-12 * a, "{emitted}");
}

#[test]
fn rectangles_discs_and_polygons() {
    assert_area("let s = m.sheet_rect(2.0, 3.0);", 6.0);
    assert_area(
        "let s = m.sheet_rect(2.0, 3.0, #{ centred: true, at: [0.0, 0.0, 5.0], axis: [0.0, 1.0, 0.0], ref: [1.0, 0.0, 0.0] });",
        6.0,
    );
    assert_area("let s = m.sheet_disc(2.0);", 4.0 * PI);
    assert_area(
        "let s = m.sheet_polygon([[0.0, 0.0], [4.0, 0.0], [0.0, 3.0]]);",
        6.0,
    );
}

/// A fill of a closed polyline, a circle and a half disc (an arc and
/// its chord as two wires' edges? one polyline and one arc: a wire per
/// edge, chained across bodies is not offered; one wire each).
#[test]
fn fills_of_closed_wires() {
    assert_area(
        "let w = m.wire_polyline([[0.0, 0.0, 1.0], [4.0, 0.0, 1.0], [4.0, 3.0, 1.0]], #{ closed: true });\nlet s = m.sheet_fill(w);",
        6.0,
    );
    assert_area(
        "let w = m.wire_circle(2.0, #{ at: [1.0, 2.0, 3.0], axis: [1.0, 1.0, 0.0], ref: [0.0, 0.0, 1.0] });\nlet s = m.sheet_fill(w);",
        4.0 * PI,
    );
    // An interpolated spline closed by nothing does not bound a face.
    let open = "let m = model();\nlet w = m.wire_line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);\nlet s = m.sheet_fill(w);\nm\n";
    match run(open) {
        Ok(_) => panic!("an open wire filled"),
        Err(e) => assert!(e.contains("do not close"), "{e}"),
    }
    let skew = "let m = model();\nlet w = m.wire_polyline([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 1.0]], #{ closed: true });\nlet s = m.sheet_fill(w);\nm\n";
    match run(skew) {
        Ok(_) => panic!("a skew wire filled"),
        Err(e) => assert!(e.contains("one plane"), "{e}"),
    }
}

/// A sheet thickens like one stripped from a solid: outward along the
/// face's normal (the placement's axis), inward against it, or both.
#[test]
fn a_sheet_thickens() {
    for (side, z) in [
        ("outward", (0.0, 2.0)),
        ("inward", (-2.0, 0.0)),
        ("both", (-1.0, 1.0)),
    ] {
        let source = format!(
            "let m = model();\nlet s = m.sheet_rect(10.0, 8.0);\ns = s.thicken(2.0, #{{ side: \"{side}\" }});\nm\n"
        );
        let mut model = run(&source).unwrap_or_else(|e| panic!("{side}: {e}"));
        let parts = model.to_parts_all().expect("parts");
        let brep = parts[0].brep();
        let v = popkern_measure::body_volume(brep, &brep.tolerance).expect("volume");
        assert!((v - 160.0).abs() < 1e-9 * 160.0, "{side}: {v}");
        let (lo, hi) = brep
            .vertices
            .values()
            .fold((f64::MAX, f64::MIN), |(lo, hi), v| {
                (lo.min(v.position.z), hi.max(v.position.z))
            });
        assert!(
            (lo - z.0).abs() < 1e-9 && (hi - z.1).abs() < 1e-9,
            "{side}: z {lo}..{hi}"
        );
    }
}
