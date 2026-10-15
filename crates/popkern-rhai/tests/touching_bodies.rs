//! Bodies that only touch (PK-756): meeting on coincident faces whose
//! normals oppose, their interiors are disjoint, so their intersection
//! is empty and A − B is A; their union merges into one solid.
//!
//! The cycloid module's output pin stands on its flange: the
//! intersection served the pin's two coincident end discs as a 301 mm³
//! "overlap", and the interference check reported the pin colliding
//! with its own flange. A hub in a race's bore, their bottoms on one
//! plane, adds a cooriented coincidence that only abuts.

use popkern_core::topology::BRep;
use popkern_rhai::ScriptHost;

/// The result of `a <op> b` for two bodies given as script lines.
fn run(a: &str, b: &str, op: &str) -> Result<BRep, String> {
    let source = format!(
        "let m = model();\nlet a = {a};\nlet b = m.tool(|t| {b});\na.{op}(b, #{{ label: \"op\" }});\nm"
    );
    let mut host = ScriptHost::new();
    let mut outcome = host.run(&source, "touch.rhai").map_err(|e| e.render())?;
    let model = outcome.model.as_mut().ok_or("no model")?;
    model
        .to_part()
        .map(|p| p.brep().clone())
        .map_err(|e| e.to_string())
}

fn volume(b: &BRep) -> f64 {
    popkern_measure::body_volume(b, &b.tolerance).expect("volume")
}

/// (body A, body B, A's volume, B's volume), each pair touching only.
fn cases() -> Vec<(&'static str, &'static str, f64, f64)> {
    let pi = std::f64::consts::PI;
    vec![
        // A pin (r 1.8, z 44..89) standing on a plate whose underside is
        // z = 89.
        (
            "m.profile(xy_at(89.0)).rect(80.0, 80.0).sketch(\"plate\").extrude(blind(10.0))",
            "t.profile(xy_at(44.0)).circle(1.8).at(16.0, 0.0).sketch(\"pin\").extrude(blind(45.0))",
            6400.0 * 10.0,
            pi * 1.8 * 1.8 * 45.0,
        ),
        // Two boxes sharing part of a face.
        (
            "m.profile(xy()).rect(10.0, 10.0).sketch(\"lo\").extrude(blind(5.0))",
            "t.profile(xy_at(5.0)).rect(6.0, 6.0).at(3.0, 1.0).sketch(\"hi\").extrude(blind(4.0))",
            500.0,
            144.0,
        ),
        // A block standing on the floor of a pocket (z = 2), clear of
        // its walls.
        (
            "m.profile(xy()).rect(20.0, 20.0).sketch(\"cup\").extrude(blind(10.0))\
             .subtract(m.tool(|t| t.profile(xy_at(2.0)).rect(16.0, 16.0).sketch(\"pocket\").extrude(blind(9.0))))",
            "t.profile(xy_at(2.0)).rect(6.0, 6.0).sketch(\"block\").extrude(blind(4.0))",
            4000.0 - 16.0 * 16.0 * 8.0,
            144.0,
        ),
        // A hub (r 27..29) in a bearing race's bore (r 29..34), their
        // bottom annuli on one plane: a cooriented coincidence that only
        // abuts at r 29 (the cycloid M module's flange and inner race,
        // served as a 16587 mm³ two-face "overlap" before the contact
        // test set abutting pairs aside).
        (
            "m.profile(xz()).path(27.0, 0.0).line_to(29.0, 0.0).line_to(29.0, 20.0).line_to(27.0, 20.0).close()\
             .sketch(\"hub\").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())",
            "t.profile(xz()).path(29.0, 0.0).line_to(34.0, 0.0).line_to(34.0, 9.0).line_to(29.0, 9.0).close()\
             .sketch(\"race\").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())",
            pi * (29.0 * 29.0 - 27.0 * 27.0) * 20.0,
            pi * (34.0 * 34.0 - 29.0 * 29.0) * 9.0,
        ),
    ]
}

#[test]
fn touching_bodies_intersect_to_nothing() {
    for (a, b, _, _) in cases() {
        match run(a, b, "intersect") {
            // An empty body is refused as a part: that is the empty result.
            Err(e) => assert!(e.contains("empty body"), "{a} ∩ {b}: {e}"),
            Ok(brep) => panic!(
                "{a} ∩ {b}: served {} faces, volume {}",
                brep.faces.len(),
                volume(&brep)
            ),
        }
    }
}

#[test]
fn subtracting_a_touching_body_leaves_the_target_whole() {
    for (a, b, va, _) in cases() {
        let brep = run(a, b, "subtract").unwrap_or_else(|e| panic!("{a} - {b}: {e}"));
        assert!(brep.validate().is_ok());
        assert!(
            ((volume(&brep) - va) / va).abs() < 1e-9,
            "{a} - {b}: volume {}",
            volume(&brep)
        );
    }
}

#[test]
fn touching_bodies_unite_into_one_solid() {
    for (a, b, va, vb) in cases() {
        let brep = run(a, b, "union").unwrap_or_else(|e| panic!("{a} + {b}: {e}"));
        assert!(brep.validate().is_ok());
        assert_eq!(
            brep.solids.len(),
            1,
            "{a} + {b}: {} solids",
            brep.solids.len()
        );
        let v = volume(&brep);
        assert!(
            ((v - va - vb) / (va + vb)).abs() < 1e-3,
            "{a} + {b}: volume {v}"
        );
    }
}
