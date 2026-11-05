//! Thickening sheets with circular free rims (brep-bench thicken/tube,
//! half-pipe, hemisphere, spherical cap): the wall between a rim and
//! its offset is an exact annulus, cylinder or cone. The planar quad
//! the wall used to be fitted to collapsed on a closed rim (its two
//! ends one vertex), and every such thicken declined.

use popkern_rhai::ScriptHost;
use std::f64::consts::PI;

fn volume(source: &str) -> f64 {
    let mut host = ScriptHost::with_options(Default::default());
    let mut model = host
        .run(source, "thicken")
        .unwrap_or_else(|e| panic!("{e}"))
        .model
        .expect("model");
    let parts = model.to_parts_all().expect("parts");
    parts
        .iter()
        .map(|p| popkern_measure::body_volume(p.brep(), &p.brep().tolerance).expect("volume"))
        .sum()
}

/// A tube: the side of a cylinder of radius 5, 10 tall, thickened 1
/// outward: π (6² − 5²) 10.
#[test]
fn a_tube_thickens_outward() {
    let v = volume(
        "let m = model();\nlet c = m.cylinder(5.0, 10.0);\nc = c.remove_faces(faces().kind(\"plane\"));\nc = c.thicken(1.0, #{ side: \"outward\" });\nm\n",
    );
    let exact = PI * 11.0 * 10.0;
    assert!((v - exact).abs() < 1e-9 * exact, "{v} vs {exact}");
}

/// A hemisphere of radius 5 thickened 1 outward: (2/3) π (6³ − 5³).
#[test]
fn a_hemisphere_thickens_outward() {
    let v = volume(
        "let m = model();\nlet s = m.sphere(5.0, #{ lat: [0.0, 90.0] });\ns = s.remove_faces(faces().kind(\"plane\"));\ns = s.thicken(1.0, #{ side: \"outward\" });\nm\n",
    );
    let exact = 2.0 / 3.0 * PI * (216.0 - 125.0);
    assert!((v - exact).abs() < 1e-7 * exact, "{v} vs {exact}");
}

/// The tube inward and both ways: π (5² − 4²) 10 and π (5.5² − 4.5²) 10.
#[test]
fn a_tube_thickens_inward_and_both_ways() {
    for (side, exact) in [("inward", PI * 9.0 * 10.0), ("both", PI * 10.0 * 10.0)] {
        let v = volume(&format!(
            "let m = model();\nlet c = m.cylinder(5.0, 10.0);\nc = c.remove_faces(faces().kind(\"plane\"));\nc = c.thicken(1.0, #{{ side: \"{side}\" }});\nm\n"
        ));
        assert!((v - exact).abs() < 1e-9 * exact, "{side}: {v} vs {exact}");
    }
}

/// A cap of a ball of radius 5 above latitude 30° (its rim 60° from the
/// pole), thickened 1 outward: its rim wall is a cone through the
/// centre; (2π/3)(6³ − 5³)(1 − cos 60°).
#[test]
fn a_spherical_cap_thickens_outward() {
    let v = volume(
        "let m = model();\nlet s = m.sphere(5.0, #{ lat: [30.0, 90.0] });\ns = s.remove_faces(faces().kind(\"plane\"));\ns = s.thicken(1.0, #{ side: \"outward\" });\nm\n",
    );
    let exact = 2.0 * PI / 3.0 * 91.0 * 0.5;
    assert!((v - exact).abs() < 1e-7 * exact, "{v} vs {exact}");
}

#[test]
fn a_hemisphere_thickens_inward() {
    let v = volume(
        "let m = model();\nlet s = m.sphere(5.0, #{ lat: [0.0, 90.0] });\ns = s.remove_faces(faces().kind(\"plane\"));\ns = s.thicken(1.0, #{ side: \"inward\" });\nm\n",
    );
    let exact = 2.0 / 3.0 * PI * (125.0 - 64.0);
    assert!((v - exact).abs() < 1e-7 * exact, "{v} vs {exact}");
}

/// brep-bench thicken/torus-inner-half: the inner half of a torus
/// (tube 2 about the middle circle 6) thickened 1 outward, away from
/// the middle circle: the sheet, the inner half of the tube-3 torus and
/// two cylinders of radius 6, 1 tall. The side came from a wedge
/// integral, meaningless on an open saddle, and it grew toward the
/// middle circle (148.33 for 216.50).
#[test]
fn an_inner_half_torus_thickens_away_from_its_middle_circle() {
    let v = volume(
        "let m = model();\nlet t = m.profile(xz()).path(6.0, -2.0).arc_to([4.0, 0.0], [6.0, 2.0]).line_to(6.0, -2.0).close().sketch(\"h\").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());\nt = t.remove_faces(faces().kind(\"cylinder\"));\nt = t.thicken(1.0, #{ side: \"outward\" });\nm\n",
    );
    let exact = 216.501118142;
    assert!((v - exact).abs() < 1e-8 * exact, "{v} vs {exact}");
}
