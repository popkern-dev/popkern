//! Ellipses, conics and splines in profiles (brep-bench §4.6): each
//! body serves its closed-form volume, and its emitted script replays
//! to the same volume.

use popkern_rhai::{ScriptHost, to_rhai};
use std::f64::consts::{FRAC_1_SQRT_2, PI};

fn volume(model: &mut popkern_script::Model) -> f64 {
    model
        .to_parts_all()
        .expect("parts")
        .iter()
        .map(|p| popkern_measure::body_volume(p.brep(), &p.brep().tolerance).expect("volume"))
        .sum()
}

fn run(source: &str, name: &str) -> popkern_script::Model {
    let mut host = ScriptHost::with_options(Default::default());
    host.run(source, name)
        .unwrap_or_else(|e| panic!("{name}: {e}\n{source}"))
        .model
        .expect("model")
}

/// The volume `body` (a Rhai expression on `m`) serves, checked
/// against `exact` and against its emitted script.
fn assert_volume(name: &str, body: &str, exact: f64) {
    let source = format!("let m = model();\nlet b = {body};\nm\n");
    let mut model = run(&source, name);
    let v = volume(&mut model);
    assert!(
        (v - exact).abs() <= 1e-6 * exact,
        "{name}: {v} against {exact}"
    );
    let emitted = to_rhai(&mut model).unwrap_or_else(|e| panic!("{name}: emit: {e}"));
    let v2 = volume(&mut run(&emitted, name));
    assert!(
        (v - v2).abs() <= 1e-9 * v,
        "{name}: emitted script serves {v2}, the script {v}\n{emitted}"
    );
}

#[test]
fn an_elliptic_cylinder() {
    assert_volume(
        "elliptic cylinder",
        "m.profile(xy()).ellipse(5.0, 3.0).at(1.0, 2.0).sketch(\"oval\").extrude(blind(10.0))",
        150.0 * PI,
    );
}

#[test]
fn an_ellipsoid_revolved_from_two_conic_quarters() {
    // The half ellipse x²/16 + z²/9 = 1, x ≥ 0, closed along the axis.
    let w = FRAC_1_SQRT_2;
    assert_volume(
        "ellipsoid",
        &format!(
            "m.profile(xz()).path(0.0, -3.0).conic_to([4.0, -3.0], [4.0, 0.0], {w:?}).conic_to([4.0, 3.0], [0.0, 3.0], {w:?}).close().sketch(\"half\").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())"
        ),
        4.0 / 3.0 * PI * 16.0 * 3.0,
    );
}

#[test]
fn a_paraboloid_cup_from_a_parabola() {
    // z = 6 (x / 2)² from the vertex to the rim at x = 2, closed across
    // the top and down the axis: the solid above the paraboloid, π r² h / 2.
    assert_volume(
        "paraboloid",
        "m.profile(xz()).path(0.0, 0.0).conic_to([1.0, 0.0], [2.0, 6.0], 1.0).line_to(0.0, 6.0).close().sketch(\"cup\").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())",
        PI * 4.0 * 6.0 / 2.0,
    );
}

#[test]
fn a_closed_quadratic_spline_over_a_square() {
    // Uniform quadratic over the square (±s, ±s): four parabolic arcs
    // between the edge midpoints, the corners their shoulders. The
    // midpoint diamond (2 s²) plus four parabolic segments of 2/3 of
    // their s²/2 triangles: (10/3) s².
    let s = 3.0;
    assert_volume(
        "square spline",
        "m.profile(xy()).spline([[3.0, 3.0], [-3.0, 3.0], [-3.0, -3.0], [3.0, -3.0]], #{ degree: 2 }).sketch(\"pillow\").extrude(blind(4.0))",
        10.0 / 3.0 * s * s * 4.0,
    );
}

#[test]
fn a_cubic_closes_a_slot() {
    // A cubic Bézier across the top from (-2, 0) to (2, 0), a line down,
    // a half circle of radius 2 round the bottom, and a line back up.
    let cps = [(-2.0, 0.0), (-1.0, 2.0), (1.0, -1.0), (2.0, 0.0)];
    let bez = |t: f64| {
        let u = 1.0 - t;
        let b = [u * u * u, 3.0 * u * u * t, 3.0 * u * t * t, t * t * t];
        (
            (0..4).map(|i| b[i] * cps[i].0).sum::<f64>(),
            (0..4).map(|i| b[i] * cps[i].1).sum::<f64>(),
        )
    };
    // ∫ y dx along the cubic (its signed area over the chord y = 0).
    let n = 200_000;
    let above: f64 = (0..n)
        .map(|i| {
            let (a, b) = (bez(i as f64 / n as f64), bez((i + 1) as f64 / n as f64));
            (b.0 - a.0) * (a.1 + b.1) / 2.0
        })
        .sum();
    let area = 2.0 * PI + 8.0 + above;
    assert_volume(
        "slot",
        "m.profile(xy()).path(-2.0, 0.0).spline_to([[-1.0, 2.0], [1.0, -1.0], [2.0, 0.0]]).line_to(2.0, -2.0).arc_to([0.0, -4.0], [-2.0, -2.0]).close().sketch(\"slot\").extrude(blind(3.0))",
        area * 3.0,
    );
}

#[test]
fn a_bad_spline_is_named() {
    let mut host = ScriptHost::with_options(Default::default());
    let source = "let m = model();\nlet b = m.profile(xy()).path(0.0, 0.0).spline_to([[1.0, 1.0], [2.0, 0.0]], #{ degree: 5 }).close().sketch(\"bad\").extrude(blind(1.0));\nm\n";
    let err = match host.run(source, "bad") {
        Ok(_) => panic!("degree 5 over three control points"),
        Err(e) => e.to_string(),
    };
    assert!(err.contains("degree 5"), "{err}");
}

fn area(body: &str) -> f64 {
    let source = format!("let m = model();\nlet b = {body};\nm\n");
    run(&source, "area")
        .to_parts_all()
        .expect("parts")
        .iter()
        .map(|p| popkern_measure::body_surface_area(p.brep(), &p.brep().tolerance).expect("area"))
        .sum()
}

/// Curved NURBS faces measure their exact area (the facets of their
/// mesh fell 4.5e-6 to 5.8e-6 short), against elementary closed forms.
#[test]
fn curved_faces_measure_their_exact_area() {
    // The paraboloid x² + y² = 8z to z = 8: side (32π/3)(5√5 − 1), top 64π.
    let got = area(
        "m.profile(xz()).path(0.0, 0.0).conic_to([4.0, 0.0], [8.0, 8.0], 1.0).line_to(0.0, 8.0).close().sketch(\"cup\").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())",
    );
    let exact = 32.0 * PI / 3.0 * (5.0 * 5f64.sqrt() - 1.0) + 64.0 * PI;
    assert!(
        (got - exact).abs() <= 1e-7 * exact,
        "paraboloid {got} vs {exact}"
    );
    // The square spline (s = 3) prism 4 tall: two ends (10/3) s² and
    // a side of four parabolic arcs, each s (1 + ln(1 + √2)/√2) long.
    let got = area(
        "m.profile(xy()).spline([[3.0, 3.0], [-3.0, 3.0], [-3.0, -3.0], [3.0, -3.0]], #{ degree: 2 }).sketch(\"pillow\").extrude(blind(4.0))",
    );
    let arc = 3.0 * (1.0 + (1.0 + 2f64.sqrt()).ln() / 2f64.sqrt());
    let exact = 2.0 * 30.0 + 4.0 * 4.0 * arc;
    assert!(
        (got - exact).abs() <= 1e-7 * exact,
        "prism {got} vs {exact}"
    );
}
