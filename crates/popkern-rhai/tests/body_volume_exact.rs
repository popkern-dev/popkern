//! PK-807: `popkern_measure::body_volume` is exact on scripted bodies.
//! A tube bored through a sketched rod carries annular caps whose bore
//! rims wind LIKE their outer rims; the planar term must subtract each
//! hole by magnitude (it added the bore and the six-axis arm's forearm
//! tube read 3.9x its volume).

fn volume_of(source: &str) -> f64 {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let mut model = host
        .run(source, "volume")
        .expect("runs")
        .model
        .expect("model");
    let parts = model.to_parts_all().expect("parts");
    parts
        .iter()
        .map(|p| popkern_measure::body_volume(p.brep(), &p.brep().tolerance).expect("volume"))
        .sum()
}

#[test]
fn a_bored_through_tube_reads_its_analytic_volume() {
    // The six-axis arm's forearm tube, as its script builds it: the tube
    // and the bore lengths are parameter expressions, and the profile
    // plane faces -x at x = -ua - 80 (ua = 425).
    let v = volume_of(
        r#"
let m = model();
m.param("forearm", 392.0, "m");
m.param("tube_ends", 157.0, "m");
m.param("bore_ends", 155.0, "m");
let x1 = -505.0;
let body = m.profile(plane_at([x1, -50.0, 0.0], [-1.0, 0.0, 0.0])).circle(33.0).sketch("tube")
    .extrude(blind(expr("forearm - tube_ends")));
body = body.subtract(m.tool(|t| t.profile(plane_at([x1 + 1.0, -50.0, 0.0], [-1.0, 0.0, 0.0]))
    .circle(29.0).sketch("bore").extrude(blind(expr("forearm - bore_ends")))), #{ label: "bore" });
m
"#,
    );
    let exact = std::f64::consts::PI * (33.0f64.powi(2) - 29.0f64.powi(2)) * 235.0;
    assert!((v - exact).abs() < 1e-9 * exact, "{v} vs {exact}");
}
