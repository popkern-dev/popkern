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

/// PK-819: the bored tube's annular caps wind their bore holes opposite
/// their outer loops (the start cap's hole ran with its outer loop; the
/// canonical-forward repair at the boolean's exit convicts the lying
/// hole pcurve), and the tube validates.
#[test]
fn a_bored_tube_cap_winds_its_hole_opposite_its_rim() {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let mut model = host
        .run(
            r#"
let m = model();
let x1 = -505.0;
let body = m.profile(plane_at([x1, -50.0, 0.0], [-1.0, 0.0, 0.0])).circle(33.0).sketch("tube")
    .extrude(blind(235.0));
body = body.subtract(m.tool(|t| t.profile(plane_at([x1 + 1.0, -50.0, 0.0], [-1.0, 0.0, 0.0]))
    .circle(29.0).sketch("bore").extrude(blind(237.0))), #{ label: "bore" });
m
"#,
            "tube",
        )
        .expect("runs")
        .model
        .expect("model");
    let parts = model.to_parts_all().expect("parts");
    let brep = parts[0].brep();
    assert!(brep.validate_all().is_empty(), "{:?}", brep.validate_all());
    let signed_area = |l: popkern_core::topology::LoopHandle| {
        let start = brep.loops[l].half_edge;
        let (mut cur, mut pts) = (start, Vec::new());
        loop {
            let he = &brep.half_edges[cur];
            let pc = he.pcurve.as_ref().expect("pcurve");
            for k in 0..16 {
                let q = pc.point_at(pc.t_start + (pc.t_end - pc.t_start) * f64::from(k) / 16.0);
                pts.push((q.x, q.y));
            }
            cur = he.next;
            if cur == start {
                break;
            }
        }
        let n = pts.len();
        (0..n)
            .map(|i| pts[i].0 * pts[(i + 1) % n].1 - pts[(i + 1) % n].0 * pts[i].1)
            .sum::<f64>()
    };
    let mut caps = 0;
    for (_, face) in &brep.faces {
        for &hole in &face.inner_loops {
            caps += 1;
            assert!(
                signed_area(face.outer_loop) * signed_area(hole) < 0.0,
                "a cap hole winds with its outer loop"
            );
        }
    }
    assert_eq!(caps, 2, "two annular caps");
}
