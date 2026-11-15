//! PK-1215: a NURBS ball offsets to the concentric ball exactly: its net
//! scaled about the centre (the fitted offset stopped at 8.6e-4 across
//! the net's C1 knot rows, its poles torn by a near-zero chart normal),
//! and its seam recomputed as the offset surface's own iso-curve.

use popkern_core::tolerance::Tolerance;

const BALL: &str = r#"fn ball_net(c, r, s) {
    let w = sqrt(0.5);
    let ci = [[1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [-1.0, 1.0], [-1.0, 0.0], [-1.0, -1.0], [0.0, -1.0], [1.0, -1.0], [1.0, 0.0]];
    let cw = [1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
    let mer = [[0.0, -r], [r, -r], [r, 0.0], [r, r], [0.0, r]];
    let mw = [1.0, w, 1.0, w, 1.0];
    let net = [];
    let wts = [];
    for i in 0..9 {
        let row = [];
        let wr = [];
        for j in 0..5 {
            row.push([c[0] + s[0] * mer[j][0] * ci[i][0], c[1] + s[1] * mer[j][0] * ci[i][1], c[2] + s[2] * mer[j][1]]);
            wr.push(cw[i] * mw[j]);
        }
        net.push(row);
        wts.push(wr);
    }
    #{ net: net, weights: wts }
}
fn ball(m, c, r, s) {
    let b = ball_net(c, r, s);
    m.sheet_nurbs(b.net, #{ degree: [2, 2], weights: b.weights,
        knots: [[0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0], [0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0]] })
}"#;

#[test]
fn a_nurbs_ball_grows_and_shrinks_to_concentric_balls() {
    for (d, r) in [(2.0_f64, 7.0_f64), (-1.5, 3.5)] {
        let source = format!(
            "{BALL}\nlet m = model();\nlet b = ball(m, [1.0, -2.0, 0.5], 5.0, [1.0, 1.0, 1.0]);\n\
             b = b.offset_solid({d});\nm"
        );
        let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
        let run = host
            .run(&source, "nurbs_ball")
            .unwrap_or_else(|e| panic!("{d}: {e}"));
        let body = run.model.expect("model").to_parts_all().expect("parts")[0]
            .brep()
            .clone();
        let errors = body.validate_all();
        assert!(errors.is_empty(), "{d}: {errors:?}");
        let v = popkern_measure::body_volume(&body, &Tolerance::default()).unwrap();
        let want = 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);
        assert!((v - want).abs() < 1e-6 * want, "{d}: {v} want {want}");
    }
}
