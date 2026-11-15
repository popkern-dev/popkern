//! PK-1200: `m.sheet_nurbs` makes the sheet of a B-spline surface from
//! its control net, and a solid split by it keeps the exact surface.
//! A Bézier dome over a block (its boundary curves flush on the block's
//! sides) and the bilinear saddle over a drum read their closed forms;
//! a net with a collapsed edge declines (only open patches are sheets).
//! PK-1209: a NURBS ellipsoid cut square to its axis keeps the cap above
//! the cut at every height (it declined at most of them).

use popkern_core::tolerance::Tolerance;

fn volume(source: &str) -> f64 {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(source, "sheet_nurbs")
        .unwrap_or_else(|e| panic!("{e}"));
    let b = run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone();
    let errors = b.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
    popkern_measure::body_volume(&b, &Tolerance::default()).unwrap()
}

#[test]
fn a_block_under_a_bezier_dome_reads_its_closed_form() {
    // 400 · Σ z / 16: each cubic Bernstein polynomial integrates to 1/4.
    let v = volume(
        r#"let m = model();
let z = [[4.0, 6.0, 7.0, 5.0], [6.0, 14.0, 12.0, 7.0], [5.0, 13.0, 16.0, 8.0], [3.0, 6.0, 8.0, 6.0]];
let net = [];
for i in 0..4 { let row = []; for j in 0..4 { row.push([20.0 * i / 3.0, 20.0 * j / 3.0, z[i][j]]); } net.push(row); }
let dome = m.sheet_nurbs(net, #{ degree: 3 });
let b = m.box(20.0, 20.0, 30.0).split([dome], containing([10.0, 10.0, 1.0]))[0];
m"#,
    );
    assert!((v - 3150.0).abs() < 1e-6 * 3150.0, "{v}");
}

#[test]
fn a_drum_under_a_bilinear_saddle_reads_its_closed_form() {
    let v = volume(
        r#"let m = model();
let saddle = m.sheet_nurbs([[[-10.0, -10.0, 22.5], [-10.0, 10.0, -2.5]], [[10.0, -10.0, -2.5], [10.0, 10.0, 22.5]]], #{ degree: 1 });
let d = m.cylinder(8.0, 30.0).split([saddle], containing([0.0, 0.0, 1.0]))[0];
m"#,
    );
    let want = 640.0 * std::f64::consts::PI;
    assert!((v - want).abs() < 1e-6 * want, "{v} want {want}");
}

#[test]
fn a_net_collapsed_to_a_pole_declines() {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let Err(e) = host.run(
            r#"let m = model();
m.sheet_nurbs([[[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]], #{ degree: 1 });
m"#,
            "sheet_nurbs",
    ) else {
        panic!("a collapsed edge declines");
    };
    assert!(e.to_string().contains("collapses to a point"), "{e}");
}

/// The ellipsoid x²/36 + y²/16 + z²/9 ≤ 1 as the sphere's biquadratic
/// net scaled along each axis.
const ELLIPSOID_NET: &str = r#"fn ball_net(c, r, s) {
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
}
let m = model();"#;

#[test]
fn an_ellipsoid_cut_square_to_its_axis_keeps_its_cap_at_every_height() {
    // The cut face is the ellipse x²/a² + y²/b² ≤ 1 − z²/c² at the cut,
    // its area π a b (1 − z²/c²) exact (a plane face integrates its edge
    // curves). body_volume reads the cap but takes a minute a cut and
    // reads it ~1e-3 low off the knot row (PK-1210); the STEP export
    // measures exact.
    let (a, b, c) = (6.0_f64, 4.0_f64, 3.0_f64);
    for cut in [0.0_f64, -0.1, 0.5, -0.7] {
        let source = format!(
            "{ELLIPSOID_NET}\nlet m = model();\nlet e = ball(m, [0.0, 0.0, 0.0], 1.0, [6.0, 4.0, 3.0]);\n\
             e = e.subtract(m.tool(|m| m.box(20.0, 20.0, 10.0, #{{ at: [-10.0, -10.0, {cut} - 10.0] }})));\nm"
        );
        let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
        let run = host
            .run(&source, "ellipsoid")
            .unwrap_or_else(|e| panic!("cut {cut}: {e}"));
        let body = run.model.expect("model").to_parts_all().expect("parts")[0]
            .brep()
            .clone();
        let errors = body.validate_all();
        assert!(errors.is_empty(), "cut {cut}: {errors:?}");
        assert_eq!(body.faces.len(), 2, "cut {cut}: the cap and the cut face");
        let flat = body
            .faces
            .iter()
            .find(|(_, f)| {
                matches!(
                    f.surface.as_ref().map(|s| &s.geometry),
                    Some(popkern_core::geometry::surface::SurfaceGeometry::Plane { .. })
                )
            })
            .map(|(h, _)| h)
            .expect("a cut face");
        let area = popkern_measure::face_area(&body, flat, &Tolerance::default()).unwrap();
        let want = std::f64::consts::PI * a * b * (1.0 - (cut / c).powi(2));
        assert!(
            (area - want).abs() < 1e-6 * want,
            "cut {cut}: {area} want {want}"
        );
    }
}

/// PK-1211: a block under a clamped bicubic B-spline whose boundary lies
/// flush on the block's sides. A fragment of the dome walked an edge
/// whose graph halves both belonged to a side plane, and the assembly
/// handed it the plane's chart (4e7 off the dome); it now leaves that
/// pcurve to the projection recompute on the dome itself.
#[test]
fn a_block_under_a_flush_bspline_dome_reads_its_closed_form() {
    let v = volume(
        r#"let m = model();
let xs = [0.0, 3.0, 7.0, 11.0, 16.0, 20.0];
let ys = [0.0, 5.0, 10.0, 15.0, 20.0];
let z = [[3.0, 5.0, 6.0, 5.0, 4.0], [5.0, 9.0, 11.0, 8.0, 5.0], [6.0, 12.0, 15.0, 10.0, 6.0],
    [4.0, 10.0, 18.0, 9.0, 5.0], [5.0, 8.0, 9.0, 7.0, 4.0], [2.0, 4.0, 5.0, 4.0, 3.0]];
let net = [];
for i in 0..6 { let row = []; for j in 0..5 { row.push([xs[i], ys[j], z[i][j]]); } net.push(row); }
let dome = m.sheet_nurbs(net, #{ degree: 3, knots: [[0.0, 0.0, 0.0, 0.0, 0.4, 0.4, 1.0, 1.0, 1.0, 1.0],
    [0.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 1.0, 1.0, 1.0, 1.0]] });
let b = m.box(20.0, 20.0, 30.0).split([dome], containing([10.0, 10.0, 1.0]))[0];
m"#,
    );
    let want = 554581.0 / 180.0;
    assert!((v - want).abs() < 1e-6 * want, "{v} want {want}");
}
