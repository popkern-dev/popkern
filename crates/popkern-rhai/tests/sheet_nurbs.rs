//! PK-1200: `m.sheet_nurbs` makes the sheet of a B-spline surface from
//! its control net, and a solid split by it keeps the exact surface.
//! A Bézier dome over a block (its boundary curves flush on the block's
//! sides) and the bilinear saddle over a drum read their closed forms;
//! a net with a collapsed edge declines (only open patches are sheets).

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
