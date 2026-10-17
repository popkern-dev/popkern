//! The scale sweep (PK-843): one script at k = 1e-3, 1 and 1e3 metres
//! must serve the same body — validate_all clean, the same face kinds,
//! and the same volume over k³. An absolute length constant shows as a
//! scale where one of these breaks (measured before the fixes here: the
//! volume quadrature's 0.001 mesh chord cost the drilled dome 28x its
//! time at k = 1e3 and 2.4e-7 of its volume at k = 1e-3).

use popkern_core::topology::BRep;

fn served(src: &str) -> BRep {
    let mut run = popkern_rhai::ScriptHost::new()
        .run(src, "scaled.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    let parts = run
        .model
        .take()
        .expect("model")
        .to_parts_all()
        .expect("parts");
    parts[0].brep().clone()
}

fn kinds(b: &BRep) -> Vec<&'static str> {
    let mut k: Vec<_> = b
        .faces
        .values()
        .filter_map(|f| f.surface.as_ref().map(|s| s.geometry.kind_name()))
        .collect();
    k.sort_unstable();
    k
}

fn sweep(name: &str, script: &str) {
    let mut base: Option<(f64, Vec<&'static str>)> = None;
    for k in [1.0, 1e-3, 1e3] {
        let b = served(&script.replace("K", &format!("{k:.1e}")));
        let errors = b.validate_all();
        assert!(errors.is_empty(), "{name} k={k:e}: {errors:?}");
        let v = popkern_measure::body_volume(&b, &b.tolerance).expect("volume") / (k * k * k);
        let (v0, k0) = base.get_or_insert_with(|| (v, kinds(&b)));
        assert_eq!(&kinds(&b), k0, "{name} k={k:e}: face kinds");
        assert!(
            ((v - *v0) / *v0).abs() < 1e-9,
            "{name} k={k:e}: volume/k^3 {v} vs {v0}"
        );
    }
}

#[test]
fn a_drilled_dome_is_the_same_body_at_every_scale() {
    sweep(
        "dome",
        r#"let m = model();
let body = m.profile(xz()).path(0.0, 0.0).line_to(K, 0.0).line_to(K, 2.0*K)
    .arc_to([0.70710678118654752*K, 2.70710678118654752*K], [0.0, 3.0*K]).close()
    .sketch("d").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
let hole = m.tool(|t| t.profile(xy_at(K)).circle(0.3*K).sketch("h").extrude(blind(3.0*K)));
body = body.subtract(hole);
m"#,
    );
}

#[test]
fn a_fully_filleted_block_is_the_same_body_at_every_scale() {
    // Every edge rounded: the corners serve sphere caps through the
    // surface-fix sphere gate, which read 0.5 m on a millimetre corner.
    sweep(
        "rounded block",
        r#"let m = model();
let b = m.profile(xy()).rect(2.0*K, 2.0*K).sketch("s").extrude(blind(K));
b = b.fillet(edges(), 0.2*K);
m"#,
    );
}
