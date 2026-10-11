//! PK-757: a plate sketched as ONE outline (a major arc round a centre
//! plus a rectangle out to one side), with a bore through the arc's
//! centre. The plate faces' loops carry only the rectangle's corners
//! while the arc bulges far past them; with vertex-only plane boxes the
//! bore's pierces were rejected AND the classifier never searched its
//! crossings there, and the body meshed open (384 open mesh edges; the
//! arm's elbow casting tore the same way with 412). Either exact
//! in-plane box alone mends it (the edge-face trim box, PK-751, and the
//! classifier window, PK-752); with both reverted it tears again. Each
//! variant serves watertight with OpenCASCADE's own cut volume.

mod union_oracle;

fn plate(cx: f64, cz: f64, ro: f64, hw: f64, reach: f64) -> String {
    let xi = cx + (ro * ro - hw * hw).sqrt();
    format!(
        "m.profile(xz_at(16.0)).path({:?}, {:?}).line_to({:?}, {:?}).line_to({xi:?}, {:?}).arc_to([{:?}, {cz:?}], [{xi:?}, {:?}]).close().sketch(\"plate\").extrude(blind(10.0))",
        cx + reach,
        cz - hw,
        cx + reach,
        cz + hw,
        cz + hw,
        cx - ro,
        cz - hw
    )
}

fn bore(cx: f64, cz: f64, r: f64) -> String {
    format!(
        "m.profile(plane_at([{cx:?}, -15.0, {cz:?}], [0.0, -1.0, 0.0])).circle({r:?}).sketch(\"b\").extrude(blind(12.0))"
    )
}

#[test]
fn an_arc_plate_bored_through_its_centre_meshes_watertight() {
    let mut failures = Vec::new();
    for (k, (cx, cz, ro, hw, reach, rb)) in [
        (-817.0, 162.5, 49.0, 25.0, 77.0, 34.05),
        (0.0, 0.0, 30.0, 12.0, 60.0, 20.0),
        (10.0, -5.0, 40.0, 30.0, 45.0, 36.0),
    ]
    .into_iter()
    .enumerate()
    {
        let name = format!("arc_plate_{k}");
        let (p, b) = (plate(cx, cz, ro, hw, reach), bore(cx, cz, rb));
        let body = format!(
            "let body = {p};\nbody.subtract(m.tool(|t| {}), #{{ label: \"bore\" }});\nm.assert_watertight();",
            b.replacen("m.", "t.", 1)
        );
        let mut model = match union_oracle::run("", &body, &format!("{name}.rhai")) {
            Ok(m) => m,
            Err(e) => {
                failures.push(e);
                continue;
            }
        };
        if let Some((served, cut)) =
            union_oracle::occ_volumes_op("", &[p.as_str(), b.as_str()], &mut model, &name, "cut")
        {
            let rel = (served - cut).abs() / cut;
            if rel > 1e-9 {
                failures.push(format!(
                    "{name}: volume {served} vs OCC cut {cut} (rel {rel:.2e})"
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}
