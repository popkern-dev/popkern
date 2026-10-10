//! PK-751: a pocket cut into a disc face whose outer loop carries only a
//! few vertices (a slot opened from the rim) while its circular rim
//! bulges far beyond them. The pocket's seam line pierces the disc face
//! inside the disc but outside the vertices' parameter box; the
//! edge-face trim test rejected that pierce, the pocket's record ring
//! crossed its seam unsplit, the wall's chart overwrapped and the cut
//! declined. The trim box of a plane now bounds every edge's exact
//! extent. Each variant serves, and its volume matches OpenCASCADE's own
//! cut of the same operands.

mod union_oracle;

const PRELUDE: &str = r#"
fn rod_at(t, p, n, r, length, label) {
    t.profile(plane_at(p, n)).circle(r).sketch(label).extrude(blind(length))
}
"#;

/// The seat disc (r 67, axis y, y 0..10), a slot rising from its bottom
/// rim whose corner pokes into the pocket circle at `(cx, cz)` radius
/// `r`, then the pocket along y.
fn variant(cx: f64, cz: f64, r: f64, inset: f64) -> (String, [String; 3]) {
    let seat = "m.profile(plane_at([0.0, 10.0, 162.5], [0.0, -1.0, 0.0])).circle(67.0).sketch(\"seat\").extrude(blind(10.0))".to_string();
    let (x1, z1) = (cx - inset, cz - inset);
    let slot = format!(
        "m.profile(xy_at(80.0)).rect(17.0, 30.0).at({:?}, 0.0).sketch(\"slot\").extrude(blind({:?}))",
        x1 - 8.5,
        z1 - 80.0
    );
    let pocket =
        format!("rod_at(m, [{cx:?}, -30.0, {cz:?}], [0.0, 1.0, 0.0], {r:?}, 36.0, \"pocket\")");
    let tool = |e: &str| e.replacen("m.", "t.", 1).replacen("(m,", "(t,", 1);
    let body = format!(
        "let b = {seat};\n\
         b = b.subtract(m.tool(|t| {}), #{{ label: \"slot\" }});\n\
         b = b.subtract(m.tool(|t| {}), #{{ label: \"pocket\" }});\n\
         m.assert_watertight();",
        tool(&slot),
        tool(&pocket)
    );
    (body, [seat, slot, pocket])
}

#[test]
fn pockets_cut_past_a_rim_slot_corner_serve_with_the_cut_volume() {
    let mut failures = Vec::new();
    for (k, (cx, cz, r, inset)) in [
        (0.0, 222.0, 5.0, 3.0),
        (10.0, 215.0, 4.0, 2.0),
        (-5.0, 210.0, 6.0, 3.5),
    ]
    .into_iter()
    .enumerate()
    {
        let name = format!("rim_slot_{k}");
        let (body, operands) = variant(cx, cz, r, inset);
        let mut model = match union_oracle::run(PRELUDE, &body, &format!("{name}.rhai")) {
            Ok(m) => m,
            Err(e) => {
                failures.push(e);
                continue;
            }
        };
        let report = popkern_rhai::report_value(&model);
        if report["failures"] != 0 {
            failures.push(format!("{name}: {report}"));
            continue;
        }
        let refs: Vec<&str> = operands.iter().map(String::as_str).collect();
        if let Some((served, cut)) =
            union_oracle::occ_volumes_op(PRELUDE, &refs, &mut model, &name, "cut")
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
