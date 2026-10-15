//! Positions are lengths in the model's unit, wherever a script gives
//! them: a `plane_at` origin, an `axis_at` origin (revolves, hinge
//! seams) and a selector's heights scale like an `xy_at` offset or a
//! `translate`. They were stored raw, so a mm model
//! sketched `plane_at([0, 0, 5], …)` 5 m up while `xy_at(5.0)` sat at
//! 5 mm, and revolved about `axis_at([10, 0, 0], …)` 10 m out.

use popkern_core::topology::BRep;
use popkern_rhai::ScriptHost;

fn body(source: &str) -> BRep {
    let mut host = ScriptHost::new();
    let mut outcome = host
        .run(source, "units.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    let model = outcome.model.as_mut().expect("model");
    model.to_part().expect("part").brep().clone()
}

/// (min, max) over the body's vertex coordinates.
fn bounds(brep: &BRep) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for (_, v) in brep.vertices.iter() {
        let p = [v.position.x, v.position.y, v.position.z];
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    (lo, hi)
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-12
}

#[test]
fn a_plane_at_origin_is_in_the_model_unit() {
    for (unit, metres) in [("m", 1.0), ("mm", 1e-3)] {
        let custom = body(&format!(
            "let m = model(#{{ unit: \"{unit}\" }});\n\
             m.profile(plane_at([0.0, 0.0, 5.0], [0.0, 0.0, 1.0])).rect(10.0, 10.0).sketch(\"a\").extrude(blind(1.0));\nm"
        ));
        let offset = body(&format!(
            "let m = model(#{{ unit: \"{unit}\" }});\n\
             m.profile(xy_at(5.0)).rect(10.0, 10.0).sketch(\"a\").extrude(blind(1.0));\nm"
        ));
        let (lo, hi) = bounds(&custom);
        assert!(
            close(lo[2], 5.0 * metres) && close(hi[2], 6.0 * metres),
            "{unit}: z {lo:?}..{hi:?}"
        );
        let (olo, ohi) = bounds(&offset);
        for k in 0..3 {
            assert!(
                close(olo[k], lo[k]) && close(ohi[k], hi[k]),
                "{unit}: plane_at {lo:?}..{hi:?} vs xy_at {olo:?}..{ohi:?}"
            );
        }
    }
}

#[test]
fn an_axis_at_origin_is_in_the_model_unit() {
    // A 1 mm square section at x 11..12 mm, revolved about x = 10 mm.
    let source = "let m = model(#{ unit: \"mm\" });\n\
        m.profile(xz()).path(11.0, 0.0).line_to(12.0, 0.0).line_to(12.0, 1.0).line_to(11.0, 1.0).close()\
        .sketch(\"a\").revolve(axis_at([10.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());\nm";
    let brep = body(source);
    let volume = popkern_measure::body_volume(&brep, &brep.tolerance).expect("volume");
    // The section is 1..2 mm from the axis: pi (2^2 - 1^2) * 1 mm^3,
    // in m^3. About an axis 10 m out it sweeps a ring ~10^4 larger.
    let expected = std::f64::consts::PI * (4.0 - 1.0) * 1e-9;
    assert!(
        ((volume - expected) / expected).abs() < 1e-6,
        "volume {volume} m^3, expected {expected}"
    );
}

/// Each body's volume of the hinged box example authored in `unit`, in
/// that unit cubed.
fn hinged_box_volumes(unit: &str) -> Vec<f64> {
    let source = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/hinged_box.rhai"),
    )
    .expect("example");
    let source = source.replacen(
        "let m = model();",
        &format!("let m = model(#{{ unit: \"{unit}\" }});"),
        1,
    );
    let mut host = ScriptHost::new();
    let mut outcome = host
        .run(&source, "hinged_box.rhai")
        .unwrap_or_else(|e| panic!("{unit}: {}", e.render()));
    let model = outcome.model.as_mut().expect("model");
    let per_unit3 = if unit == "mm" { 1e-9 } else { 1.0 };
    model
        .to_parts_all()
        .expect("parts")
        .iter()
        .map(|p| {
            popkern_measure::body_volume(p.brep(), &p.brep().tolerance).expect("volume") / per_unit3
        })
        .collect()
}

#[test]
fn a_millimetre_hinged_box_is_the_metre_one_scaled() {
    // The example's face selectors (`at_height(30)`) and its seam hinge's
    // `axis_at` origin are positions in the model's unit: in a mm model
    // they matched nothing / placed the seam 1000x out.
    let metre = hinged_box_volumes("m");
    let milli = hinged_box_volumes("mm");
    assert_eq!(metre.len(), milli.len());
    for (a, b) in metre.iter().zip(&milli) {
        assert!(((a - b) / a).abs() < 1e-6, "volume {b} (mm) vs {a} (m)");
    }
}
