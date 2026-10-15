//! Positions are lengths in the model's unit, wherever a script gives
//! them: a `plane_at` origin and an `axis_at` origin scale like an
//! `xy_at` offset or a `translate`. They were stored raw, so a mm model
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
