//! Engraving around an island: a box minus a ring-shaped cutter (a
//! square ring poking through the top face, so the ring's hole leaves an
//! island standing in the pocket) serves at every scale with exactly the
//! ring's area times the cut depth removed. This is the shape of every
//! engraved 0, 4, 6, 8 and 9.
//!
//! The top face's ring-shaped fragment lies inside the cutter. Its
//! offset-pair probes at the fine scale said so, but the coarse retry
//! probed past the cutter's 0.2 top and 0.5 floor, read Outside on both
//! sides and kept the fragment; the cut then declined on curveless
//! edges. Below ~0.02 m it always did; at 0.03 the coarse probe still
//! fit inside the cutter and it served.

use popkern_core::topology::BRep;
use popkern_rhai::ScriptHost;

fn engraved(unit: &str, k: f64) -> BRep {
    let source = format!(
        "let m = model(#{{ unit: \"{unit}\" }});\n\
         let k = {k};\n\
         let die = m.profile(xy_at(-8.0 * k)).rect(16.0 * k, 16.0 * k).sketch(\"block\").extrude(blind(16.0 * k));\n\
         let tool = m.tool(|t| {{\n\
             let ink = t.profile(xy_at(7.5 * k)).rect(6.0 * k, 6.0 * k).sketch(\"ink\").extrude(blind(0.7 * k));\n\
             ink.subtract(t.profile(xy_at(7.3 * k)).rect(3.0 * k, 3.0 * k).sketch(\"counter\").extrude(blind(1.1 * k)))\n\
         }});\n\
         die = die.subtract(tool, #{{ label: \"engrave\" }});\nm"
    );
    let mut host = ScriptHost::new();
    let mut outcome = host
        .run(&source, "engrave.rhai")
        .unwrap_or_else(|e| panic!("k {k} {unit}: {}", e.render()));
    let model = outcome.model.as_mut().expect("model");
    model
        .to_part()
        .unwrap_or_else(|e| panic!("k {k} {unit}: {e}"))
        .brep()
        .clone()
}

fn check(unit: &str, metres_per_unit: f64, k: f64) {
    let brep = engraved(unit, k);
    assert!(brep.validate().is_ok(), "k {k} {unit}: invalid");
    let volume = popkern_measure::body_volume(&brep, &brep.tolerance).expect("volume");
    // The cut reaches from the top (8k) down to the cutter floor (7.5k).
    let s = k * metres_per_unit;
    let expected = (16.0 * s).powi(3) - (36.0 - 9.0) * s * s * 0.5 * s;
    assert!(
        ((volume - expected) / expected).abs() < 1e-9,
        "k {k} {unit}: volume {volume} vs {expected}"
    );
}

#[test]
fn a_ring_engraving_leaves_its_island_at_every_scale() {
    for k in [1.0, 0.03, 0.01, 0.001] {
        check("m", 1.0, k);
    }
    check("mm", 1e-3, 1.0);
}
