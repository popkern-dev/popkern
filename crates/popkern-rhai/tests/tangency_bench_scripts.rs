//! The tangent and singular boolean cases as scripts (PK-1124): every
//! case of the tangency corpus written the way a user writes it, run
//! through the Rhai host and judged against its closed-form volume and
//! area (the cases' own derivations; see each script's comment).
//!
//! A ratchet keeps the table honest, as in popkern-examples'
//! `tangency_corpus`: a `Gate` row must pass, and a `Target` row (the
//! exit of a later stage of the tangency foundations work) fails the
//! test when it starts passing, so it is promoted to a gate.

use popkern_rhai::ScriptHost;

/// Relative volume and area error a served body may carry.
const MEASURE_REL: f64 = 1e-6;

#[derive(Debug, Clone, Copy)]
enum Expect {
    /// One valid body with this volume and area.
    Body(f64, f64),
    /// A decline naming a non-manifold result (the correct result is
    /// not a manifold solid; the default policy declines).
    NonManifold,
    /// A valid body, its measures not known in closed form.
    Valid,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Status {
    Gate,
    /// A row the stage named has yet to pass (the ratchet's slot: none
    /// are open today, every bench row a gate).
    #[allow(dead_code)]
    Target(&'static str),
}

use Expect::*;
use Status::*;

/// (script, expectation, status).
const CASES: &[(&str, Expect, Status)] = &[
    (
        "ball_cut_aslant_through_a_pole",
        Body(493.199_540_828_230_97, 307.421_614_303_321_2),
        Gate,
    ),
    (
        "ball_halved_through_its_poles",
        Body(261.799_387_799_149_43, 235.619_449_019_234_48),
        Gate,
    ),
    (
        "ball_in_a_paraboloid_cup",
        Body(654.498_469_497_873_6, 455.329_476_094_537_3),
        Gate,
    ),
    (
        "ball_in_band",
        Body(580.147_443_362_915_2, 370.707_933_123_595_6),
        Gate,
    ),
    (
        "balls_meeting_through_a_pole",
        Body(938.289_005_872_151_5, 502.654_824_574_366_9),
        Gate,
    ),
    (
        "coaxial_overlap",
        Body(1_256.637_061_435_917_3, 659.734_457_253_856_6),
        Gate,
    ),
    (
        "cone_cut_aslant_through_its_apex",
        Body(210.617_317_023_784_58, 224.899_368_302_873_6),
        Gate,
    ),
    (
        "cone_halved_through_its_apex",
        Body(130.899_693_899_574_71, 177.080_092_307_881_5),
        Gate,
    ),
    (
        "cone_tip_on_face",
        Body(4_276.343_798_232_436, 1_781.613_623_880_359_5),
        Gate,
    ),
    (
        "cube_turned_a_hair",
        Body(1_000.008_726_493_952_7, 600.005_235_896_371_7),
        Gate,
    ),
    (
        "cylinder_inscribed_through_a_cube",
        Body(1_785.398_163_397_448_3, 914.159_265_358_979_3),
        Gate,
    ),
    ("grazing_cylinders", Valid, Gate),
    // The saw is a slab 1e-6 thick: no inside at the tolerance's
    // resolution, so the cut removes nothing and the cube serves whole
    // (PK-1131; the slotted body was the target before).
    ("hairline_slot", Body(1_000.0, 600.0), Gate),
    (
        "half_balls_overlapping",
        Body(392.699_081_698_724_15, 314.159_265_358_979_3),
        Gate,
    ),
    (
        "half_tori_overlapping",
        Body(1_332.396_594_147_063_4, 1_105.034_289_195_738),
        Gate,
    ),
    (
        "ice_cream_cone",
        Body(15_079.644_737_231_007, 3_015.928_947_446_201),
        Gate,
    ),
    (
        "inscribed_ball_common",
        Body(523.598_775_598_298_9, 314.159_265_358_979_3),
        Gate,
    ),
    ("inscribed_ball_union", Body(1_000.0, 600.0), Gate),
    (
        "knife_edge_dish",
        Body(523.598_775_598_298_9, 549.778_714_378_213_7),
        Gate,
    ),
    (
        "midsphere_common",
        Body(965.068_858_214_997_2, 547.466_718_050_061_6),
        Gate,
    ),
    (
        "napkin_ring_poles_on_its_rims",
        Body(268.082_573_106_329, 402.123_859_659_493_5),
        Gate,
    ),
    (
        "paraboloid_tip_a_ten_thousandth_high",
        Body(1.256_637_061_435_917_3e-7, 0.005_026_579_661_408_410_6),
        Gate,
    ),
    (
        "paraboloid_tip_on_face",
        Body(4_706.858_347_057_704, 1_824.495_160_552_916_8),
        Gate,
    ),
    (
        "pole_cap_a_ten_thousandth_high",
        Body(1.570_785_854_819_384_7e-7, 0.006_283_153_891_253_051),
        Gate,
    ),
    (
        "sectors_sharing_an_apex",
        Body(146.516_406_055_211_18, 179.057_093_587_771_33),
        Gate,
    ),
    (
        "tangent_tee_common",
        Body(239.746_428_826_596_14, 214.855_941_054_945_7),
        Gate,
    ),
    (
        "tangent_tee_union",
        Body(2_964.678_077_834_993, 1_506.736_833_112_261),
        Gate,
    ),
    (
        "thin_lens",
        Body(7.853_981_372_175_095e-12, 3.141_592_653_589_793_5e-5),
        Gate,
    ),
    (
        "torus_in_drum",
        Body(2_098.609_984_319_263, 1_442.481_399_767_068_8),
        Gate,
    ),
    (
        "torus_through_a_square_slab",
        Body(2_388.751_667_405_912_7, 1_777.260_264_867_049_4),
        Gate,
    ),
    ("villarceau_half_torus", NonManifold, Gate),
    (
        "viviani_common",
        Body(1_205.506_213_504_306_5, 628.318_530_717_958_7),
        Gate,
    ),
    ("viviani_cut", NonManifold, Gate),
    (
        "viviani_union",
        Body(4_868.239_583_435_96, 1_539.380_400_258_998_6),
        Gate,
    ),
    (
        "wedge_a_millionth_steep",
        Body(0.000_125, 100.000_075_000_025),
        Gate,
    ),
];

fn script(name: &str) -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/tangency_bench")
        .join(format!("{name}.rhai"));
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// Run one case: `Ok(worst relative measure error)` or why it failed.
fn run(name: &str, expect: Expect) -> Result<f64, String> {
    let mut host = ScriptHost::new();
    let outcome = host.run(&script(name), &format!("{name}.rhai"));
    let outcome = match (outcome, expect) {
        (Err(e), NonManifold) => {
            let text = format!("{e:?}");
            return if text.contains("non-manifold") {
                Ok(0.0)
            } else {
                Err(format!("declined, not as non-manifold: {text}"))
            };
        }
        (Err(e), _) => return Err(format!("declined: {e:?}")),
        (Ok(o), _) => o,
    };
    let mut model = outcome.model.ok_or("no model")?;
    let part = model.to_part().map_err(|e| format!("part: {e:?}"))?;
    let brep = part.brep();
    if brep.faces.is_empty() {
        return Err("served an empty body".into());
    }
    if let Some(first) = brep.validate_all().first() {
        return Err(format!("validate_all: {first}"));
    }
    let (want_v, want_a) = match expect {
        NonManifold => return Err("served, expected a non-manifold decline".into()),
        Valid => return Ok(0.0),
        Body(v, a) => (v, a),
    };
    let v = popkern_measure::body_volume(brep, &brep.tolerance).map_err(|e| format!("{e}"))?;
    let a =
        popkern_measure::body_surface_area(brep, &brep.tolerance).map_err(|e| format!("{e}"))?;
    let (ev, ea) = ((v - want_v).abs() / want_v, (a - want_a).abs() / want_a);
    // The boundary is known to the tolerance, so the volume is known to
    // about tolerance × area beyond the measure's own precision.
    let bound_v = MEASURE_REL * want_v + brep.tolerance.linear() * want_a;
    // And the area to about tolerance × perimeter (a square's for that
    // area).
    let bound_a = MEASURE_REL * want_a + 4.0 * brep.tolerance.linear() * want_a.sqrt();
    if (v - want_v).abs() > bound_v || (a - want_a).abs() > bound_a {
        return Err(format!(
            "volume {v} (want {want_v}, rel {ev:.1e}), area {a} (want {want_a}, rel {ea:.1e})"
        ));
    }
    Ok(ev.max(ea))
}

/// Seconds a script may take before it counts as failed.
const RUN_LIMIT_S: u64 = 120;

/// [`run`] in a thread with a deadline, so one stuck case cannot hide
/// the rest of the table (the thread is left to the process's end).
fn run_with_watchdog(name: &'static str, expect: Expect) -> Result<f64, String> {
    let (send, receive) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = send.send(run(name, expect));
    });
    receive
        .recv_timeout(std::time::Duration::from_secs(RUN_LIMIT_S))
        .unwrap_or_else(|_| Err(format!("did not finish in {RUN_LIMIT_S} s")))
}

#[test]
fn tangency_bench_scripts() {
    let mut problems = Vec::new();
    // `TANGENCY_BENCH_CASE=<name>` runs one case, for reproducing it.
    let only = std::env::var("TANGENCY_BENCH_CASE").ok();
    for &(name, expect, status) in CASES {
        if only.as_deref().is_some_and(|o| o != name) {
            continue;
        }
        let result = run_with_watchdog(name, expect);
        eprintln!(
            "{name:<40} {status:>16?} {}",
            match &result {
                Ok(e) => format!("pass (worst {e:.1e})"),
                Err(why) => format!("FAIL {why}"),
            }
        );
        match (status, &result) {
            (Gate, Err(why)) => problems.push(format!("gate {name} failed: {why}")),
            (Target(stage), Ok(_)) => problems.push(format!(
                "target {name} ({stage}) now passes: promote it to a gate"
            )),
            _ => {}
        }
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}
