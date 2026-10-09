//! Drive-train truth for the cycloidal joint modules (six-axis-arm
//! §9.6): sweep each module's output through a full turn and check the
//! parts that move against each other — the two discs, the ring pins,
//! the output pins and the cam shaft. No pair may interfere, and the
//! closest disc-to-ring-pin gap must stay within [b/2, 2b]: the flank
//! is the pin circle offset by the backlash b, fitted within b/4.
//!
//! Negative controls prove the check can fail: disc B phased like disc
//! A, or the rotor geared at −N instead of −(N−1), interferes.

use popkern_assembly::{Assembly, Document, InterferenceOptions, PairStatus};
use popkern_rhai::ScriptHost;

const SCRIPT: &str = include_str!("../examples/cycloidal_joint.rhai");

fn build(source: &str) -> Assembly {
    let mut host = ScriptHost::new();
    let mut outcome = host
        .run(source, "cycloidal_joint.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    match outcome.take_document() {
        Some(Document::Assembly(a)) => *a,
        _ => panic!("the script returns an assembly"),
    }
}

/// Only the S module (the controls rebuild the script).
fn s_only(source: &str) -> String {
    let edited = source.replace("for s in sizes {", "for s in [sizes[2]] {");
    assert_ne!(edited, source, "the module loop moved");
    edited
}

struct Module {
    size: &'static str,
    pins: usize,
    backlash: f64,
}

const MODULES: [Module; 3] = [
    Module {
        size: "L",
        pins: 30,
        backlash: 0.02,
    },
    Module {
        size: "M",
        pins: 26,
        backlash: 0.02,
    },
    Module {
        size: "S",
        pins: 22,
        backlash: 0.02,
    },
];

/// The drive-train occurrences of one module.
fn drive_train(m: &Module) -> Vec<String> {
    let s = m.size;
    let mut paths = vec![
        format!("disc_{s}_a"),
        format!("disc_{s}_b"),
        format!("rotor_{s}/shaft"),
    ];
    paths.extend((0..m.pins).map(|k| format!("static_{s}/ring_pin_{k}")));
    paths.extend((0..6).map(|k| format!("output_{s}/pin_{k}")));
    paths
}

/// At output angle `deg`: the interfering pairs and the smallest
/// disc-to-ring-pin distance.
/// A flagged pair: both occurrence paths and the status.
type Flagged = Vec<(String, String, String)>;

fn sample(asm: &mut Assembly, m: &Module, deg: f64) -> (Flagged, f64) {
    asm.set_assembly_param(&format!("j_{}", m.size), deg)
        .expect("slider");
    let report = asm
        .interference(&InterferenceOptions {
            clearance: 4.0 * m.backlash,
            among: drive_train(m),
            ..InterferenceOptions::default()
        })
        .expect("interference");
    let mut bad = Vec::new();
    let mut min_gap = f64::INFINITY;
    for pair in &report.pairs {
        let disc_pin = (pair.a.starts_with("disc_") && pair.b.contains("ring_pin"))
            || (pair.b.starts_with("disc_") && pair.a.contains("ring_pin"));
        match &pair.status {
            PairStatus::Interfering { .. } | PairStatus::Touching { .. } => {
                bad.push((pair.a.clone(), pair.b.clone(), format!("{:?}", pair.status)));
            }
            PairStatus::Clear { distance, .. } if disc_pin => min_gap = min_gap.min(*distance),
            PairStatus::Clear { .. } => {}
        }
    }
    (bad, min_gap)
}

/// Sweep every module through a full output turn in `step`° steps,
/// one thread per module (each on its own copy of the assembly), and
/// check the drive train at every step.
fn sweep_all(step: f64) {
    let asm = build(SCRIPT);
    let steps = (360.0 / step).round() as usize;
    std::thread::scope(|scope| {
        for m in &MODULES {
            let mut asm = asm.clone();
            scope.spawn(move || {
                let mut worst = f64::INFINITY;
                for k in 0..steps {
                    let deg = step * k as f64;
                    let (bad, gap) = sample(&mut asm, m, deg);
                    assert!(bad.is_empty(), "module {} at {deg}°: {bad:?}", m.size);
                    worst = worst.min(gap);
                }
                assert!(
                    worst >= m.backlash / 2.0 && worst <= 2.0 * m.backlash,
                    "module {}: closest disc-pin gap {worst} outside [b/2, 2b]",
                    m.size
                );
            });
        }
    });
}

/// The design's resolution (§9.6): 2° steps, 540 checks.
#[test]
fn every_module_turns_a_full_revolution_clear_of_its_pins() {
    sweep_all(2.0);
}

/// The first interfering sample of the S module over a coarse sweep.
fn first_contact(source: &str) -> Option<(f64, Flagged)> {
    let mut asm = build(source);
    let m = &MODULES[2];
    (0..36).map(|k| 10.0 * k as f64).find_map(|deg| {
        let (bad, _) = sample(&mut asm, m, deg);
        (!bad.is_empty()).then_some((deg, bad))
    })
}

#[test]
fn disc_b_phased_like_disc_a_interferes() {
    let source = s_only(SCRIPT).replace(
        "follow: j, ratio: 1.0, phase: 180.0 }",
        "follow: j, ratio: 1.0, phase: 0.0 }",
    );
    assert!(source.contains("phase: 0.0 }"), "the disc B phase moved");
    // The disc is BUILT at its phase; re-place it unturned too.
    let source = source.replace(
        ".rotate(z_axis, 180.0).translate(at - s.e",
        ".translate(at - s.e",
    );
    assert!(
        first_contact(&source).is_some(),
        "the misphased disc must hit the pins"
    );
}

#[test]
fn a_rotor_geared_at_minus_n_interferes() {
    let source = s_only(SCRIPT).replace("ratio: -(s.n - 1.0) }", "ratio: -(s.n * 1.0) }");
    assert!(
        source.contains("ratio: -(s.n * 1.0) }"),
        "the rotor ratio moved"
    );
    assert!(
        first_contact(&source).is_some(),
        "the wrong ratio must hit the pins"
    );
}

/// The script's own part functions (layout, flank fit, disc and rod
/// bodies), sliced at its section headers.
fn part_functions() -> String {
    let slice = |from: &str, to: &str| -> &str {
        let a = SCRIPT.find(from).expect("section start");
        let b = SCRIPT[a..].find(to).expect("section end") + a;
        &SCRIPT[a..b]
    };
    [
        slice("// ---- Derived layout", "// ---- Part bodies"),
        slice("fn rod(t, z", "// A shaft with the magnet hub"),
        slice(
            "// The cam bore and six output-pin windows",
            "// Cross-roller races",
        ),
    ]
    .concat()
}

#[test]
fn an_output_pin_through_its_window_never_intersects_the_disc() {
    // PK-748: with no face-face records the containment probes counted
    // ray hits inside the disc's windows, and the intersection served
    // the whole pin at output angles 6°–10°.
    for theta in [0.0, 6.0, 8.0, 10.0] {
        let source = format!(
            r#"{}
let pi = 3.141592653589793;
let s = #{{ name: "L", od: 110.0, len: 118.0, pcd: 84.0, bolts: 12, tap: 4.2, screw: "M5",
    n: 30, big_r: 40.0, rp: 2.5, e: 0.9, b: 0.02, disc_t: 6.0, cam_r: 11.0, rop: 26.0,
    op_r: 3.0, shaft_r: 8.0, d: 5.0 }};
let g = layout(s);
let fit = cycloid_fit(s, pi);
let th = {theta};
let cam = -(s.n - 1.0) * th * pi / 180.0;
let z = axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
let m = model();
let disc = disc_body(m, s, fit, pi).rotate(z, th).translate(s.e * cos(cam), s.e * sin(cam), g.z_d0);
let pin = m.tool(|t| rod(t, g.z_d0 - 1.0, s.rop, 0.0, s.op_r, g.len + 2.0 - g.z_d0, "pin")
    .rotate(z, 60.0 + th));
disc.intersect(pin, #{{ label: "overlap" }});
m"#,
            part_functions()
        );
        let mut host = ScriptHost::new();
        let mut outcome = host
            .run(&source, "window.rhai")
            .unwrap_or_else(|e| panic!("{}", e.render()));
        let model = outcome.model.as_mut().expect("model");
        assert!(
            model.to_part().is_err(),
            "at {theta}° the pin through its window intersected the disc"
        );
    }
}
