//! Kinematic truth for the six-axis arm (design
//! `docs/architecture/six-axis-arm.md` §5.1, §9.5): every link placement
//! the mate solve finds matches the Denavit–Hartenberg oracle at 200
//! seeded poses inside the joint limits and at the wrist singularity,
//! and a table with d4 off by 1 mm does not. Also the slider budget
//! (§8: a joint re-solve ≤ 5 ms), the tube-length sliders, and the
//! script sharing its joint modules with `cycloidal_joint.rhai`.

use popkern_assembly::{Assembly, Document};
use popkern_rhai::ScriptHost;

const SCRIPT: &str = include_str!("../examples/six_axis_arm.rhai");
const MODULES: &str = include_str!("../examples/cycloidal_joint.rhai");

/// The arm at home (the displayed-pose line set to `home`).
fn build() -> Assembly {
    let source = SCRIPT.replace("\nready(a);\na\n", "\nhome(a);\na\n");
    assert_ne!(source, SCRIPT, "the displayed-pose line moved");
    let mut host = ScriptHost::new();
    let mut outcome = host
        .run(&source, "six_axis_arm.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    match outcome.take_document() {
        Some(Document::Assembly(a)) => *a,
        _ => panic!("the script returns an assembly"),
    }
}

type M4 = [[f64; 4]; 4];

/// The §4.1 table: d, a, α and the joint limits (degrees).
struct Dh {
    d: [f64; 6],
    a: [f64; 6],
    alpha: [f64; 6],
}

const LIMITS: [[f64; 2]; 6] = [
    [-360.0, 360.0],
    [-360.0, 360.0],
    [-170.0, 170.0],
    [-360.0, 360.0],
    [-360.0, 360.0],
    [-360.0, 360.0],
];

fn table() -> Dh {
    let h = std::f64::consts::FRAC_PI_2;
    Dh {
        d: [162.5, 0.0, 0.0, 133.3, 99.7, 99.6],
        a: [0.0, -425.0, -392.0, 0.0, 0.0, 0.0],
        alpha: [h, 0.0, 0.0, h, -h, 0.0],
    }
}

fn step(theta: f64, d: f64, a: f64, alpha: f64) -> M4 {
    let (c, s, ca, sa) = (theta.cos(), theta.sin(), alpha.cos(), alpha.sin());
    [
        [c, -s * ca, s * sa, a * c],
        [s, c * ca, -c * sa, a * s],
        [0.0, sa, ca, d],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

fn mul(a: &M4, b: &M4) -> M4 {
    let mut r = [[0.0; 4]; 4];
    for (i, row) in r.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v = (0..4).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    r
}

fn inverse(m: &M4) -> M4 {
    let mut r = [[0.0; 4]; 4];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = m[j][i];
        }
        r[i][3] = -(0..3).map(|k| m[k][i] * m[k][3]).sum::<f64>();
    }
    r[3][3] = 1.0;
    r
}

/// DH frames 0..6 at joint angles `q` (degrees).
fn frames(t: &Dh, q: &[f64; 6]) -> Vec<M4> {
    let mut f = vec![step(0.0, 0.0, 0.0, 0.0)];
    for k in 0..6 {
        let next = mul(&f[k], &step(q[k].to_radians(), t.d[k], t.a[k], t.alpha[k]));
        f.push(next);
    }
    f
}

/// Each root and the DH frame it rides: a link authored at home moves
/// by T_k(q) · T_k(0)⁻¹. The tube roots ride the same frame as their
/// far end castings (the tube length is fixed here).
const ROOTS: [(&str, usize); 9] = [
    ("link_0", 0),
    ("link_1", 1),
    ("link_2_tube", 2),
    ("link_2", 2),
    ("link_3_tube", 3),
    ("link_3", 3),
    ("link_4", 4),
    ("link_5", 5),
    ("link_6", 6),
];

/// The largest disagreement between the solved root placements and the
/// oracle: translation in mm, rotation entries scaled by the reach.
fn worst(asm: &Assembly, t: &Dh, q: &[f64; 6], reach: f64) -> f64 {
    let home = frames(&table(), &[0.0; 6]);
    let now = frames(t, q);
    let rows = asm.tree();
    let mut worst: f64 = 0.0;
    for (name, k) in ROOTS {
        let want = mul(&now[k], &inverse(&home[k]));
        let got = rows
            .iter()
            .find(|r| r.path.to_string() == name)
            .map(|r| r.world.transform().matrix)
            .unwrap_or_else(|| panic!("no root {name}"));
        for i in 0..3 {
            for j in 0..4 {
                let scale = if j == 3 { 1.0 } else { reach };
                worst = worst.max((got[i][j] - want[i][j]).abs() * scale);
            }
        }
    }
    worst
}

fn pose(asm: &mut Assembly, q: &[f64; 6]) {
    for (k, v) in q.iter().enumerate() {
        asm.set_assembly_param(&format!("j{}", k + 1), *v)
            .unwrap_or_else(|e| panic!("j{} = {v}: {e}", k + 1));
    }
}

/// 200 poses from a fixed linear congruential sequence, inside the
/// limits, then the wrist singularity (j5 = 0: j4 and j6 collinear).
fn poses() -> Vec<[f64; 6]> {
    let mut state: u64 = 0x5eed_a12e;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut out: Vec<[f64; 6]> = (0..200)
        .map(|_| {
            let mut q = [0.0; 6];
            for (k, v) in q.iter_mut().enumerate() {
                let [lo, hi] = LIMITS[k];
                *v = lo + (hi - lo) * next();
            }
            q
        })
        .collect();
    out.push([37.0, -71.0, 104.0, -122.0, 0.0, 215.0]);
    out
}

#[test]
fn the_solved_links_match_the_dh_oracle_and_a_wrong_table_does_not() {
    let mut asm = build();
    let reach = 850.0;
    let gate = 1e-9 * reach;
    let mut wrong = table();
    wrong.d[3] += 1.0;
    let mut times = Vec::new();
    for q in poses() {
        let t = std::time::Instant::now();
        pose(&mut asm, &q);
        times.push(t.elapsed().as_secs_f64() / 6.0);
        let err = worst(&asm, &table(), &q, reach);
        assert!(
            err <= gate,
            "pose {q:?}: links off the oracle by {err:e} mm"
        );
        // Negative control: d4 + 1 mm moves the wrist links by ~1 mm.
        let off = worst(&asm, &wrong, &q, reach);
        assert!(
            off > 0.5,
            "pose {q:?}: a d4 + 1 mm table still matched ({off})"
        );
    }
    // A joint move re-solves within the slider budget (§8: ≤ 5 ms).
    times.sort_by(f64::total_cmp);
    let median = times[times.len() / 2] * 1e3;
    eprintln!(
        "joint re-solve: median {median:.2} ms over {} poses",
        times.len()
    );
    assert!(median <= 5.0, "joint re-solve median {median:.2} ms");
}

#[test]
fn the_tube_sliders_move_the_far_joints_and_rebuild_only_the_tubes() {
    let mut asm = build();
    let at = |asm: &Assembly, name: &str| {
        asm.tree()
            .iter()
            .find(|r| r.path.to_string() == name)
            .map(|r| r.world.origin())
            .unwrap_or_else(|| panic!("no {name}"))
    };
    let (elbow, wrist) = (at(&asm, "link_2"), at(&asm, "link_3"));
    let t = std::time::Instant::now();
    let rebuilt = asm
        .set_assembly_param("upper_arm", 450.0)
        .expect("upper_arm");
    let elapsed = t.elapsed();
    assert_eq!(rebuilt, vec!["upper_arm_tube".to_string()]);
    // §8: one tube replay and a re-solve within 1.5 s.
    assert!(
        elapsed.as_secs_f64() <= 1.5,
        "upper_arm slider took {elapsed:?}"
    );
    // The far end casting and everything beyond move 25 mm out along
    // the arm (−x at home).
    let (e2, w2) = (at(&asm, "link_2"), at(&asm, "link_3"));
    assert!((e2.x - elbow.x + 25.0).abs() < 1e-9 && (e2.y - elbow.y).abs() < 1e-9);
    assert!((w2.x - wrist.x + 25.0).abs() < 1e-9);
    let rebuilt = asm.set_assembly_param("forearm", 380.0).expect("forearm");
    assert_eq!(rebuilt, vec!["forearm_tube".to_string()]);
    assert!((at(&asm, "link_3").x - w2.x - 12.0).abs() < 1e-9);
    // Below the end castings' reach the slider refuses.
    assert!(asm.set_assembly_param("upper_arm", 200.0).is_err());
}

#[test]
fn the_arm_shares_its_joint_modules_with_the_module_showcase() {
    let block = |s: &'static str| {
        let a = s.find("// ---- Size table").expect("size table");
        let b = s[a..]
            .find("// ---- One module: sub-assemblies")
            .map(|b| a + b);
        (a, b)
    };
    let (a, b) = block(MODULES);
    let shared = &MODULES[a..b.expect("module block end")];
    assert!(
        SCRIPT.contains(shared),
        "six_axis_arm.rhai must carry cycloidal_joint.rhai's module block verbatim"
    );
}

/// Joint k's static and output link roots (1-based joints).
const JOINT_LINKS: [(&str, &str); 6] = [
    ("link_0", "link_1"),
    ("link_1", "link_2_tube"),
    ("link_2", "link_3_tube"),
    ("link_3", "link_4"),
    ("link_4", "link_5"),
    ("link_5", "link_6"),
];

/// Ring pins per joint (module sizes L, L, M, S, S, S).
const PINS: [usize; 6] = [30, 30, 26, 22, 22, 22];

#[test]
fn every_joint_drive_train_sits_as_in_the_module_showcase() {
    use popkern_assembly::{InterferenceOptions, PairStatus};
    let mut asm = build();
    let backlash = 0.02;
    for q in [[0.0; 6], [37.0, -71.0, 104.0, -122.0, 0.0, 215.0]] {
        pose(&mut asm, &q);
        for k in 1..=6 {
            let (st, out) = JOINT_LINKS[k - 1];
            let mut among = vec![
                format!("j{k}_disc_a"),
                format!("j{k}_disc_b"),
                format!("j{k}_rotor/shaft"),
            ];
            among.extend((0..PINS[k - 1]).map(|i| format!("{st}/j{k}_ring_pin_{i}")));
            among.extend((0..6).map(|i| format!("{out}/j{k}_output_pin_{i}")));
            let report = asm
                .interference(&InterferenceOptions {
                    clearance: 4.0 * backlash,
                    among,
                    ..InterferenceOptions::default()
                })
                .expect("interference");
            let mut gap = f64::INFINITY;
            for p in &report.pairs {
                match &p.status {
                    PairStatus::Clear { distance, .. } => {
                        if p.a.contains("disc") != p.b.contains("disc")
                            && (p.a.contains("ring_pin") || p.b.contains("ring_pin"))
                        {
                            gap = gap.min(*distance);
                        }
                    }
                    other => panic!("pose {q:?} j{k}: {} x {}: {other:?}", p.a, p.b),
                }
            }
            assert!(
                (backlash / 2.0..=2.0 * backlash).contains(&gap),
                "pose {q:?} j{k}: disc-to-pin gap {gap} outside [b/2, 2b]"
            );
        }
    }
}

#[test]
fn every_joint_swings_its_full_range_and_the_base_yaw_turns_clear() {
    use popkern_assembly::SwingOptions;
    let mut asm = build();
    for k in 1..=6 {
        let (st, out) = JOINT_LINKS[k - 1];
        let [from, to] = LIMITS[k - 1];
        let t = std::time::Instant::now();
        let report = asm
            .swing_with(
                &format!("j{k}"),
                &SwingOptions {
                    from,
                    to,
                    step: 5.0,
                    resolution: 0.5,
                    clearance: Some(0.5),
                    // The joint's own races, flange and housing turn
                    // about their common axis: the module's business.
                    exclude: vec![[format!("{st}/j{k}_"), format!("{out}/j{k}_")]],
                },
            )
            .unwrap_or_else(|e| panic!("j{k} swing: {e}"));
        let elapsed = t.elapsed();
        eprintln!(
            "j{k} swing [{from}, {to}]: {elapsed:?}, {} samples, first contact {:?} {:?}, closest {:?} at {:?}",
            report.samples.len(),
            report.first_contact,
            report.first_contact_pair,
            report.min_clearance,
            report.min_clearance_at
        );
        // The swing restores the joint.
        let p = asm
            .params()
            .iter()
            .find(|p| p.name == format!("j{k}"))
            .expect("param");
        assert!(p.display_value().abs() < 1e-9);
        if k == 1 {
            // Turning the whole arm about the base axis meets nothing.
            assert_eq!(
                report.first_contact, None,
                "the base yaw collided: {:?}",
                report.first_contact_pair
            );
        }
        if k == 2 {
            // §8: a full shoulder swing within 20 s.
            assert!(elapsed.as_secs_f64() <= 20.0, "j2 swing took {elapsed:?}");
        }
    }
}

/// The first colliding pair at pose `q`, over every pair the arm's
/// joints separate: a swing of each joint through 0.01° from its pose
/// checks all pairs across that joint (upstream against downstream),
/// and every pair of parts in different links is across some joint.
fn collision_at(asm: &mut Assembly, q: &[f64; 6]) -> Option<(String, String)> {
    use popkern_assembly::SwingOptions;
    pose(asm, q);
    for k in 1..=6 {
        let (st, out) = JOINT_LINKS[k - 1];
        let from = q[k - 1];
        let to = if from + 0.01 <= LIMITS[k - 1][1] {
            from + 0.01
        } else {
            from - 0.01
        };
        let report = asm
            .swing_with(
                &format!("j{k}"),
                &SwingOptions {
                    from,
                    to,
                    step: 0.01,
                    resolution: 0.01,
                    clearance: Some(0.0),
                    exclude: vec![[format!("{st}/j{k}_"), format!("{out}/j{k}_")]],
                },
            )
            .unwrap_or_else(|e| panic!("j{k} at {q:?}: {e}"));
        if report.first_contact.is_some() {
            return report.first_contact_pair;
        }
    }
    None
}

#[test]
fn the_pose_presets_are_collision_free_and_a_folded_elbow_is_not() {
    let mut asm = build();
    for (name, q) in [
        ("home", [0.0; 6]),
        ("ready", [0.0, -60.0, 100.0, -130.0, -90.0, 0.0]),
        ("reach", [45.0, -20.0, 30.0, -100.0, -90.0, 0.0]),
    ] {
        // The script's own preset sets exactly these angles.
        if name != "home" {
            let line = format!(
                "fn {name}(a) {{\n    let q = [{}];",
                q.map(|v| format!("{v:.1}")).join(", ")
            );
            assert!(SCRIPT.contains(&line), "preset {name} is not {q:?}");
        }
        assert_eq!(collision_at(&mut asm, &q), None, "preset {name} collides");
    }
    // Negative control: the elbow folded to its limit drives the
    // forearm's end into the shoulder.
    let folded = [0.0, 0.0, -170.0, 0.0, 0.0, 0.0];
    assert!(
        collision_at(&mut asm, &folded).is_some(),
        "the folded elbow must collide"
    );
}
