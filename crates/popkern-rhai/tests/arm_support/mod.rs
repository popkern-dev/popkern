//! Shared by the six-axis arm tests: the script, its build at home, the
//! Denavit–Hartenberg oracle (design §5.1), joint poses, and the
//! cross-link collision probe.
#![allow(dead_code)]

use popkern_assembly::{Assembly, Document};
use popkern_rhai::ScriptHost;

pub const SCRIPT: &str = include_str!("../../examples/six_axis_arm.rhai");
pub const MODULES: &str = include_str!("../../examples/cycloidal_joint.rhai");

/// The arm at home (the displayed-pose line set to `home`).
pub fn build() -> Assembly {
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

pub type M4 = [[f64; 4]; 4];

/// The §4.1 table: d, a, α and the joint limits (degrees).
pub struct Dh {
    pub d: [f64; 6],
    pub a: [f64; 6],
    pub alpha: [f64; 6],
}

pub const LIMITS: [[f64; 2]; 6] = [
    [-360.0, 360.0],
    [-360.0, 360.0],
    [-170.0, 170.0],
    [-360.0, 360.0],
    [-360.0, 360.0],
    [-360.0, 360.0],
];

pub fn table() -> Dh {
    let h = std::f64::consts::FRAC_PI_2;
    Dh {
        d: [162.5, 0.0, 0.0, 133.3, 99.7, 99.6],
        a: [0.0, -425.0, -392.0, 0.0, 0.0, 0.0],
        alpha: [h, 0.0, 0.0, h, -h, 0.0],
    }
}

pub fn step(theta: f64, d: f64, a: f64, alpha: f64) -> M4 {
    let (c, s, ca, sa) = (theta.cos(), theta.sin(), alpha.cos(), alpha.sin());
    [
        [c, -s * ca, s * sa, a * c],
        [s, c * ca, -c * sa, a * s],
        [0.0, sa, ca, d],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

pub fn mul(a: &M4, b: &M4) -> M4 {
    let mut r = [[0.0; 4]; 4];
    for (i, row) in r.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v = (0..4).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    r
}

pub fn inverse(m: &M4) -> M4 {
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
pub fn frames(t: &Dh, q: &[f64; 6]) -> Vec<M4> {
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
pub const ROOTS: [(&str, usize); 9] = [
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
pub fn worst(asm: &Assembly, t: &Dh, q: &[f64; 6], reach: f64) -> f64 {
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

pub fn pose(asm: &mut Assembly, q: &[f64; 6]) {
    for (k, v) in q.iter().enumerate() {
        asm.set_assembly_param(&format!("j{}", k + 1), *v)
            .unwrap_or_else(|e| panic!("j{} = {v}: {e}", k + 1));
    }
}

/// 200 poses from a fixed linear congruential sequence, inside the
/// limits, then the wrist singularity (j5 = 0: j4 and j6 collinear).
pub fn poses() -> Vec<[f64; 6]> {
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

/// Joint k's static and output link roots (1-based joints).
pub const JOINT_LINKS: [(&str, &str); 6] = [
    ("link_0", "link_1"),
    ("link_1", "link_2_tube"),
    ("link_2", "link_3_tube"),
    ("link_3", "link_4"),
    ("link_4", "link_5"),
    ("link_5", "link_6"),
];

/// Ring pins per joint (module sizes L, L, M, S, S, S).
pub const PINS: [usize; 6] = [30, 30, 26, 22, 22, 22];

/// The first colliding pair at pose `q`, over every pair the arm's
/// joints separate: a swing of each joint through 0.01° from its pose
/// checks all pairs across that joint (upstream against downstream),
/// and every pair of parts in different links is across some joint.
pub fn collision_at(asm: &mut Assembly, q: &[f64; 6]) -> Option<(String, String)> {
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
