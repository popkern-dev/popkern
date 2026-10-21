//! The arm's kinematic skeleton for target posing (six-axis-arm §11
//! E4): the links on the DH axes, the target mates, seeded reachable
//! targets and the joint-space oracle. Shared by `arm_target_posing`
//! and `arm_target_limits`.
#![allow(dead_code)]

use crate::arm_support::*;
use popkern_assembly::{Assembly, Mate, MateEntity, MateKind, MateRef, Parent, Placement};
use popkern_script::Plane;
use popkern_script::prelude::*;

pub fn block() -> Model {
    let mut m = Model::new();
    m.profile(Plane::XY)
        .rect(20.0, 20.0)
        .sketch("b")
        .expect("sketch")
        .extrude(&mut m, End::Blind(20.0.into()))
        .expect("extrude");
    m
}

pub fn col(m: &M4, j: usize) -> [f64; 3] {
    [m[0][j], m[1][j], m[2][j]]
}

pub fn entity(i: &str, entity: MateEntity) -> MateRef {
    MateRef {
        instance: i.into(),
        entity,
    }
}

pub fn axis(i: &str, origin: [f64; 3], direction: [f64; 3]) -> MateRef {
    entity(i, MateEntity::Axis { origin, direction })
}

pub fn plane(i: &str, origin: [f64; 3], normal: [f64; 3]) -> MateRef {
    entity(i, MateEntity::Plane { origin, normal })
}

pub fn mate(kind: MateKind, a: MateRef, b: MateRef) -> Mate {
    Mate {
        kind,
        a,
        b,
        drive: None,
    }
}

/// The skeleton with the target `delta` (world mm) from the home
/// flange once `s` = 1.
pub fn skeleton(delta: [f64; 3]) -> Assembly {
    skeleton_with(delta, false)
}

/// [`skeleton`]; with `joints`, each revolute is a joint `j1..j6` with
/// the arm's limits, released (PK-808), instead of bare mates.
pub fn skeleton_with(delta: [f64; 3], joints: bool) -> Assembly {
    skeleton_with_limits(delta, joints, &LIMITS)
}

/// [`skeleton_with`], its released joints held to `limits`.
pub fn skeleton_with_limits(delta: [f64; 3], joints: bool, limits: &[[f64; 2]; 6]) -> Assembly {
    let f = frames(&table(), &[0.0; 6]);
    let mut asm = Assembly::new("skeleton");
    let d = asm.add_part("block", block()).expect("part");
    for name in (0..=6)
        .map(|k| format!("link_{k}"))
        .chain(["target".into()])
    {
        asm.place(Parent::Root, d, Some(&name), Placement::identity())
            .expect("place");
    }
    asm.set_ground("link_0").expect("ground");
    for k in 1..=6 {
        let (o, z) = (col(&f[k - 1], 3), col(&f[k - 1], 2));
        let (a, b) = (format!("link_{}", k - 1), format!("link_{k}"));
        if joints {
            let name = format!("j{k}");
            asm.add_revolute(
                &name,
                axis(&a, o, z),
                axis(&b, o, z),
                &popkern_assembly::RevoluteOptions {
                    angle: Some(0.0),
                    limits: Some(limits[k - 1]),
                    ..Default::default()
                },
            )
            .expect("revolute");
            asm.release_joint(&name).expect("release");
            continue;
        }
        asm.add_mate(mate(MateKind::Concentric, axis(&a, o, z), axis(&b, o, z)))
            .expect("concentric");
        asm.add_mate(mate(
            MateKind::AxialOffset { distance: 0.0 },
            axis(&a, o, z),
            axis(&b, o, z),
        ))
        .expect("seat");
    }
    asm.add_param("s", 0.0, Unit::Ratio).expect("s");
    for (i, n) in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        .into_iter()
        .enumerate()
    {
        let index = asm.mates().len();
        asm.add_mate(mate(
            MateKind::Distance {
                distance: 0.0,
                aligned: true,
            },
            plane("link_0", [0.0; 3], n),
            plane("target", [0.0; 3], n),
        ))
        .expect("target plane");
        asm.drive_mate("s", index, delta[i], 0.0).expect("drive");
    }
    let (o6, z6) = (col(&f[6], 3), col(&f[6], 2));
    asm.add_mate(mate(
        MateKind::Concentric,
        axis("link_6", o6, z6),
        axis("target", o6, z6),
    ))
    .expect("flange axis");
    asm.add_mate(mate(
        MateKind::Coincident { aligned: true },
        plane("link_6", o6, z6),
        plane("target", o6, z6),
    ))
    .expect("flange face");
    asm
}

pub fn placement(asm: &Assembly, name: &str) -> M4 {
    asm.tree()
        .iter()
        .find(|r| r.path.to_string() == name)
        .map(|r| r.world.transform().matrix)
        .unwrap_or_else(|| panic!("no {name}"))
}

pub fn apply(m: &M4, p: [f64; 3]) -> [f64; 3] {
    [0, 1, 2].map(|i| m[i][0] * p[0] + m[i][1] * p[1] + m[i][2] * p[2] + m[i][3])
}

pub fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Joint angles (degrees) from the solved placements: link k turns
/// about its home axis relative to link k − 1.
pub fn joints(asm: &Assembly) -> [f64; 6] {
    let f = frames(&table(), &[0.0; 6]);
    let mut q = [0.0; 6];
    for k in 1..=6 {
        let rel = mul(
            &inverse(&placement(asm, &format!("link_{}", k - 1))),
            &placement(asm, &format!("link_{k}")),
        );
        let (u, z) = (col(&f[k - 1], 0), col(&f[k - 1], 2));
        let ru = [0, 1, 2].map(|i| rel[i][0] * u[0] + rel[i][1] * u[1] + rel[i][2] * u[2]);
        let cross = [
            u[1] * ru[2] - u[2] * ru[1],
            u[2] * ru[0] - u[0] * ru[2],
            u[0] * ru[1] - u[1] * ru[0],
        ];
        let sin = cross[0] * z[0] + cross[1] * z[1] + cross[2] * z[2];
        let cos = u[0] * ru[0] + u[1] * ru[1] + u[2] * ru[2];
        q[k - 1] = sin.atan2(cos).to_degrees();
    }
    q
}

/// Seeded LCG in [0, 1).
pub struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// The flange's home axis and its error at `q`, scaled to millimetres
/// of a 400 mm lever (so position and direction weigh alike).
pub fn axis_error(q: &[f64; 6]) -> [f64; 3] {
    let t = table();
    let home = col(&frames(&t, &[0.0; 6])[6], 2);
    let z = col(&frames(&t, q)[6], 2);
    [0, 1, 2].map(|i| 400.0 * (z[i] - home[i]))
}

/// A reachable target: random q1..q5, then (q4, q5) by Gauss–Newton so
/// the flange axis is the home axis; None when that does not converge.
pub fn reachable(rng: &mut Lcg) -> Option<[f64; 3]> {
    let mut q = [0.0; 6];
    for (k, span) in [360.0, 360.0, 320.0, 360.0, 360.0].into_iter().enumerate() {
        q[k] = span * (rng.next() - 0.5);
    }
    reachable_from(q)
}

/// The flange position at `q` once (q4, q5) bring its axis to the home
/// axis by Gauss-Newton; None when that does not converge.
pub fn reachable_from(mut q: [f64; 6]) -> Option<[f64; 3]> {
    for _ in 0..60 {
        let e = axis_error(&q);
        if e.iter().map(|v| v * v).sum::<f64>().sqrt() < 1e-9 {
            return Some(col(&frames(&table(), &q)[6], 3));
        }
        let mut j = [[0.0; 2]; 3];
        for (c, idx) in [3, 4].into_iter().enumerate() {
            let mut qp = q;
            qp[idx] += 1e-6;
            let ep = axis_error(&qp);
            for r in 0..3 {
                j[r][c] = (ep[r] - e[r]) / 1e-6;
            }
        }
        let dot = |a: usize, b: usize| (0..3).map(|r| j[r][a] * j[r][b]).sum::<f64>();
        let (a11, a12, a22) = (dot(0, 0), dot(0, 1), dot(1, 1));
        let (b1, b2) = (
            (0..3).map(|r| j[r][0] * e[r]).sum::<f64>(),
            (0..3).map(|r| j[r][1] * e[r]).sum::<f64>(),
        );
        let det = a11 * a22 - a12 * a12;
        if det.abs() < 1e-14 {
            return None;
        }
        q[3] -= (a22 * b1 - a12 * b2) / det;
        q[4] -= (a11 * b2 - a12 * b1) / det;
    }
    None
}

/// Damped least squares on the DH chain from `q` (joint space, the
/// independent oracle): whether it puts the flange at `p` on its home
/// axis.
pub fn oracle_reaches(q: &[f64; 6], p: [f64; 3]) -> bool {
    let residual = |q: &[f64; 6]| {
        let o = col(&frames(&table(), q)[6], 3);
        let a = axis_error(q);
        [o[0] - p[0], o[1] - p[1], o[2] - p[2], a[0], a[1], a[2]]
    };
    let norm = |r: &[f64; 6]| r.iter().map(|v| v * v).sum::<f64>().sqrt();
    let (mut q, mut lambda) = (*q, 1.0);
    for _ in 0..400 {
        let r = residual(&q);
        if norm(&r) < 1e-9 {
            return true;
        }
        let mut j = [[0.0; 6]; 6];
        for c in 0..6 {
            let mut qp = q;
            qp[c] += 1e-6;
            let rp = residual(&qp);
            for i in 0..6 {
                j[i][c] = (rp[i] - r[i]) / 1e-6;
            }
        }
        // (JᵀJ + λI) dq = −Jᵀr by Gauss–Jordan with partial pivoting.
        let mut a = [[0.0; 7]; 6];
        for i in 0..6 {
            for k in 0..6 {
                a[i][k] = (0..6).map(|m| j[m][i] * j[m][k]).sum::<f64>();
            }
            a[i][i] += lambda;
            a[i][6] = -(0..6).map(|m| j[m][i] * r[m]).sum::<f64>();
        }
        for c in 0..6 {
            let pivot = (c..6)
                .max_by(|x, y| a[*x][c].abs().total_cmp(&a[*y][c].abs()))
                .unwrap_or(c);
            a.swap(c, pivot);
            let row = a[c];
            for (_, target) in a.iter_mut().enumerate().filter(|(i, _)| *i != c) {
                let f = target[c] / row[c];
                for (t, r) in target.iter_mut().zip(row).skip(c) {
                    *t -= f * r;
                }
            }
        }
        let mut next = q;
        for c in 0..6 {
            next[c] += a[c][6] / a[c][c];
        }
        if norm(&residual(&next)) < norm(&r) {
            q = next;
            lambda = (lambda * 0.3).max(1e-9);
        } else {
            lambda *= 10.0;
        }
    }
    false
}

/// The 100 seeded reachable targets (world flange positions).
pub fn targets() -> Vec<[f64; 3]> {
    let mut rng = Lcg(0x05ee_d726);
    let mut out = Vec::new();
    while out.len() < 100 {
        out.extend(reachable(&mut rng));
    }
    out
}

pub fn offset(p: [f64; 3], home: [f64; 3]) -> [f64; 3] {
    [p[0] - home[0], p[1] - home[1], p[2] - home[2]]
}
