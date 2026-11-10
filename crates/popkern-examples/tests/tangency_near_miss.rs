//! Near-miss sweep of the tangency corpus (PK-1128, the plan's §8.2
//! probe): each touch case with one operand moved along its contact's
//! normal by δ ∈ {±ε/2, ±2ε, ±10ε}, ε the linear tolerance.
//!
//! Within the tolerance band (|δ| = ε/2) the geometry is the touch to
//! every downstream consumer, and the outcome must be the δ = 0 one
//! (the case's closed form, judged as the corpus judges it). Past the
//! band the operands genuinely overlap or part, and the outcome may
//! change, but only honestly: a served body passes `validate_all` and
//! its volume stays within |δ| times the operands' area of the δ = 0
//! answer (a body cannot gain or lose more by moving a face |δ|), or
//! the boolean declines with a typed reason about the contact
//! (`NonManifoldResult`, `RelationUndecided`). An internal failure, a
//! validation failure or a hang fails the row. The sweep runs at the
//! identity placement and a rigid turn.

mod tangency_corpus_support;

use popkern_bool::{BoolError, BoolOptions, NonManifoldPolicy, boolean_with_options};
use popkern_core::geometry::{Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use tangency_corpus_support::{Expect, cases};

/// The touch cases swept: the operand moved and the contact's normal.
const SWEPT: &[(&str, usize, (f64, f64, f64))] = &[
    ("inscribed-ball-common", 1, (1.0, 0.0, 0.0)),
    ("inscribed-ball-union", 1, (1.0, 0.0, 0.0)),
    ("inscribed-ball-cut", 1, (1.0, 0.0, 0.0)),
    ("ball-resting-common", 1, (0.0, 0.0, 1.0)),
    ("ball-resting-union", 1, (0.0, 0.0, 1.0)),
    ("dimple-touching-edge", 1, (1.0, 0.0, 0.0)),
    ("ice-cream-cone", 1, (0.0, 0.0, 1.0)),
    (
        "ball-resting-on-an-edge-union",
        1,
        (0.707_106_781_186_547_5, 0.0, 0.707_106_781_186_547_5),
    ),
    (
        "ball-resting-on-an-edge-common",
        1,
        (0.707_106_781_186_547_5, 0.0, 0.707_106_781_186_547_5),
    ),
];

/// Offsets in units of the linear tolerance.
const OFFSETS: [f64; 6] = [0.5, -0.5, 2.0, -2.0, 10.0, -10.0];

/// Runs known to fail, with the follow-up that owns them: they must
/// keep failing until it lands (a target that starts passing fails the
/// test, asking for it to be removed here). The ball two tolerances
/// clear of the cube's edge serves the two shells, but validate_all
/// reads the plane x sphere circle touching the top edge from outside
/// as a crossing (PK-1153).
const KNOWN: &[(&str, f64, &str)] = &[("ball-resting-on-an-edge-union", 2.0, "PK-1153")];

/// Seconds a run may take: a hang, not a slow run.
const RUN_LIMIT_S: u64 = 300;

fn run(
    op: popkern_bool::BoolOperation,
    a: BRep,
    b: BRep,
    tol: Tolerance,
) -> Option<Result<BRep, BoolError>> {
    let (send, receive) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let options = BoolOptions {
            non_manifold: NonManifoldPolicy::Decline,
        };
        let _ = send.send(boolean_with_options(op, a, b, tol, None, &options));
    });
    receive
        .recv_timeout(std::time::Duration::from_secs(RUN_LIMIT_S))
        .ok()
}

fn area(b: &BRep, tol: &Tolerance) -> f64 {
    popkern_measure::body_surface_area(b, tol).unwrap_or(0.0)
}

/// One run's verdict.
fn judge(
    outcome: Option<Result<BRep, BoolError>>,
    expect: &Expect,
    delta: f64,
    reach: f64,
    tol: &Tolerance,
) -> Result<(), String> {
    let Some(outcome) = outcome else {
        return Err(format!("did not finish in {RUN_LIMIT_S} s"));
    };
    let within_band = delta.abs() <= tol.linear();
    match (&outcome, expect) {
        (Ok(body), _) if !body.faces.is_empty() => {
            if let Some(first) = body.validate_all().first() {
                return Err(format!("validate_all: {first}"));
            }
            let v = popkern_measure::body_volume(body, tol).map_err(|e| e.to_string())?;
            match expect {
                Expect::Body { volume, area } => {
                    let bound = if within_band {
                        1e-6 * volume + tol.linear() * area
                    } else {
                        1e-6 * volume + delta.abs() * reach
                    };
                    ((v - volume).abs() <= bound)
                        .then_some(())
                        .ok_or(format!("volume {v} for {volume} (bound {bound:.2e})"))
                }
                // The touch's answer was no body or a decline: past the
                // band a valid overlap or two parted bodies are honest,
                // within it they are not.
                _ if within_band => Err(format!("served V{v} for {expect:?}")),
                _ => Ok(()),
            }
        }
        (Ok(_), Expect::Empty) => Ok(()),
        (Ok(_), _) if !within_band => Ok(()),
        (Err(BoolError::NonManifoldResult { .. }), Expect::NonManifold(_)) => Ok(()),
        (Err(BoolError::NonManifoldResult { .. } | BoolError::RelationUndecided { .. }), _)
            if !within_band =>
        {
            Ok(())
        }
        (Ok(_), e) => Err(format!("served an empty body, expected {e:?}")),
        (Err(err), e) => Err(format!("declined `{err}`, expected {e:?}")),
    }
}

#[test]
fn touches_moved_within_and_past_the_tolerance_band() {
    let base = Tolerance::default();
    let turned = Transform3::rotation_around_axis(
        popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.7, 0.3),
        0.645_771_823_237_9,
    )
    .expect("axis");
    let all = cases();
    let mut failures = Vec::new();
    let mut known_failed = vec![false; KNOWN.len()];
    for &(name, moved, n) in SWEPT {
        let case = all.iter().find(|c| c.name == name).expect("case");
        for place in [Transform3::identity(), turned] {
            for k in OFFSETS {
                let delta = k * base.linear();
                let (a, b) = (case.build)(base);
                let shift = Transform3::translation(n.0 * delta, n.1 * delta, n.2 * delta);
                let (a, b) = if moved == 0 {
                    (a.transform(&shift).expect("move"), b)
                } else {
                    (a, b.transform(&shift).expect("move"))
                };
                let reach = area(&a, &base) + area(&b, &base);
                let (a, b) = (
                    a.transform(&place).expect("place"),
                    b.transform(&place).expect("place"),
                );
                let verdict = judge(run(case.op, a, b, base), &case.decline, delta, reach, &base);
                if let Err(why) = verdict {
                    match KNOWN.iter().position(|&(n, d, _)| n == name && d == k) {
                        Some(i) => known_failed[i] = true,
                        None => failures.push(format!("{name} delta {k} eps: {why}")),
                    }
                }
            }
        }
    }
    for (i, &(name, k, owner)) in KNOWN.iter().enumerate() {
        if !known_failed[i] {
            failures.push(format!(
                "{name} delta {k} eps ({owner}) now passes everywhere: remove it from KNOWN"
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
