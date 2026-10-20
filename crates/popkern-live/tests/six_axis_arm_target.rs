//! Target posing on the six-axis arm (PK-808, six-axis-arm §11 E4):
//! with `target_posing` on, the flange follows a target block moved by
//! tx, ty, tz while the released joints j1..j5 turn to reach it, and a
//! target past the joint limits is refused with nothing moved.

use popkern_live::{Envelope, HostOptions, Session};

const ARM: &str = include_str!("../../popkern-rhai/examples/six_axis_arm.rhai");

fn world(e: &Envelope, path: &str) -> [f32; 16] {
    e.assembly
        .as_ref()
        .expect("assembly")
        .instances
        .iter()
        .find(|i| i.path == path)
        .unwrap_or_else(|| panic!("no instance {path}"))
        .world
}

fn param(e: &Envelope, name: &str) -> f64 {
    e.params
        .iter()
        .find(|p| p.name == name)
        .unwrap_or_else(|| panic!("no param {name}"))
        .value
}

/// The translation column of a column-major world matrix.
fn origin(m: &[f32; 16]) -> [f64; 3] {
    [m[12] as f64, m[13] as f64, m[14] as f64]
}

#[test]
fn the_flange_follows_a_moved_target_within_the_joint_limits() {
    let script = ARM.replace("let target_posing = false;", "let target_posing = true;");
    assert_ne!(script, ARM, "the arm script lost its target_posing switch");
    let mut session = Session::new(HostOptions {
        max_operations: None,
        timeout: None,
        ..HostOptions::default()
    });
    let built = session.run(&script, "six_axis_arm.rhai");
    assert!(built.outcome.is_served(), "{:?}", built.outcome);

    // Move the target along x, then diagonally: the target sits at its
    // offset, every released joint stays within its limits, and the
    // arm has turned to follow.
    let limits = [
        ("j1", -360.0, 360.0),
        ("j2", -360.0, 360.0),
        ("j3", -170.0, 170.0),
        ("j4", -360.0, 360.0),
        ("j5", -360.0, 360.0),
    ];
    let mut served = 0;
    for offset in [[150.0, 0.0, 0.0], [150.0, 80.0, 0.0], [150.0, 80.0, -60.0]] {
        let mut last = None;
        for (k, name) in ["tx", "ty", "tz"].into_iter().enumerate() {
            let t = std::time::Instant::now();
            let e = session.set_param(name, offset[k]);
            let took = t.elapsed();
            assert!(e.outcome.is_served(), "{offset:?} {name}: {:?}", e.outcome);
            // A slider move: 0.1-0.25 s at the test profile once the
            // reducer followers solve alone after the released chain
            // (13-14 s as one 180-column cluster).
            assert!(
                took < std::time::Duration::from_secs(5),
                "{name} = {} took {took:?}",
                offset[k]
            );
            last = Some(e);
        }
        let e = last.expect("a move");
        let t = origin(&world(&e, "target"));
        for k in 0..3 {
            assert!(
                (t[k] - offset[k]).abs() < 1e-3,
                "{offset:?}: target at {t:?}"
            );
        }
        for (name, lo, hi) in limits {
            let q = param(&e, name);
            assert!(q >= lo - 1e-6 && q <= hi + 1e-6, "{offset:?}: {name} = {q}");
        }
        assert!(
            param(&e, "j2").abs() > 1e-3 || param(&e, "j3").abs() > 1e-3,
            "{offset:?}: the arm did not move"
        );
        served += 1;
    }
    assert_eq!(served, 3);

    // A released joint's slider is the solver's to set.
    let refused = session.set_param("j2", 10.0);
    assert!(
        !refused.outcome.is_served(),
        "a released joint took a slider value"
    );

    // Far past the arm's reach: refused, nothing moved.
    let before = session.set_param("tx", 150.0);
    let flange = world(&before, "link_6");
    let far = session.set_param("tx", -2000.0);
    assert!(!far.outcome.is_served(), "an unreachable target was served");
    let after = session.set_param("ty", 80.0);
    assert!(after.outcome.is_served(), "{:?}", after.outcome);
    let moved = world(&after, "link_6")
        .iter()
        .zip(&flange)
        .fold(0.0_f32, |m, (a, b)| m.max((a - b).abs()));
    assert!(moved < 1e-3, "a refused target moved the flange by {moved}");
}
