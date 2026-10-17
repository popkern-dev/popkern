//! A film frame's pose in one request (PK-853): `set_joints` places the
//! six-axis arm exactly where the same joints set one by one place it,
//! with one solve and one envelope — the per-frame cost a film export
//! pays. A non-joint name falls back to setting each in turn.

use popkern_live::{Envelope, HostOptions, Session};

const ARM: &str = include_str!("../../popkern-rhai/examples/six_axis_arm.rhai");

fn placed(e: &Envelope) -> Vec<(String, [f32; 16])> {
    e.assembly
        .as_ref()
        .expect("assembly")
        .instances
        .iter()
        .map(|i| (i.path.clone(), i.world))
        .collect()
}

fn pose(k: f64) -> Vec<(String, f64)> {
    let mut values: Vec<(String, f64)> = (1..=6)
        .map(|j| (format!("j{j}"), 12.0 * k + 5.0 * j as f64))
        .collect();
    values.push(("grip".to_string(), 10.0 + 4.0 * k));
    values
}

#[test]
fn a_pose_in_one_request_places_the_arm_as_its_joints_one_by_one() {
    let mut session = Session::new(HostOptions {
        max_operations: None,
        timeout: None,
        ..HostOptions::default()
    });
    assert!(session.run(ARM, "six_axis_arm.rhai").outcome.is_served());

    // One by one, then away, then the same pose in one request.
    let mut sequential = None;
    for (name, value) in pose(2.0) {
        sequential = Some(session.set_param(&name, value));
    }
    let sequential = placed(&sequential.expect("a pose"));
    for (name, value) in pose(0.0) {
        assert!(session.set_param(&name, value).outcome.is_served());
    }
    let batched = session.set_joints(&pose(2.0));
    assert!(batched.outcome.is_served(), "{:?}", batched.outcome);
    assert!(batched.bodies_unchanged && batched.bodies.is_empty());
    let rows: Vec<_> = batched
        .params
        .iter()
        .map(|p| (p.name.clone(), p.value))
        .collect();
    for (name, value) in pose(2.0) {
        assert!(
            rows.iter()
                .any(|(n, v)| *n == name && (v - value).abs() < 1e-9),
            "{name} row"
        );
    }
    let mut worst = 0.0_f32;
    for ((path, a), (other, b)) in sequential.iter().zip(placed(&batched)) {
        assert_eq!(*path, other);
        for (x, y) in a.iter().zip(&b) {
            worst = worst.max((x - y).abs());
        }
    }
    // f32 world matrices in model units (mm): 1e-3 is a micron.
    assert!(
        worst < 1e-3,
        "pose differs from its joints one by one by {worst}"
    );

    // A frame's pose costs one solve: measured against seven
    // `set_param`s of the same moves (63 ms native at the test
    // profile, 7 × 270 kB envelopes) it must be well under half.
    let time = |f: &mut dyn FnMut(f64)| {
        let mut best = std::time::Duration::MAX;
        for k in 3..6 {
            let t = std::time::Instant::now();
            f(k as f64);
            best = best.min(t.elapsed());
        }
        best
    };
    let one_by_one = time(&mut |k| {
        for (name, value) in pose(k) {
            let _ = session.set_param(&name, value).to_json();
        }
    });
    let one_request = time(&mut |k| {
        let _ = session.set_joints(&pose(k + 0.5)).to_json();
    });
    eprintln!("[pk853] pose: one by one {one_by_one:?}, one request {one_request:?}");
    assert!(
        one_request.as_secs_f64() < 0.5 * one_by_one.as_secs_f64(),
        "one request {one_request:?} vs one by one {one_by_one:?}"
    );

    // A part parameter in the pose falls back to setting each in turn.
    let mut mixed = pose(1.0);
    mixed.push(("not_a_joint".to_string(), 1.0));
    let declined = session.set_joints(&mixed);
    assert!(!declined.outcome.is_served(), "an unknown name served");
}
