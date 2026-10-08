//! Mate-graph decomposition and driven continuation (six-axis-arm
//! §5.3): a serial chain solves link by link from the ground, a loop
//! solves as one cluster, and a large slider jump lands where small
//! steps do.

use popkern_assembly::{
    Assembly, Datum, MateEntity, MateKind, MateRef, Parent, Placement, RevoluteOptions,
};
use popkern_script::Plane;
use popkern_script::prelude::*;

fn bar() -> Model {
    let mut m = Model::new();
    m.profile(Plane::XY)
        .rect(2.0, 0.4)
        .sketch("bar")
        .expect("sketch")
        .extrude(&mut m, End::Blind(0.2.into()))
        .expect("extrude");
    m
}

fn axis(instance: &str, x: f64) -> MateRef {
    MateRef {
        instance: instance.into(),
        entity: MateEntity::Axis {
            origin: [x, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
        },
    }
}

fn datum(instance: &str, datum: Datum) -> MateRef {
    MateRef {
        instance: instance.into(),
        entity: MateEntity::Datum { datum },
    }
}

/// `links` bars in a serial chain: link k turns about link k−1's tip
/// (x = 1) with joint `j<k>`, its own pivot at x = −1.
fn chain(links: usize) -> Assembly {
    let mut asm = Assembly::new("chain");
    let d = asm.add_part("bar", bar()).expect("part");
    asm.place(Parent::Root, d, Some("link_0"), Placement::identity())
        .expect("base");
    for k in 1..=links {
        asm.place(
            Parent::Root,
            d,
            Some(&format!("link_{k}")),
            Placement::translation(2.0 * k as f64 + 0.3, 0.5, 0.0),
        )
        .expect("link");
        asm.add_revolute(
            &format!("j{k}"),
            axis(&format!("link_{}", k - 1), 1.0),
            axis(&format!("link_{k}"), -1.0),
            &RevoluteOptions {
                zero: Some((
                    datum(&format!("link_{}", k - 1), Datum::X),
                    datum(&format!("link_{k}"), Datum::X),
                )),
                angle: Some(0.0),
                ..RevoluteOptions::default()
            },
        )
        .expect("joint");
    }
    asm
}

/// The tip of the last link (its x = 1 point), world.
fn tip(asm: &Assembly, links: usize) -> [f64; 3] {
    let row = asm
        .leaves()
        .into_iter()
        .find(|r| r.path.to_string() == format!("link_{links}"))
        .expect("tip link");
    let o = row.world.origin();
    let x = row.world.x_axis();
    [o.x + x.x, o.y + x.y, o.z + x.z]
}

/// Forward kinematics: every link is length 2, turned by the sum of the
/// joint angles before it.
fn tip_fk(angles: &[f64]) -> [f64; 3] {
    let (mut x, mut y, mut turn) = (1.0, 0.0, 0.0_f64);
    for a in angles {
        turn += a.to_radians();
        x += 2.0 * turn.cos();
        y += 2.0 * turn.sin();
    }
    [x, y, 0.0]
}

#[test]
fn a_serial_chain_solves_one_link_at_a_time() {
    let links = 31; // 32 roots with the base
    let mut asm = chain(links);
    let report = asm.solve_mates().expect("solve");
    assert_eq!(report.clusters.len(), links, "one cluster per link");
    assert!(report.clusters.iter().all(|c| c.roots.len() == 1));
    let mut angles = vec![0.0; links];
    for (k, a) in angles.iter_mut().enumerate() {
        *a = 7.0 * ((k % 5) as f64 - 2.0);
        asm.set_assembly_param(&format!("j{}", k + 1), *a)
            .expect("slider");
    }
    let (t, f) = (tip(&asm, links), tip_fk(&angles));
    for i in 0..3 {
        assert!((t[i] - f[i]).abs() < 1e-9, "tip {t:?} vs FK {f:?}");
    }
    // A slider move re-solves the whole chain fast: measured 2.7 ms at
    // the test profile (opt-level 2; budget 5 ms native, design §8). The
    // bound leaves headroom for a loaded machine.
    let start = std::time::Instant::now();
    asm.set_assembly_param("j1", 25.0).expect("slider");
    let ms = start.elapsed().as_secs_f64() * 1e3;
    assert!(ms < 50.0, "a 32-root chain slider re-solve took {ms:.1} ms");
}

#[test]
fn a_loop_solves_as_one_cluster() {
    // Two struts, one pinned at the base tip (1, 0), one at link_1's tip
    // (3, 0), their own tips joined: a triangle of sides 2, 2, 2 that no
    // single strut's mates to placed links fix alone.
    let mut asm = chain(1);
    let d = asm.definition_named("bar").expect("bar");
    for (name, x) in [("strut_a", 0.0), ("strut_b", 3.0)] {
        asm.place(
            Parent::Root,
            d,
            Some(name),
            Placement::translation(x, 3.0, 0.0),
        )
        .expect("strut");
    }
    let concentric = |a: MateRef, b: MateRef| popkern_assembly::Mate {
        kind: MateKind::Concentric,
        a,
        b,
        drive: None,
    };
    asm.add_mate(concentric(axis("link_0", 1.0), axis("strut_a", -1.0)))
        .expect("pivot a");
    asm.add_mate(concentric(axis("link_1", 1.0), axis("strut_b", -1.0)))
        .expect("pivot b");
    let report = asm
        .add_mate(concentric(axis("strut_a", 1.0), axis("strut_b", 1.0)))
        .expect("loop closure");
    assert!(report.converged);
    assert!(
        report
            .clusters
            .iter()
            .any(|c| c.roots.contains(&"strut_a".to_string())
                && c.roots.contains(&"strut_b".to_string())),
        "the loop solves together: {:?}",
        report.clusters
    );
}

#[test]
fn a_large_slider_jump_lands_where_small_steps_do() {
    let mut jump = chain(3);
    let mut stepped = chain(3);
    for (j, target) in [("j1", 179.0), ("j2", -170.0), ("j3", 135.0)] {
        jump.set_assembly_param(j, target).expect("jump");
        let mut v = 0.0;
        while (target - v).abs() > 1e-12 {
            v += (target - v).clamp(-5.0, 5.0);
            stepped.set_assembly_param(j, v).expect("step");
        }
    }
    let (a, b) = (tip(&jump, 3), tip(&stepped, 3));
    for i in 0..3 {
        assert!((a[i] - b[i]).abs() < 1e-9, "jump {a:?} vs stepped {b:?}");
    }
}

/// The triangle loop of `a_loop_solves_as_one_cluster`, driven by j1.
fn triangle() -> Assembly {
    let mut asm = chain(1);
    let d = asm.definition_named("bar").expect("bar");
    for (name, x) in [("strut_a", 0.0), ("strut_b", 3.0)] {
        asm.place(
            Parent::Root,
            d,
            Some(name),
            Placement::translation(x, 3.0, 0.0),
        )
        .expect("strut");
    }
    let concentric = |a: MateRef, b: MateRef| popkern_assembly::Mate {
        kind: MateKind::Concentric,
        a,
        b,
        drive: None,
    };
    for (a, b) in [
        (axis("link_0", 1.0), axis("strut_a", -1.0)),
        (axis("link_1", 1.0), axis("strut_b", -1.0)),
        (axis("strut_a", 1.0), axis("strut_b", 1.0)),
    ] {
        asm.add_mate(concentric(a, b)).expect("loop mate");
    }
    asm
}

fn apex(asm: &Assembly) -> [f64; 3] {
    let row = asm
        .leaves()
        .into_iter()
        .find(|r| r.path.to_string() == "strut_a")
        .expect("strut");
    let (o, x) = (row.world.origin(), row.world.x_axis());
    [o.x + x.x, o.y + x.y, o.z + x.z]
}

#[test]
fn a_driven_loop_keeps_its_branch_through_a_large_jump() {
    let mut jump = triangle();
    let mut stepped = triangle();
    let start = apex(&jump);
    jump.set_assembly_param("j1", 150.0).expect("jump");
    let mut v = 0.0;
    while v < 150.0 {
        v = (v + 5.0_f64).min(150.0);
        stepped.set_assembly_param("j1", v).expect("step");
    }
    let (a, b) = (apex(&jump), apex(&stepped));
    for i in 0..3 {
        assert!(
            (a[i] - b[i]).abs() < 1e-9,
            "the jump took another branch: {a:?} vs stepped {b:?} (start {start:?})"
        );
    }
}
