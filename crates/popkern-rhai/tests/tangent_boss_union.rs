//! PK-758: the six-axis arm's end castings. A ring with a web through
//! it and a boss whose top generator is TANGENT to the web's top face
//! (a), or whose end cap cuts the ring wall a millimetre inside its
//! extreme (b). Both unions serve watertight, and their volume matches
//! OpenCASCADE's own fuse of the same operands (an independent boolean
//! kernel). Sweeps run the boss through exact tangency and the cap
//! through the rim extreme.
//!
//! The mechanisms, each measured on these inputs: the plane × arc
//! roots of a notch rim stored a turn away from its window (seeded
//! once per Bézier segment, one of two crossings found) now come from
//! a closed form; the web's notch through the ring wall straddles the
//! wall's u seam, so the hole filter unwraps hole contours and ray
//! classification drops hits inside curved faces' holes; and the
//! duplicate-cover test no longer calls a 0.5-thick lens one curve.
//!
//! PK-780 added four, each negative-controlled on these sweeps: the
//! pierce override votes planar faces with EXACT rim containment (an
//! 8-station octagon cut 3.2 inside the boss cap's rim, and the cap
//! 0.1-0.3 inside the ring's extreme stranded its top sliver); record
//! consolidation drops a "duplicate" only when it lies ON the other
//! record (two exact lines 0.51 apart, from a boss 1e-3 proud of the
//! web, sat inside the 1.4 slop gate); the dense record projection
//! polishes a closed record from both window ends (a junction 1.4e-3
//! short of the rim record's seam read 0.26 off it); and sliver
//! absorption checks the sliver's edge curves, not just its corners
//! (a 1e-4-proud boss strip was merged into the web plane after the
//! serve gate). OCC's own fuse is wrong 1e-6 from the ring extreme;
//! that point is judged against fuses 1e-3 and 2e-3 away.

mod union_oracle;

struct Casting {
    ring: String,
    web: String,
    boss: String,
}

/// (a): ring r 49 (y −16..−26), web y −25..−17, boss r 33 along x from
/// `cap_x` to −740 at height `boss_y` (top generator at boss_y + 33).
fn tangent_boss(boss_y: f64, cap_x: f64) -> Casting {
    Casting {
        ring: "m.profile(plane_at([-817.0, -16.0, 162.5], [0.0, -1.0, 0.0])).circle(49.0).sketch(\"ring\").extrude(blind(10.0))".into(),
        web: "m.profile(xy_at(137.5)).rect(42.95, 8.0).at(-761.475, -21.0).sketch(\"web_sketch\").extrude(blind(50.0))".into(),
        boss: format!(
            "m.profile(plane_at([{cap_x:?}, {boss_y:?}, 162.5], [1.0, 0.0, 0.0])).circle(33.0).sketch(\"boss\").extrude(blind({:?}))",
            -740.0 - cap_x
        ),
    }
}

/// (b): ring r 58 (y −2..−12), web y −11..−3, boss r 42.5 along x
/// from `cap_x` to −335 at y −45; the ring's extreme is x = −367.
fn near_tangent_cap(cap_x: f64) -> Casting {
    Casting {
        ring: "m.profile(plane_at([-425.0, -2.0, 162.5], [0.0, -1.0, 0.0])).circle(58.0).sketch(\"ring\").extrude(blind(10.0))".into(),
        web: "m.profile(xy_at(132.5)).rect(48.95, 8.0).at(-359.475, -7.0).sketch(\"web_sketch\").extrude(blind(60.0))".into(),
        boss: format!(
            "m.profile(plane_at([{cap_x:?}, -45.0, 162.5], [1.0, 0.0, 0.0])).circle(42.5).sketch(\"boss\").extrude(blind({:?}))",
            -335.0 - cap_x
        ),
    }
}

fn union_script(c: &Casting) -> String {
    format!(
        "let d = {};\nd = d.union(m.tool(|t| {}), #{{ label: \"web\" }});\nd = d.union(m.tool(|t| {}), #{{ label: \"boss\" }});\nm.assert_watertight();",
        c.ring,
        c.web.replacen("m.", "t.", 1),
        c.boss.replacen("m.", "t.", 1)
    )
}

/// Serve the union: `(served, fused)` volumes from the OCC venv (our
/// STEP's adaptive volume, OCC's own fuse of the operands), `None`
/// without it.
fn served_and_fused(c: &Casting, name: &str) -> Result<Option<(f64, f64)>, String> {
    let mut model = union_oracle::run("", &union_script(c), &format!("{name}.rhai"))?;
    let report = popkern_rhai::report_value(&model);
    if report["failures"] != 0 {
        return Err(format!("{name}: {report}"));
    }
    let operands = [c.ring.as_str(), c.web.as_str(), c.boss.as_str()];
    Ok(union_oracle::occ_volumes("", &operands, &mut model, name))
}

/// Serve the union; with the OCC venv, its volume must match OCC's own
/// fuse of the operands within 1e-9 relative (measured: 1e-13..2e-10).
fn serves(c: &Casting, name: &str) -> Result<(), String> {
    if let Some((served, fused)) = served_and_fused(c, name)? {
        let rel = (served - fused).abs() / fused;
        if rel > 1e-9 {
            return Err(format!(
                "{name}: volume {served} vs OCC fuse {fused} (rel {rel:.2e})"
            ));
        }
    }
    Ok(())
}

/// [`serves`] for a cap within 1e-5 of the ring's extreme, where OCC's
/// own fuse is not an oracle: at a cap 1e-6 past the extreme it reads
/// 135 below the fuse 1e-3 further out, while the volume is linear in
/// the cap position to ~1e-3 there (measured; the cap section exposed
/// beyond the web is ~5386, the ring lens it gains is O(√Δ)). The
/// reference is the linear extrapolation of OCC's fuses at 1e-3 and
/// 2e-3 from the extreme on the same side, within 1e-8 relative.
fn serves_near_extreme(cap_x: f64) -> Result<(), String> {
    let name = format!("cap_{cap_x}");
    let Some((served, _)) = served_and_fused(&near_tangent_cap(cap_x), &name)? else {
        return Ok(());
    };
    let side = if cap_x < -367.0 { -1.0 } else { 1.0 };
    let (x1, x2) = (-367.0 + side * 1e-3, -367.0 + side * 2e-3);
    let fused_at = |x: f64| -> Result<f64, String> {
        served_and_fused(&near_tangent_cap(x), &format!("cap_ref_{x}"))?
            .map(|(_, fused)| fused)
            .ok_or_else(|| "OCC oracle vanished".to_string())
    };
    let (f1, f2) = (fused_at(x1)?, fused_at(x2)?);
    let want = f1 + (cap_x - x1) * (f1 - f2) / (x1 - x2);
    let rel = (served - want).abs() / want;
    if rel > 1e-8 {
        return Err(format!(
            "{name}: volume {served} vs extrapolated OCC fuse {want} (rel {rel:.2e})"
        ));
    }
    Ok(())
}

/// A sweep point PK-780 owns: it must decline honestly (a boolean
/// decline, never an invalid serve); the pin flips when it serves.
fn pinned_to_pk780(c: &Casting, name: &str) {
    match serves(c, name) {
        Ok(()) => panic!("PK-780: {name} now serves; promote it to the served set"),
        Err(e) => assert!(
            e.contains("merge failed") || e.contains("result validation failed"),
            "{name}: {e}"
        ),
    }
}

#[test]
fn the_tangent_boss_casting_serves_through_tangency() {
    let mut failures = Vec::new();
    for delta in [0.0, 1e-6, 1e-4, 1e-3, 3e-3, 1e-2, -1e-3, -0.1, 0.1] {
        let name = format!("boss_{delta:e}");
        if let Err(e) = serves(&tangent_boss(-50.0 + delta, -770.0), &name) {
            failures.push(e);
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
    for delta in [-1e-6, 1e-7] {
        pinned_to_pk780(
            &tangent_boss(-50.0 + delta, -770.0),
            &format!("boss_{delta:e}"),
        );
    }
}

#[test]
fn the_near_tangent_cap_casting_serves_through_the_rim_extreme() {
    let mut failures = Vec::new();
    for cap_x in [
        -368.0,
        -367.5,
        -367.3,
        -367.1,
        -367.05,
        -367.0 + 1e-6,
        -367.001,
        -366.999,
        -366.9,
    ] {
        if let Err(e) = serves(&near_tangent_cap(cap_x), &format!("cap_{cap_x}")) {
            failures.push(e);
        }
    }
    if let Err(e) = serves_near_extreme(-367.0 - 1e-6) {
        failures.push(e);
    }
    assert!(failures.is_empty(), "{failures:#?}");
    for cap_x in [-367.0, -367.00001, -367.0001] {
        pinned_to_pk780(&near_tangent_cap(cap_x), &format!("cap_{cap_x}"));
    }
}
