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

use popkern_rhai::ScriptHost;

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

fn run(body: &str, name: &str) -> Result<popkern_script::Model, String> {
    let source = format!("let m = model();\n{body}\nm\n");
    let mut host = ScriptHost::new();
    let outcome = host
        .run(&source, name)
        .map_err(|e| format!("{name}: {}", e.render()))?;
    outcome.model.ok_or_else(|| "no model".into())
}

fn union_script(c: &Casting) -> String {
    format!(
        "let d = {};\nd = d.union(m.tool(|t| {}), #{{ label: \"web\" }});\nd = d.union(m.tool(|t| {}), #{{ label: \"boss\" }});\nm.assert_watertight();",
        c.ring,
        c.web.replacen("m.", "t.", 1),
        c.boss.replacen("m.", "t.", 1)
    )
}

/// OCC's oracle on popkern's union and on its own fuse of the three
/// operands: (served volume, fused volume), or `None` without the venv.
fn occ_volumes(c: &Casting, served: &mut popkern_script::Model, name: &str) -> Option<(f64, f64)> {
    let home = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default());
    let python = home.join("venvs/ocp/bin/python");
    let tools = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle");
    if !python.exists() || !tools.exists() {
        eprintln!("[pk758] OCC oracle absent; volumes not cross-checked");
        return None;
    }
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("pk758");
    std::fs::create_dir_all(&dir).expect("dir");
    let step = dir.join(format!("{name}.step"));
    served.to_step(&step).expect("STEP export");
    let oracle = |script: &str, args: &[&std::path::Path], extra: &[&str]| -> serde_json::Value {
        let out = std::process::Command::new(&python)
            .arg(tools.join(script))
            .args(args)
            .args(extra)
            .output()
            .expect("oracle");
        serde_json::from_slice(&out.stdout).expect("oracle json")
    };
    let check = oracle("occ_check.py", &[&step], &[]);
    assert_eq!(
        check["passed"], true,
        "{name}: OCC rejects the STEP: {check}"
    );
    let solids = oracle("occ_solids.py", &[&step], &["--adaptive"]);
    let rows = solids["solids"].as_array().expect("solids");
    assert_eq!(rows.len(), 1, "{name}: one solid {solids}");
    let served_volume = rows[0]["volume"].as_f64().expect("volume");
    let mut parts = Vec::new();
    for (k, body) in [&c.ring, &c.web, &c.boss].into_iter().enumerate() {
        let mut model =
            run(&format!("let d = {body};"), &format!("{name}_{k}.rhai")).expect("operand");
        let path = dir.join(format!("{name}_operand_{k}.step"));
        model.to_step(&path).expect("operand STEP");
        parts.push(path);
    }
    let refs: Vec<&std::path::Path> = parts.iter().map(|p| p.as_path()).collect();
    let fused = oracle("occ_fuse.py", &refs, &[]);
    assert_eq!(fused["done"], true, "{name}: OCC fuse {fused}");
    Some((
        served_volume,
        fused["volume"].as_f64().expect("fused volume"),
    ))
}

/// Serve the union; with the OCC venv, its volume must match OCC's own
/// fuse of the operands within 1e-9 relative (measured: 1e-13..2e-10).
fn serves(c: &Casting, name: &str) -> Result<(), String> {
    let mut model = run(&union_script(c), &format!("{name}.rhai"))?;
    let report = popkern_rhai::report_value(&model);
    if report["failures"] != 0 {
        return Err(format!("{name}: {report}"));
    }
    if let Some((served, fused)) = occ_volumes(c, &mut model, name) {
        let rel = (served - fused).abs() / fused;
        if rel > 1e-9 {
            return Err(format!(
                "{name}: volume {served} vs OCC fuse {fused} (rel {rel:.2e})"
            ));
        }
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
    for delta in [0.0, 1e-6, -1e-3, -0.1, 0.1] {
        let name = format!("boss_{delta:e}");
        if let Err(e) = serves(&tangent_boss(-50.0 + delta, -770.0), &name) {
            failures.push(e);
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
    for delta in [-1e-6, 1e-3] {
        pinned_to_pk780(
            &tangent_boss(-50.0 + delta, -770.0),
            &format!("boss_{delta:e}"),
        );
    }
}

#[test]
fn the_near_tangent_cap_casting_serves_through_the_rim_extreme() {
    let mut failures = Vec::new();
    for cap_x in [-368.0, -367.0 + 1e-6, -367.001, -366.999, -366.9] {
        if let Err(e) = serves(&near_tangent_cap(cap_x), &format!("cap_{cap_x}")) {
            failures.push(e);
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
    for cap_x in [-367.0 - 1e-6, -367.1] {
        pinned_to_pk780(&near_tangent_cap(cap_x), &format!("cap_{cap_x}"));
    }
}
