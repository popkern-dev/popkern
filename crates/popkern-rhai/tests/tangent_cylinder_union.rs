//! Unions where a box face is tangent to a cylinder along a line
//! (PK-750). The six-axis arm's base slots are stadiums: two rods with a
//! box between them whose long faces are tangent to both rods along the
//! box's own vertical edges, and whose caps are flush with the rods'.
//! Three defects stacked on this input: the rod's far rim arc was
//! removed as a "duplicate" of the near-rim record (same end nodes,
//! same carrier circle), the box's tangent face sampled its own tangent
//! edge (on the rod) and classified Inside, and the rod's seam, replaced
//! by a cut line lying on it, lost its seam role in the chart.

use popkern_rhai::ScriptHost;

const PI: f64 = std::f64::consts::PI;
const R: f64 = 4.5;
const H: f64 = 16.0;

const ROD: &str = r#"
fn rod_at(t, p, n, r, length, label) {
    t.profile(plane_at(p, n)).circle(r).sketch(label).extrude(blind(length))
}
"#;

/// A served union: meshed volume, face count, and the model.
struct Served {
    volume: f64,
    faces: usize,
    model: popkern_script::Model,
}

fn serve(body: &str, name: &str) -> Result<Served, String> {
    let source = format!("{ROD}\nlet m = model();\n{body}\nm.assert_watertight();\nm\n");
    let mut host = ScriptHost::new();
    let outcome = host
        .run(&source, name)
        .map_err(|e| format!("{name}: {}\n{}", e.render(), e.payload))?;
    let mut model = outcome.model.ok_or("no model")?;
    let report = popkern_rhai::report_value(&model);
    assert_eq!(report["failures"], 0, "{name}: {report}");
    for row in report["features"].as_array().into_iter().flatten() {
        assert_eq!(
            row["fallbacks"].as_array().map(Vec::len).unwrap_or(0),
            0,
            "{name}: a fallback served {}: {row}",
            row["label"]
        );
    }
    let part = model.to_part().map_err(|e| e.to_string())?;
    let brep = part.brep();
    brep.validate()
        .map_err(|e| format!("{name}: invalid: {e}"))?;
    let volume = popkern_measure::body_volume(brep, &brep.tolerance).map_err(|e| e.to_string())?;
    let faces = brep.faces.len();
    Ok(Served {
        volume,
        faces,
        model,
    })
}

/// The meshed volume within 1e-4 (the kernel meshes curved faces), and
/// OCC's adaptive volume of the STEP within 1e-9 when the OCP venv is
/// present.
fn assert_volume(name: &str, served: &mut Served, want: f64, solids: usize) {
    assert_volume_within(name, served, want, solids, 1e-9 * want);
}

/// As [`assert_volume`], with OCC's volume within `allowed` (absolute).
fn assert_volume_within(name: &str, served: &mut Served, want: f64, solids: usize, allowed: f64) {
    let rel = (served.volume - want).abs() / want;
    assert!(
        rel < 1e-4,
        "{name}: volume {} vs {want} (rel {rel:.2e})",
        served.volume
    );
    let home = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default());
    let python = home.join("venvs/ocp/bin/python");
    let tools = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle");
    if !python.exists() || !tools.exists() {
        eprintln!("[pk750] OCC oracle absent; exact volume not cross-checked");
        return;
    }
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("pk750");
    std::fs::create_dir_all(&dir).expect("dir");
    let step = dir.join(format!("{name}.step"));
    served.model.to_step(&step).expect("STEP export");
    let run = |script: &str, extra: &[&str]| -> serde_json::Value {
        let out = std::process::Command::new(&python)
            .arg(tools.join(script))
            .arg(&step)
            .args(extra)
            .output()
            .expect("oracle");
        serde_json::from_slice(&out.stdout).expect("oracle json")
    };
    let check = run("occ_check.py", &[]);
    assert_eq!(
        check["passed"], true,
        "{name}: OCC rejects the STEP: {check}"
    );
    let rows = run("occ_solids.py", &["--adaptive"]);
    let rows = rows["solids"].as_array().expect("solids");
    assert_eq!(rows.len(), solids, "{name}: solids {rows:?}");
    let occ: f64 = rows.iter().filter_map(|r| r["volume"].as_f64()).sum();
    let err = (occ - want).abs();
    assert!(
        err <= allowed,
        "{name}: OCC volume {occ} vs {want} (|err| {err:.2e} > {allowed:.2e})"
    );
}

/// The stadium's cross-section area: a 12 × 2w box between two discs of
/// radius R centred on its short sides. Below tangency (w < R) each
/// disc's inner half also pokes out past the box by one circular
/// segment beyond |y| = w.
fn stadium_area(w: f64) -> f64 {
    let segment = if w < R {
        R * R * (w / R).acos() - w * (R * R - w * w).sqrt()
    } else {
        0.0
    };
    12.0 * 2.0 * w + PI * R * R + 2.0 * segment
}

fn stadium(w: f64) -> String {
    format!(
        r#"
let a = rod_at(m, [-6.0, 0.0, 0.0], [0.0, 0.0, 1.0], {R:?}, {H:?}, "a");
let b = m.tool(|t| rod_at(t, [6.0, 0.0, 0.0], [0.0, 0.0, 1.0], {R:?}, {H:?}, "b"));
a = a.union(b, #{{ label: "ab" }});
let c = m.tool(|t| t.profile(xy_at(0.0)).rect(12.0, {:?}).sketch("mid").extrude(blind({H:?})));
a.union(c, #{{ label: "stadium" }});
"#,
        2.0 * w
    )
}

#[test]
fn a_stadium_box_tangent_to_two_rods_serves() {
    let mut s = serve(&stadium(R), "stadium.rhai").unwrap_or_else(|e| panic!("{e}"));
    // Top, bottom, two flat sides, two half walls.
    assert_eq!(s.faces, 6, "stadium faces");
    assert_volume("stadium", &mut s, stadium_area(R) * H, 1);
}

/// The box half-width swept through exact tangency: each side serves
/// with the analytic volume, or the sweep names the offset that fails.
/// Within the linear tolerance of tangency (|δ| = 1e-6 at the default
/// 1e-6) the contact is ambiguous: the box's corner edge lies within
/// tolerance of the rod and is welded onto the rod's boundary (measured:
/// δ = −1e-6 serves the rod wall through the box corner 1e-6 inside the
/// rod). Every boundary point then sits within tolerance of the exact
/// body, so |ΔV| ≤ tolerance × surface area bounds the volume.
#[test]
fn the_stadium_serves_through_tangency() {
    let linear = popkern_core::tolerance::Tolerance::default().linear();
    let mut failures = Vec::new();
    for delta in [-0.1, -1e-3, -1e-6, 1e-6, 1e-3, 0.1] {
        // δ = −1e-3: the box face sits 1e-3 inside tangency over a
        // coplanar cap overlap; a rim arc dipping 1e-3 below the box
        // face lands on it and validation declines honestly (the
        // near-tangency class, PK-780); the pin flips when that lands.
        if delta == -1e-3 {
            let Err(e) = serve(&stadium(R + delta), "stadium_near.rhai") else {
                panic!("PK-780: the near-tangent stadium now serves; promote it");
            };
            assert!(e.contains("result validation failed"), "{e}");
            continue;
        }
        let w = R + delta;
        let name = format!("stadium_{delta:e}");
        let want = stadium_area(w) * H;
        let allowed = if delta.abs() <= linear * (1.0 + 1e-9) {
            let surface = 2.0 * stadium_area(w) + 2.0 * 12.0 * H + 2.0 * PI * R * H;
            linear * surface
        } else {
            1e-9 * want
        };
        match serve(&stadium(w), &format!("{name}.rhai")) {
            Ok(mut s) => assert_volume_within(&name, &mut s, want, 1, allowed),
            Err(e) => failures.push(format!("δ = {delta:e}: {}", e.lines().next().unwrap_or(""))),
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

#[test]
fn a_box_tangent_to_one_rod_serves() {
    let body = format!(
        r#"
let a = rod_at(m, [-6.0, 0.0, 0.0], [0.0, 0.0, 1.0], {R:?}, {H:?}, "a");
let c = m.tool(|t| t.profile(xy_at(0.0)).rect(12.0, 9.0).sketch("mid").extrude(blind({H:?})));
a.union(c, #{{ label: "half" }});
"#
    );
    let mut s = serve(&body, "half.rhai").unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(s.faces, 6, "half stadium faces");
    let want = (12.0 * 2.0 * R + 0.5 * PI * R * R) * H;
    assert_volume("half", &mut s, want, 1);
}

/// A box touching a rod along one line from outside: the union is two
/// lumps sharing a line, or a declared refusal — never an invalid serve.
#[test]
fn a_box_touching_a_rod_along_a_line_never_serves_invalid() {
    let body = format!(
        r#"
let a = rod_at(m, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], {R:?}, {H:?}, "a");
let c = m.tool(|t| t.profile(xy_at(0.0)).rect(6.0, 9.0).sketch("side").extrude(blind({H:?})).translate(7.5, 0.0, 0.0));
a.union(c, #{{ label: "touch" }});
"#
    );
    match serve(&body, "touch.rhai") {
        Ok(mut s) => {
            assert_eq!(s.faces, 9, "a rod and a box, untouched");
            let want = (6.0 * 9.0 + PI * R * R) * H;
            assert_volume("touch", &mut s, want, 2);
        }
        Err(e) => assert!(
            !e.contains("invalid:"),
            "a touching union served an invalid body: {e}"
        ),
    }
}
