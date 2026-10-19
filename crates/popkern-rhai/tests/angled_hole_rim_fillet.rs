//! The rim of a hole drilled at an angle fillets, G1 and G2 (PK-888,
//! PK-889, PK-884): the top rim of an r 5 hole tilted off a block's
//! face normal is a plane x oblique cylinder ellipse, and its contacts
//! are not circles.
//!
//! The strip is a periodic NURBS ring (PK-890): a chord-length periodic
//! cubic through the walker's stations, so every pcurve pass treats its
//! seam like a cylinder's. Before it the strip was a clamped chart that
//! overlapped its closure, every pass projected each seam use and
//! closed-rail end to the nearest of two images, and OpenCASCADE read
//! the strip's wire self-intersecting. Beneath that, each
//! negative-controlled on this sweep: the strip-side rail pcurves ran
//! their curve rather than their half-edge; a non-circular rail kept
//! the stations' winding instead of its old rim half-edge's; the rail
//! search clamped its bracket at a periodic boundary's end (0.07 off at
//! a seam 0.0017 of a turn past u = 0); a uniform-parameter strip on
//! the unevenly spaced stations of a 25 deg rim parted 3e-5 from its
//! rails; the hole's seam was shortened off its generator; the rails
//! were polylines (5e-3 off the wall); and the station refinement
//! stopped early on full-period seed windows.
//!
//! Each served body is validate_all-clean, OCC-valid as one solid, on
//! the ideal rolling ball, and matches OpenCASCADE's own fillet in
//! volume to within that fillet's error. Volumes are BRepGProp's,
//! which itself errs ~0.02 here on a cylinder face trimmed by a
//! B-spline rail: a fine-mesh volume reads the roll 0 and roll 7
//! bodies at 25 deg equal to 7.5e-5, where BRepGProp reads them 0.022
//! apart.

use popkern_rhai::ScriptHost;

/// The hole's drilling direction: `tilt_deg` off -z, rolled by
/// `roll` twelfths of a turn about z.
fn hole_axis(tilt_deg: f64, roll: usize) -> (f64, f64, f64) {
    let t = tilt_deg.to_radians();
    let phi = roll as f64 * std::f64::consts::TAU / 12.0;
    (t.sin() * phi.cos(), t.sin() * phi.sin(), -t.cos())
}

fn block(tilt_deg: f64, roll: usize, fillet: Option<&str>) -> String {
    let (ax, ay, az) = hole_axis(tilt_deg, roll);
    let rim = fillet
        .map(|c| {
            format!(
                "b = b.fillet(edges().between(\"plane\", \"cylinder\").at_height(20.0), 1.5, \
                 #{{ continuity: \"{c}\", label: \"rim\" }});\n"
            )
        })
        .unwrap_or_default();
    format!(
        "let m = model();\n\
         let b = m.profile(xy()).rect(40.0, 40.0).sketch(\"s\").extrude(blind(20.0));\n\
         let d = m.tool(|t| t.profile(plane_at([0.0, 0.0, 30.0], [{ax}, {ay}, {az}])).circle(5.0).sketch(\"h\").extrude(blind(60.0)));\n\
         b = b.subtract(d, #{{ label: \"hole\" }});\n{rim}m"
    )
}

fn model(source: &str) -> Result<popkern_script::Model, String> {
    ScriptHost::new()
        .run(source, "angled_hole.rhai")
        .map_err(|e| format!("{} {}", e.render(), e.payload))?
        .model
        .ok_or_else(|| "no model".to_string())
}

/// The OCC oracle's JSON for `script` run on `model`'s STEP; `None`
/// without the OCP venv.
fn oracle(
    model: &mut popkern_script::Model,
    name: &str,
    script: &str,
    args: &[&str],
) -> Option<serde_json::Value> {
    let home = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default());
    let python = home.join("venvs/ocp/bin/python");
    let tools = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle");
    if !python.exists() {
        eprintln!("[pk888] OCC oracle absent; not cross-checked");
        return None;
    }
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("pk888");
    std::fs::create_dir_all(&dir).expect("dir");
    let step = dir.join(format!("{name}.step"));
    model.to_step(&step).expect("STEP export");
    let out = std::process::Command::new(&python)
        .arg(tools.join(script))
        .arg(&step)
        .args(args)
        .output()
        .expect("oracle");
    Some(serde_json::from_slice(&out.stdout).expect("oracle json"))
}

/// The worst rolling-ball deviation of the served fillet strip: from
/// each strip point a ball centre one radius along the normal must sit
/// 1.5 from the top plane z = 20 and from the r 5 hole wall. Measured
/// against the analytic supports, independently of the construction.
fn ball_deviation(brep: &popkern_core::topology::BRep, tilt: f64, roll: usize) -> f64 {
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::geometry::{Point3, SurfaceEval, Vec3};
    let (ax, ay, az) = hole_axis(tilt, roll);
    let axis = Vec3::new(ax, ay, az);
    let origin = Point3::new(0.0, 0.0, 30.0);
    let mut worst: f64 = 0.0;
    let mut strips = 0;
    for (_, face) in brep.faces.iter() {
        let Some(SurfaceGeometry::Nurbs(strip)) = face.surface.as_ref().map(|s| &s.geometry) else {
            continue;
        };
        strips += 1;
        let ((u0, u1), (v0, v1)) = strip.parameter_range();
        for i in 0..400 {
            for j in 0..=8 {
                let u = u0 + (u1 - u0) * i as f64 / 400.0;
                let v = v0 + (v1 - v0) * j as f64 / 8.0;
                let p = strip.point_at(u, v);
                let Some(n) = strip.normal_at(u, v).normalized() else {
                    continue;
                };
                for sign in [1.0, -1.0] {
                    let c = p + n * (1.5 * sign);
                    let to_plane = 20.0 - c.z;
                    let w = c - origin;
                    let to_wall = (w - axis * w.dot(&axis)).length() - 5.0;
                    if to_plane > 0.0 && to_wall > 0.0 {
                        worst = worst.max((to_plane - 1.5).abs()).max((to_wall - 1.5).abs());
                    }
                }
            }
        }
    }
    assert_eq!(strips, 1, "one fillet strip");
    worst
}

/// Every placement serves validate_all-clean, OCC-valid as one solid,
/// on the ideal rolling ball to the linear tolerance, and agrees in
/// volume with OpenCASCADE's own fillet of the same block built
/// natively, to within that fillet's own measured error: OCC
/// approximates a non-circular fillet (9.8e-4 off the ball at 8 deg,
/// 8.2e-3 at 25 deg, exact at 0 deg) and its volume errs by about that
/// times the fillet's area.
#[test]
fn an_angled_holes_rim_fillets_on_the_rolling_ball() {
    for tilt in [8.0, 25.0, 30.0] {
        for roll in [0, 7] {
            let name = format!("tilt{tilt}_roll{roll}");
            let mut filleted =
                model(&block(tilt, roll, Some("g1"))).unwrap_or_else(|e| panic!("{name}: {e}"));
            let brep = filleted.to_brep().expect("brep");
            let errors = brep.validate_all();
            assert!(errors.is_empty(), "{name}: {errors:?}");
            let ball = ball_deviation(&brep, tilt, roll);
            assert!(
                ball <= 1e-5,
                "{name}: strip {ball:.2e} off the rolling ball"
            );
            let Some(check) = oracle(&mut filleted, &name, "occ_check.py", &[]) else {
                continue;
            };
            assert_eq!(
                check["passed"], true,
                "{name}: OCC rejects the STEP: {check}"
            );
            let solids =
                oracle(&mut filleted, &name, "occ_solids.py", &["--adaptive"]).expect("oracle");
            let solids = solids["solids"].as_array().expect("solids");
            assert_eq!(solids.len(), 1, "{name}: one solid");
            let served = solids[0]["volume"].as_f64().expect("volume");
            let (ax, ay, az) = hole_axis(tilt, roll);
            let axis = [ax, ay, az].map(|c| c.to_string());
            let mut args: Vec<&str> = vec!["--angled-hole"];
            args.extend(axis.iter().map(String::as_str));
            args.extend(["20", "1.5"]);
            let occ = oracle(&mut filleted, &name, "occ_fillet.py", &args).expect("oracle");
            assert_eq!(occ["valid"], true, "{name}: OCC's own fillet: {occ}");
            let want = occ["volume"].as_f64().expect("OCC fillet volume");
            let occ_error = occ["ball_worst"].as_f64().expect("ball")
                * occ["fillet_area"].as_f64().expect("area");
            let gap = (served - want).abs();
            eprintln!(
                "[pk888] {name}: ball {ball:.2e}; served {served} vs OCC fillet {want} \
                 (gap {gap:.3e}, OCC's own error bound {occ_error:.3e})"
            );
            assert!(
                gap <= occ_error + 1e-9 * want,
                "{name}: volume {served} vs OCC fillet {want}: gap {gap:.3e} past OCC's \
                 error bound {occ_error:.3e}"
            );
        }
    }
}

/// The G2 rim (PK-884): the strip is the periodic ring through the
/// walker's G2 cross-sections, its stations refined like G1's (a closed
/// spine has no ends to extend). Every placement serves
/// validate_all-clean and OCC-valid as one solid, and popkern's own
/// continuity evidence reads G2 on every strip-to-support seam.
#[test]
fn an_angled_holes_rim_fillets_g2() {
    use popkern_blend::{ContinuityGates, ContinuityLevel, measure_blend_continuity_evidence};
    use popkern_core::geometry::surface::SurfaceGeometry;
    for tilt in [8.0, 25.0] {
        for roll in [0, 7] {
            let name = format!("g2_tilt{tilt}_roll{roll}");
            let mut filleted =
                model(&block(tilt, roll, Some("g2"))).unwrap_or_else(|e| panic!("{name}: {e}"));
            let brep = filleted.to_brep().expect("brep");
            let errors = brep.validate_all();
            assert!(errors.is_empty(), "{name}: {errors:?}");
            let strips: Vec<_> = brep
                .faces
                .iter()
                .filter(|(_, f)| {
                    matches!(
                        f.surface.as_ref().map(|s| &s.geometry),
                        Some(SurfaceGeometry::Nurbs(_))
                    )
                })
                .map(|(h, _)| h)
                .collect();
            assert_eq!(strips.len(), 1, "{name}: one fillet strip");
            let evidence = measure_blend_continuity_evidence(
                &brep,
                &strips,
                ContinuityLevel::G2,
                &ContinuityGates::commercial(),
                33,
            );
            assert_eq!(
                evidence.achieved,
                ContinuityLevel::G2,
                "{name}: {evidence:?}"
            );
            eprintln!(
                "[pk884] {name}: g1 {:.2e} g2 {:.2e}",
                evidence.worst_g1_normal, evidence.worst_g2_residual
            );
            let Some(check) = oracle(&mut filleted, &name, "occ_check.py", &[]) else {
                continue;
            };
            assert_eq!(
                check["passed"], true,
                "{name}: OCC rejects the STEP: {check}"
            );
            let solids =
                oracle(&mut filleted, &name, "occ_solids.py", &["--adaptive"]).expect("oracle");
            let solids = solids["solids"].as_array().expect("solids");
            assert_eq!(solids.len(), 1, "{name}: one solid");
        }
    }
}
