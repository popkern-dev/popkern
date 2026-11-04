//! Sheet × solid booleans through the cell complex against
//! OpenCASCADE (PK-978 design 9): Common(solid, sheet) area,
//! Cut(sheet, solid) area and Splitter(solid | sheet) piece volumes
//! over the design's corpus, each in four placements (rigid motions
//! and scales 10 and 0.01). `occ_sheet.py` builds the same operands
//! from their parameters, so the kernels share neither boolean code nor
//! an exchange format. Every row must be served by the complex (no
//! `PartitionFallback` line).
//!
//! The OCC leg SKIPS when the OCP venv is absent; the popkern results
//! are still validated.

use popkern_bool::{BoolOperation, SheetBooleanResult, boolean_with_sheets};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_measure::{body_volume, face_area};
use serde_json::{Value, json};

/// Similarity placement: scale about the origin, rotate, translate.
struct Placement {
    axis: [f64; 3],
    angle: f64,
    translation: [f64; 3],
    scale: f64,
}

impl Placement {
    fn transform(&self) -> Transform3 {
        let [ax, ay, az] = self.axis;
        let [tx, ty, tz] = self.translation;
        let rotation = Transform3::rotation_around_axis(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(ax, ay, az),
            self.angle,
        )
        .expect("rotation");
        Transform3::translation(tx, ty, tz)
            .compose(&rotation.compose(&Transform3::scale_uniform(self.scale)))
    }

    fn json(&self) -> Value {
        json!({"axis": self.axis, "angle": self.angle, "translation": self.translation, "scale": self.scale})
    }
}

const PLACEMENTS: [Placement; 4] = [
    Placement {
        axis: [0.0, 0.0, 1.0],
        angle: 0.0,
        translation: [0.0, 0.0, 0.0],
        scale: 1.0,
    },
    Placement {
        axis: [1.0, 1.0, 0.0],
        angle: 0.7,
        translation: [3.0, -2.0, 1.0],
        scale: 1.0,
    },
    Placement {
        axis: [0.0, 1.0, 1.0],
        angle: 2.1,
        translation: [-1.0, 4.0, 2.0],
        scale: 10.0,
    },
    Placement {
        axis: [1.0, 0.0, 1.0],
        angle: 1.3,
        translation: [0.02, 0.0, -0.01],
        scale: 0.01,
    },
];

fn solid(spec: &Value, tol: Tolerance) -> BRep {
    match spec["kind"].as_str() {
        Some("box") => BRep::unit_cube(tol).expect("cube"),
        Some("sphere") => BRep::sphere(spec["radius"].as_f64().expect("r"), tol).expect("sphere"),
        other => panic!("solid {other:?}"),
    }
}

fn sheet(spec: &Value, tol: Tolerance) -> BRep {
    let f = |k: &str, i: usize| spec[k][i].as_f64().expect("number");
    let body = match spec["kind"].as_str() {
        Some("plane") => BRep::sheet_quad(
            f("x", 0),
            f("x", 1),
            f("y", 0),
            f("y", 1),
            spec["z"].as_f64().expect("z"),
            tol,
        )
        .expect("quad"),
        Some("tube") => {
            let (z0, z1) = (f("z", 0), f("z", 1));
            let r = spec["radius"].as_f64().expect("r");
            let mut c = BRep::cylinder(r, z1 - z0, tol)
                .and_then(|c| {
                    c.transform(&Transform3::translation(f("centre", 0), f("centre", 1), z0))
                })
                .expect("cylinder");
            let caps: Vec<_> = c
                .faces
                .iter()
                .filter(|(_, face)| {
                    matches!(
                        face.surface.as_ref().map(|s| &s.geometry),
                        Some(popkern_core::geometry::surface::SurfaceGeometry::Plane { .. })
                    )
                })
                .map(|(h, _)| h)
                .collect();
            c.retire_faces_to_boundary(&caps).expect("open tube");
            c
        }
        other => panic!("sheet {other:?}"),
    };
    let tilt = spec["tilt_x"].as_f64().unwrap_or(0.0);
    body.transform(&Transform3::rotation_x(tilt)).expect("tilt")
}

fn area(b: &BRep, tol: &Tolerance) -> f64 {
    b.faces
        .keys()
        .map(|f| face_area(b, f, tol).expect("area"))
        .sum()
}

fn open_area(
    r: Result<SheetBooleanResult, popkern_bool::BoolError>,
    row: &str,
    tol: &Tolerance,
) -> f64 {
    match r.unwrap_or_else(|e| panic!("{row}: {e}")) {
        SheetBooleanResult::OpenShell(b) => {
            let findings = b.validate_all();
            assert!(findings.is_empty(), "{row}: {findings:?}");
            area(&b, tol)
        }
        other => panic!("{row}: expected an open shell, got {other:?}"),
    }
}

fn volumes(
    r: Result<SheetBooleanResult, popkern_bool::BoolError>,
    row: &str,
    tol: &Tolerance,
) -> Vec<f64> {
    let pieces = match r.unwrap_or_else(|e| panic!("{row}: {e}")) {
        SheetBooleanResult::Solid(b) => vec![b],
        SheetBooleanResult::Pieces(p) => p,
        other => panic!("{row}: expected solids, got {other:?}"),
    };
    let mut out: Vec<f64> = pieces
        .iter()
        .map(|p| {
            let findings = p.validate_all();
            assert!(findings.is_empty(), "{row}: {findings:?}");
            body_volume(p, tol).expect("volume")
        })
        .collect();
    out.sort_by(f64::total_cmp);
    out
}

fn oracle(spec: &Value) -> Option<Value> {
    let python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let script = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tools/occ-oracle/occ_sheet.py");
    if !python.exists() || !script.exists() {
        return None;
    }
    let out = std::process::Command::new(&python)
        .arg(&script)
        .arg(spec.to_string())
        .output()
        .expect("oracle");
    Some(serde_json::from_slice(&out.stdout).expect("oracle json"))
}

/// Relative difference, against the larger of the two and a floor of
/// `floor` (an empty part is 0 in both kernels).
fn gap(a: f64, b: f64, floor: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs()).max(floor)
}

/// The corpus runs in a child of this test binary with the shadow log
/// set (the complex records every fallback there), so the parent can
/// assert none happened without touching its own environment.
#[test]
fn sheet_booleans_match_occ() {
    const LOG: &str = "POPKERN_COMPLEX_SHADOW_LOG";
    if std::env::var_os(LOG).is_some() {
        return corpus();
    }
    let log = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("sheet_occ_oracle.log");
    let _ = std::fs::remove_file(&log);
    let child = std::process::Command::new(std::env::current_exe().expect("test binary"))
        .args(["--exact", "sheet_booleans_match_occ", "--nocapture"])
        .env(LOG, &log)
        .output()
        .expect("child run");
    eprintln!("{}", String::from_utf8_lossy(&child.stderr));
    assert!(
        child.status.success(),
        "corpus failed:\n{}",
        String::from_utf8_lossy(&child.stdout)
    );
    let fallbacks = std::fs::read_to_string(&log).unwrap_or_default();
    assert!(
        !fallbacks.contains("PartitionFallback"),
        "complex fell back:\n{fallbacks}"
    );
}

fn corpus() {
    let tol = Tolerance::default();
    let cube = json!({"kind": "box", "corners": [0, 0, 0, 1, 1, 1]});
    let ball = json!({"kind": "sphere", "centre": [0, 0, 0], "radius": 1.0});
    let corpus = [
        (
            "plane_through",
            &cube,
            json!({"kind": "plane", "x": [-1, 2], "y": [-1, 2], "z": 0.5}),
        ),
        (
            "plane_ending_inside",
            &cube,
            json!({"kind": "plane", "x": [-1, 0.5], "y": [-1, 2], "z": 0.5}),
        ),
        (
            "plane_inside",
            &cube,
            json!({"kind": "plane", "x": [0.25, 0.75], "y": [0.25, 0.75], "z": 0.5}),
        ),
        (
            "plane_outside",
            &cube,
            json!({"kind": "plane", "x": [2, 3], "y": [0, 1], "z": 0.5}),
        ),
        (
            "plane_oblique",
            &cube,
            json!({"kind": "plane", "x": [-1, 2], "y": [-1, 2], "z": 0.5, "tilt_x": 0.3}),
        ),
        (
            "tube_through",
            &cube,
            json!({"kind": "tube", "centre": [0.5, 0.5], "radius": 0.3, "z": [-0.5, 1.5]}),
        ),
        (
            "plane_through_poles",
            &ball,
            json!({"kind": "plane", "x": [-2, 2], "y": [-2, 2], "z": 0.0, "tilt_x": std::f64::consts::FRAC_PI_2}),
        ),
    ];
    let mut worst = (0.0_f64, String::new());
    let mut checked = 0;
    for (name, solid_spec, sheet_spec) in &corpus {
        for (k, place) in PLACEMENTS.iter().enumerate() {
            let row = format!("{name}@{k}");
            let t = place.transform();
            let body = solid(solid_spec, tol).transform(&t).expect("place solid");
            let s = sheet(sheet_spec, tol).transform(&t).expect("place sheet");
            let common = open_area(
                boolean_with_sheets(BoolOperation::Intersection, &body, &s, tol),
                &row,
                &tol,
            );
            let cut = open_area(
                boolean_with_sheets(BoolOperation::Subtraction, &s, &body, tol),
                &row,
                &tol,
            );
            let split = volumes(
                boolean_with_sheets(BoolOperation::Subtraction, &body, &s, tol),
                &row,
                &tol,
            );
            let spec = json!({"solid": solid_spec, "sheet": sheet_spec, "placement": place.json()});
            let Some(occ) = oracle(&spec) else { continue };
            assert_eq!(occ["done"], json!(true), "{row}: OCC {occ}");
            let scale2 = place.scale * place.scale;
            let occ_split: Vec<f64> = occ["split_volumes"]
                .as_array()
                .expect("split")
                .iter()
                .map(|v| v.as_f64().expect("v"))
                .collect();
            assert_eq!(
                split.len(),
                occ_split.len(),
                "{row}: pieces {split:?} vs OCC {occ_split:?}"
            );
            let mut gaps = vec![
                (
                    "common",
                    gap(
                        common,
                        occ["common_area"].as_f64().expect("a"),
                        1e-9 * scale2,
                    ),
                ),
                (
                    "cut",
                    gap(cut, occ["cut_area"].as_f64().expect("a"), 1e-9 * scale2),
                ),
            ];
            for (a, b) in split.iter().zip(&occ_split) {
                gaps.push(("split", gap(*a, *b, 1e-9 * scale2 * place.scale)));
            }
            for (what, g) in gaps {
                if g > worst.0 {
                    worst = (g, format!("{row} {what}"));
                }
            }
            checked += 1;
        }
    }
    eprintln!(
        "[sheet-oracle] {checked} rows; worst relative gap {:.3e} at {}",
        worst.0, worst.1
    );
    assert!(
        worst.0 < 1e-6,
        "worst relative gap {:.3e} at {}",
        worst.0,
        worst.1
    );
}
