//! PK-323 row 468b: OCC oracle battery — every closed corpus fixture
//! and every currently-served shell/offset result exports to STEP and
//! passes OpenCASCADE's BRepCheck with zero open free bounds (the
//! interop bar, design §9 oracle #4).
//!
//! The leg SKIPS when the OCP venv is absent so CI without the oracle
//! stays green (established seam-matrix pattern).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::corpus;
use popkern_shell::offset::offset_solid;
use popkern_shell::{ShellDirection, ShellOptions, shell, thicken::thicken};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn export_step_solid(brep: &BRep, path: &std::path::Path) -> Result<(), String> {
    let name = path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "fixture".into());
    let doc = popkern_io::StepDocument {
        products: vec![popkern_io::StepProductNode {
            name: name.clone(),
            id: name.to_uppercase(),
            children: vec![],
            shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                kind: popkern_io::ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: popkern_io::RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: popkern_io::ProductMetadata::default(),
        }],
        roots: vec![0],
        units: popkern_io::UnitContext::default(),
        opaque_bag: popkern_io::OpaqueBag::default(),
    };
    popkern_io::export_step_file(&doc, path, &popkern_io::StepExportOptions::default())
        .map_err(|e| format!("{e:?}"))
}

/// Run the OCC oracle on an exported STEP file. `None` when the OCP
/// venv is absent (the leg skips; CI without the oracle stays green).
fn occ_oracle_passes(path: &std::path::Path) -> Option<bool> {
    let oracle_python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let oracle_script = std::path::Path::new("../../tools/occ-oracle/occ_check.py");
    if !oracle_python.exists() || !oracle_script.exists() {
        return None;
    }
    let output = std::process::Command::new(&oracle_python)
        .arg(oracle_script)
        .arg(path)
        .output()
        .ok()?;
    let report = String::from_utf8_lossy(&output.stdout);
    eprintln!(
        "[occ-battery] {}: {}",
        path.display(),
        report.chars().take(300).collect::<String>()
    );
    Some(report.contains("\"passed\": true"))
}

fn top_plane_face(brep: &BRep, z: f64) -> popkern_core::topology::FaceHandle {
    use popkern_core::geometry::surface::SurfaceGeometry;
    brep.faces
        .keys()
        .find(|&fh| {
            matches!(
                brep.faces[fh].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { .. })
            ) && {
                let ps = brep.face_outer_loop_positions(fh);
                !ps.is_empty() && ps.iter().all(|p| (p.z - z).abs() < 1e-9)
            }
        })
        .expect("plane face at z")
}

#[test]
fn corpus_and_served_shells_pass_the_occ_oracle() {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("pk323_occ_battery");
    std::fs::create_dir_all(&dir).expect("battery dir");

    let mut bodies: Vec<(String, BRep)> = Vec::new();
    // Closed corpus fixtures (the open sheet has free bounds by
    // design — it is exercised by the thicken-rims stage instead).
    for f in corpus(tol()).expect("corpus builds") {
        if f.watertight {
            bodies.push((f.name.to_string(), f.brep));
        }
    }
    // Served shell/offset results (the current serving envelope).
    let cube = BRep::unit_cube(tol()).unwrap();
    bodies.push((
        "thickened_cube".into(),
        thicken(&cube, 0.1, &ShellOptions::default(), tol())
            .expect("thicken cube")
            .brep,
    ));
    let top = top_plane_face(&cube, 1.0);
    bodies.push((
        "shelled_cube".into(),
        shell(
            &cube,
            0.1,
            ShellDirection::Inward,
            &[top],
            &[],
            &ShellOptions::default(),
            tol(),
            None,
        )
        .expect("shell cube")
        .brep,
    ));
    let cyl = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let cyl_top = top_plane_face(&cyl, 2.0);
    bodies.push((
        "shelled_cylinder".into(),
        shell(
            &cyl,
            0.1,
            ShellDirection::Inward,
            &[cyl_top],
            &[],
            &ShellOptions::default(),
            tol(),
            None,
        )
        .expect("shell cylinder")
        .brep,
    ));
    let cone = BRep::cone(1.0, 2.0, tol()).unwrap();
    bodies.push((
        "offset_cone".into(),
        offset_solid(&cone, 0.15, &ShellOptions::default(), tol())
            .expect("offset cone")
            .brep,
    ));
    // PK-324 serving envelope: curved recompute arms + open-sheet rims.
    let sphere = BRep::sphere(0.5, tol()).unwrap();
    bodies.push((
        "thickened_sphere".into(),
        thicken(&sphere, 0.08, &ShellOptions::default(), tol())
            .expect("thicken sphere")
            .brep,
    ));
    let torus = BRep::torus(0.6, 0.25, tol()).unwrap();
    bodies.push((
        "thickened_torus".into(),
        thicken(&torus, 0.08, &ShellOptions::default(), tol())
            .expect("thicken torus")
            .brep,
    ));
    let filleted = popkern_examples::offset_corpus::filleted_cube(tol())
        .expect("filleted cube")
        .brep;
    bodies.push((
        "thickened_filleted_cube".into(),
        thicken(&filleted, 0.05, &ShellOptions::default(), tol())
            .expect("thicken filleted cube")
            .brep,
    ));
    // PK-327 global-SI trims: the welded groove crease and the
    // pinch-pole apple are the riskiest interop classes the trim
    // stage produces.
    let groove = popkern_examples::offset_corpus::grooved_cylinder(tol())
        .expect("grooved cylinder")
        .brep;
    bodies.push((
        "trimmed_grooved_cylinder_grow".into(),
        offset_solid(&groove, 0.08, &ShellOptions::default(), tol())
            .expect("trimmed groove grow")
            .brep,
    ));
    let lens = popkern_examples::offset_corpus::torus_lens(tol())
        .expect("torus lens")
        .brep;
    bodies.push((
        "torus_lens_apple".into(),
        offset_solid(&lens, 0.2, &ShellOptions::default(), tol())
            .expect("apple grow")
            .brep,
    ));
    let open = BRep::open_box(tol()).unwrap();
    bodies.push((
        "thickened_open_box".into(),
        thicken(&open, 0.06, &ShellOptions::default(), tol())
            .expect("thicken open box")
            .brep,
    ));
    // PK-325 arc joins: the rounded box (strips + sphere octants) and
    // the wavy-plate acceptance body under the Minkowski default.
    let cube = BRep::unit_cube(tol()).unwrap();
    bodies.push((
        "rounded_cube_offset".into(),
        offset_solid(&cube, 0.15, &ShellOptions::default(), tol())
            .expect("rounded cube offset")
            .brep,
    ));
    if let Some(wavy) = popkern_examples::offset_corpus::corpus(tol())
        .ok()
        .and_then(|fs| fs.into_iter().find(|f| f.name == "wavy_plate"))
    {
        bodies.push((
            "thickened_wavy_plate_arc".into(),
            thicken(&wavy.brep, 0.04, &ShellOptions::default(), tol())
                .expect("thicken wavy plate")
                .brep,
        ));
    }

    let mut skipped = false;
    let mut failures = Vec::new();
    for (name, brep) in &bodies {
        let path = dir.join(format!("{name}.step"));
        export_step_solid(brep, &path).unwrap_or_else(|e| panic!("{name}: STEP export: {e}"));
        match occ_oracle_passes(&path) {
            None => {
                skipped = true;
                break;
            }
            Some(true) => {}
            Some(false) => failures.push(name.clone()),
        }
    }
    if skipped {
        eprintln!("[occ-battery] OCP venv absent — leg skipped");
        return;
    }
    assert!(
        failures.is_empty(),
        "OCC oracle failures (BRepCheck or open free bounds): {failures:?}"
    );
}
