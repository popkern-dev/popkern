//! PK-344 acceptance gates: double-seam closed-rim anchors (junction
//! circles between two seamed revolution faces) serve through the ARC
//! join machinery at the full geometric envelope — no Intersection
//! downgrade, validate_all-clean, volume-sane, OCC-clean where the
//! oracle is present.
//!
//! Mechanism (measured in): the corner plan retains the min-handle run
//! at a closed-rim anchor; seams ending at the anchor re-author after
//! the corner moves the vertices (their recompute declined against the
//! interim anchor position); and conic-NURBS circular spines (boolean
//! junction circles) promote to EXACT tori (the general canal fit lost
//! 4.5e-2 at the near-spindle bore-vanish grow).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::{bracket_bore_countersink, sphere_capped_rod};
use popkern_shell::offset::offset_solid;
use popkern_shell::{ShellDirection, ShellOptions, ShellWarning, max_safe_offset, shell};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn assert_no_downgrade(warnings: &[ShellWarning], label: &str) {
    assert!(
        !warnings.iter().any(|w| matches!(
            w,
            ShellWarning::JoinStyleDowngraded { .. }
                | ShellWarning::JoinStyleDowngradedAtEdge { .. }
        )),
        "{label}: arc joins must serve without downgrade: {warnings:?}"
    );
}

fn mesh_volume(brep: &BRep) -> f64 {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    let mut vol = 0.0_f64;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        vol += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (vol / 6.0).abs()
}

#[test]
fn bracket_outward_serves_arc_joins_to_the_bore_vanish_max() {
    let f = bracket_bore_countersink(tol()).expect("fixture");
    let source_volume = mesh_volume(&f.brep);
    let report = max_safe_offset(&f.brep, ShellDirection::Outward, tol()).expect("bound");
    // The geometric envelope: 0.99 × the bore radius (0.2).
    assert!(
        (report.max_distance - 0.198).abs() < 1e-9,
        "bore vanish binds: {report}"
    );
    for d in [0.09, 0.15, report.max_distance] {
        let r = offset_solid(&f.brep, d, &ShellOptions::default(), tol())
            .unwrap_or_else(|e| panic!("outward grow at {d} must serve: {e}"));
        assert!(
            r.brep.validate_all().is_empty(),
            "grow at {d}: validate_all-clean"
        );
        assert_no_downgrade(&r.warnings, "bracket outward");
        // Volume sanity: a grow strictly increases volume.
        let grown = mesh_volume(&r.brep);
        assert!(
            grown > source_volume,
            "grow at {d}: volume {grown} must exceed source {source_volume}"
        );
    }
}

#[test]
fn rod_inward_serves_arc_joins_to_the_wall_vanish_max() {
    let f = sphere_capped_rod(tol()).expect("fixture");
    let report = max_safe_offset(&f.brep, ShellDirection::Inward, tol()).expect("bound");
    assert!(
        (report.max_distance - 0.297).abs() < 1e-9,
        "cylinder vanish binds: {report}"
    );
    for d in [0.1, 0.2, report.max_distance] {
        let r = shell(
            &f.brep,
            d,
            ShellDirection::Inward,
            &[],
            &[],
            &ShellOptions::default(),
            tol(),
            None,
        )
        .unwrap_or_else(|e| panic!("inward hollow at {d} must serve: {e}"));
        assert!(
            r.brep.validate_all().is_empty(),
            "hollow at {d}: validate_all-clean"
        );
        assert_no_downgrade(&r.warnings, "rod inward");
    }
}

/// OCC oracle leg (skips without the OCP venv): the max-grow bracket
/// and max-hollow rod as STEP round-trips through OpenCASCADE.
#[test]
fn closed_rim_serves_pass_the_occ_oracle() {
    let oracle_python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let oracle_script = std::path::Path::new("../../tools/occ-oracle/occ_check.py");
    if !oracle_python.exists() || !oracle_script.exists() {
        eprintln!("[pk344-occ] OCP venv absent — leg skipped");
        return;
    }
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("pk344_occ");
    std::fs::create_dir_all(&dir).expect("dir");

    let mut bodies: Vec<(&str, BRep)> = Vec::new();
    let f = bracket_bore_countersink(tol()).expect("fixture");
    bodies.push((
        "bracket_grow_max",
        offset_solid(&f.brep, 0.198, &ShellOptions::default(), tol())
            .expect("serves")
            .brep,
    ));
    let g = sphere_capped_rod(tol()).expect("fixture");
    bodies.push((
        "rod_hollow_max",
        shell(
            &g.brep,
            0.297,
            ShellDirection::Inward,
            &[],
            &[],
            &ShellOptions::default(),
            tol(),
            None,
        )
        .expect("serves")
        .brep,
    ));
    for (name, brep) in bodies {
        let path = dir.join(format!("{name}.step"));
        let doc = popkern_io::StepDocument {
            products: vec![popkern_io::StepProductNode {
                name: name.into(),
                id: name.to_uppercase(),
                children: vec![],
                shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                    kind: popkern_io::ShapeRepKind::AdvancedBrep,
                    solids: vec![brep],
                    context: popkern_io::RepresentationContext::default(),
                    non_brep_items: vec![],
                })],
                metadata: popkern_io::ProductMetadata::default(),
            }],
            roots: vec![0],
            units: popkern_io::UnitContext::default(),
            opaque_bag: popkern_io::OpaqueBag::default(),
        };
        popkern_io::export_step_file(&doc, &path, &popkern_io::StepExportOptions::default())
            .expect("export");
        let output = std::process::Command::new(&oracle_python)
            .arg(oracle_script)
            .arg(&path)
            .output()
            .expect("oracle runs");
        let report = String::from_utf8_lossy(&output.stdout);
        assert!(
            report.contains("\"passed\": true"),
            "{name}: OCC oracle: {}",
            report.chars().take(400).collect::<String>()
        );
    }
}

/// Visual verification (geometry-affecting change rule): the bracket
/// grown to the bore-vanish max — the junction ring's exact torus
/// strip wraps the nearly-closed bore.
#[test]
fn render_pk344_bracket_max_grow() {
    use glam::Vec3 as GVec3;
    use popkern_viz::camera::Camera;
    let f = bracket_bore_countersink(tol()).expect("fixture");
    let grown = offset_solid(&f.brep, 0.198, &ShellOptions::default(), tol())
        .expect("serves")
        .brep;
    let out_dir = std::path::Path::new("../../examples/output");
    std::fs::create_dir_all(out_dir).expect("output dir");
    let camera = Camera {
        eye: GVec3::new(2.9, 2.5, 2.6),
        target: GVec3::new(1.0, 0.5, 0.2),
        up: GVec3::new(0.0, 0.0, 1.0),
        fov_y: 45.0_f32.to_radians(),
        aspect: 800.0 / 600.0,
        near: 0.01,
        far: 100.0,
    };
    let path = out_dir.join("pk344_bracket_max_grow.png");
    popkern_viz::headless::render_to_png(&grown, &camera, 800, 600, &path, true).expect("render");
}
