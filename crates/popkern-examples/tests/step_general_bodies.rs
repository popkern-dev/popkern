//! STEP import and export of general bodies (PK-1108, S9).
//!
//! The OpenCASCADE-written non-manifold fixtures (one `CLOSED_SHELL`
//! listing the faces of solids that touch along an edge or at a corner)
//! import under `StepNonManifold::General` as one general body: its
//! contact edge's ring holds every solid's faces, nothing is left
//! recorded, `validate_all` is clean, and each solid region's volume is
//! the solid OpenCASCADE reads back. Exported again it is the same
//! form: it reimports as the same general body, as the same pieces under
//! the default option, and OpenCASCADE reads the same solids from it.

use popkern_core::topology::{BRep, BodyKind};
use popkern_io::step::diagnostics::Severity;
use popkern_io::{StepExportOptions, StepImportOptions, StepNonManifold, export_step, import_step};

fn fixture(name: &str) -> Vec<u8> {
    let path = format!(
        "{}/../popkern-io/tests/fixtures/nonmanifold/{name}.step",
        env!("CARGO_MANIFEST_DIR")
    );
    std::fs::read(&path).unwrap_or_else(|e| panic!("{path}: {e}"))
}

/// The document and its one body under `non_manifold`.
fn import(
    data: &[u8],
    non_manifold: StepNonManifold,
    name: &str,
) -> (popkern_io::StepDocument, BRep) {
    let imported = import_step(
        data,
        &StepImportOptions {
            non_manifold,
            ..Default::default()
        },
    )
    .expect("import");
    let errors: Vec<_> = imported
        .diagnostics
        .entries()
        .iter()
        .filter(|e| e.severity == Severity::Error)
        .collect();
    assert!(errors.is_empty(), "{name}: {errors:?}");
    let mut bodies: Vec<BRep> = imported
        .document
        .products
        .iter()
        .flat_map(|p| &p.shapes)
        .flat_map(|s| s.solids.iter().cloned())
        .collect();
    assert_eq!(bodies.len(), 1, "{name}: one body");
    (imported.document, bodies.remove(0))
}

fn region_volumes(body: &BRep) -> Vec<f64> {
    let mut v: Vec<f64> = popkern_measure::region_volumes(body, &body.tolerance)
        .expect("volumes")
        .into_iter()
        .map(|(_, v)| v)
        .collect();
    v.sort_by(f64::total_cmp);
    v
}

fn assert_volumes(name: &str, ours: &[f64], want: &[f64]) {
    assert_eq!(ours.len(), want.len(), "{name}: {ours:?} vs {want:?}");
    for (o, w) in ours.iter().zip(want) {
        let gap = (o - w).abs() / w;
        assert!(gap < 1e-6, "{name}: {o} vs {w} ({gap:.2e})");
    }
}

/// The solids OpenCASCADE reads from `data`, ascending; `None` without
/// the OCP venv.
fn occ_volumes(data: &[u8], name: &str) -> Option<Vec<f64>> {
    let home = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default());
    let python = home.join("venvs/ocp/bin/python");
    let oracle = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tools/occ-oracle/occ_solids.py");
    if !python.exists() || !oracle.exists() {
        return None;
    }
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("step_general");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join(format!("{name}.step"));
    std::fs::write(&path, data).expect("write");
    let out = std::process::Command::new(&python)
        .arg(&oracle)
        .arg(&path)
        .output()
        .expect("oracle");
    let report: serde_json::Value = serde_json::from_slice(&out.stdout).expect("json");
    let mut v: Vec<f64> = report["solids"]
        .as_array()
        .expect("solids")
        .iter()
        .map(|s| s["volume"].as_f64().expect("volume"))
        .collect();
    v.sort_by(f64::total_cmp);
    Some(v)
}

const ROWS: &[(&str, &[f64])] = &[
    ("nm_cubes_on_an_edge", &[1.0, 1.0]),
    ("nm_cubes_at_a_corner", &[1.0, 1.0]),
    ("nm_three_blocks_on_an_edge", &[500.0, 500.0, 500.0]),
    ("nm_cylinder_on_a_block", &[18.84955592153876, 72.0]),
];

#[test]
fn self_touching_shells_import_as_general_bodies_and_export_again() {
    for &(name, want) in ROWS {
        let (doc, body) = import(&fixture(name), StepNonManifold::General, name);
        assert_eq!(body.body_kind(), BodyKind::General, "{name}");
        assert!(
            body.topology_contacts.is_empty(),
            "{name}: {:?}",
            body.topology_contacts
        );
        let findings = body.validate_all();
        assert!(findings.is_empty(), "{name}: {findings:?}");
        assert_volumes(name, &region_volumes(&body), want);

        let written = export_step(&doc, &StepExportOptions::default()).expect("export");
        let (_, again) = import(&written, StepNonManifold::General, name);
        assert_eq!(again.body_kind(), BodyKind::General, "{name}: re-read");
        let counts = |b: &BRep| (b.vertices.len(), b.edges.len(), b.faces.len());
        assert_eq!(counts(&again), counts(&body), "{name}: re-read counts");
        assert_volumes(name, &region_volumes(&again), want);

        let (_, pieces) = import(&written, StepNonManifold::Pieces, name);
        assert!(
            !pieces.topology_contacts.is_empty(),
            "{name}: pieces record contacts"
        );
        assert_eq!(pieces.solid_regions().len(), want.len(), "{name}: pieces");

        if let Some(occ) = occ_volumes(&written, name) {
            assert_volumes(name, &occ, want);
        }
    }
}

/// Shells no split into manifold pieces fits (PK-1070, PK-1071) import
/// face by face and glue into one general body whose enclosed cells are
/// its solids: a solid touching itself along an edge (one region, its
/// contact edge a ring of four), and two boxes sharing a face listed once
/// (two regions, the face between them). Both export (the internal face
/// listed once) and re-read; OpenCASCADE reads the same solids from
/// every file.
#[test]
fn shells_without_a_manifold_split_import_as_general_bodies() {
    for &(name, want, internal_faces) in &[
        ("nm_one_piece_touching_itself", &[4.0][..], 0usize),
        ("nm_boxes_sharing_a_face", &[1.0, 1.0][..], 1),
    ] {
        let data = fixture(name);
        let (doc, body) = import(&data, StepNonManifold::General, name);
        assert_eq!(body.body_kind(), BodyKind::General, "{name}");
        let findings = body.validate_all();
        assert!(findings.is_empty(), "{name}: {findings:?}");
        assert_volumes(name, &region_volumes(&body), want);
        let solid = |s| {
            body.shell_region(s).and_then(|r| body.region_kind(r))
                == Some(popkern_core::topology::RegionKind::Solid)
        };
        let between = body
            .faces
            .values()
            .filter(|f| solid(f.back) && solid(f.front))
            .count();
        assert_eq!(between, internal_faces, "{name}: internal faces");
        if let Some(occ) = occ_volumes(&data, name) {
            assert_volumes(name, &occ, want);
        }
        {
            let written = export_step(&doc, &StepExportOptions::default()).expect("export");
            let (_, again) = import(&written, StepNonManifold::General, name);
            assert_eq!(again.body_kind(), BodyKind::General, "{name}: re-read");
            assert_volumes(name, &region_volumes(&again), want);
            if let Some(occ) = occ_volumes(&written, &format!("{name}_written")) {
                assert_volumes(name, &occ, want);
            }
        }
    }
}
