//! STEP import of non-manifold bodies (PK-980).
//!
//! Each fixture is one `MANIFOLD_SOLID_BREP` whose `CLOSED_SHELL` lists
//! the faces of several solids that touch, written by OpenCASCADE
//! (`tools/occ-oracle/occ_write_nonmanifold.py`): a contact edge with
//! four face uses, or a vertex two closed sets of faces share. The
//! import serves manifold pieces meeting at recorded contacts, gated on
//! the solids OpenCASCADE reads back from the same file (per-piece
//! volume multiset, 1e-6 relative), a clean `validate_all`, the
//! contacts recorded, and every piece meshing watertight.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, TopologyContactKind};
use popkern_io::step::diagnostics::Severity;
use popkern_io::{StepImportOptions, import_step};

fn fixture(name: &str) -> Vec<u8> {
    let path = format!(
        "{}/../popkern-io/tests/fixtures/nonmanifold/{name}.step",
        env!("CARGO_MANIFEST_DIR")
    );
    std::fs::read(&path).unwrap_or_else(|e| panic!("{path}: {e}"))
}

/// The one body the file imports, and the import's errors, under the
/// default import options (healing on: its vertex sewing must leave the
/// pieces' contact copies apart).
fn import(name: &str) -> (BRep, Vec<String>) {
    let imported = import_step(&fixture(name), &StepImportOptions::default()).expect("import");
    let errors = imported
        .diagnostics
        .entries()
        .iter()
        .filter(|e| e.severity == Severity::Error)
        .map(|e| format!("{e:?}"))
        .collect();
    let mut bodies: Vec<BRep> = imported
        .document
        .products
        .iter()
        .flat_map(|p| &p.shapes)
        .flat_map(|s| s.solids.iter().cloned())
        .collect();
    assert_eq!(bodies.len(), 1, "{name}: one body ({errors:?})");
    (bodies.remove(0), errors)
}

struct Row {
    name: &'static str,
    /// The solids OpenCASCADE reads back from the file, ascending.
    occ_volumes: &'static [f64],
    edge_contacts: usize,
    vertex_contacts: usize,
}

const ROWS: &[Row] = &[
    Row {
        name: "nm_cubes_on_an_edge",
        occ_volumes: &[1.0, 1.0],
        edge_contacts: 1,
        vertex_contacts: 2,
    },
    Row {
        name: "nm_cubes_at_a_corner",
        occ_volumes: &[1.0, 1.0],
        edge_contacts: 0,
        vertex_contacts: 1,
    },
    Row {
        name: "nm_three_blocks_on_an_edge",
        occ_volumes: &[500.0, 500.0, 500.0],
        edge_contacts: 3,
        vertex_contacts: 4,
    },
    Row {
        name: "nm_cylinder_on_a_block",
        occ_volumes: &[18.84955592153876, 72.0],
        edge_contacts: 1,
        vertex_contacts: 2,
    },
];

#[test]
fn non_manifold_shells_import_as_pieces_meeting_at_recorded_contacts() {
    for row in ROWS {
        let (body, errors) = import(row.name);
        assert!(errors.is_empty(), "{}: {errors:?}", row.name);
        let findings = body.validate_all();
        assert!(findings.is_empty(), "{}: {findings:?}", row.name);
        let count = |kind| {
            body.topology_contacts
                .iter()
                .filter(|c| c.kind == kind)
                .count()
        };
        assert_eq!(
            (
                count(TopologyContactKind::Edge),
                count(TopologyContactKind::Vertex)
            ),
            (row.edge_contacts, row.vertex_contacts),
            "{}: contacts {:?}",
            row.name,
            body.topology_contacts
        );
        let pieces = popkern_bool::pieces::piece_bodies(&body).expect("pieces");
        let mut volumes: Vec<f64> = pieces
            .iter()
            .map(|p| popkern_measure::body_volume(p, &p.tolerance).expect("volume"))
            .collect();
        volumes.sort_by(f64::total_cmp);
        assert_eq!(volumes.len(), row.occ_volumes.len(), "{}", row.name);
        for (ours, occ) in volumes.iter().zip(row.occ_volumes) {
            let gap = (ours - occ).abs() / occ;
            assert!(
                gap < 1e-6,
                "{}: volume {ours} vs OCC {occ} ({gap:.2e})",
                row.name
            );
        }
        for piece in &pieces {
            let mesh = popkern_mesh::tessellate(piece, &Default::default()).expect("mesh");
            assert_eq!(popkern_mesh::welded_open_edges(&mesh), 0, "{}", row.name);
        }
    }
}

/// An imported contact is a contact a later boolean carries: a block
/// bridging the two cubes on an edge unions them into one solid.
#[test]
fn a_boolean_bridges_an_imported_contact() {
    let (body, _) = import("nm_cubes_on_an_edge");
    let tol = Tolerance::default();
    let square = [
        Point3::new(0.5, 0.5, 0.2),
        Point3::new(1.5, 0.5, 0.2),
        Point3::new(1.5, 1.5, 0.2),
        Point3::new(0.5, 1.5, 0.2),
    ];
    let bridge = BRep::prism(&square, Vec3::new(0.0, 0.0, 1.0), 0.6, tol).expect("bridge");
    let union = boolean(BoolOperation::Union, body, bridge, tol, None).expect("union");
    assert_eq!(union.solids.len(), 1);
    let volume = popkern_measure::body_volume(&union, &tol).expect("volume");
    assert!((volume - 2.3).abs() < 1e-9, "{volume}");
    assert!(union.validate_all().is_empty());
}
