//! STEP import of non-manifold bodies (PK-980).
//!
//! Each fixture is one `MANIFOLD_SOLID_BREP` whose `CLOSED_SHELL` lists
//! the faces of several solids that touch, written by OpenCASCADE
//! (`tools/occ-oracle/occ_write_nonmanifold.py`): a contact edge with
//! four face uses, or a vertex two closed sets of faces share. By
//! default the import serves its manifold pieces as separate bodies
//! (PK-1115), gated on the solids OpenCASCADE reads back from the same
//! file (per-piece volume multiset, 1e-6 relative), each body's clean
//! `validate_all`, the bodies touching where the shell did, and every
//! body meshing watertight.

use popkern_bool::{BoolOptions, NonManifoldPolicy, SetExpr, nary_boolean};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_io::step::diagnostics::Severity;
use popkern_io::{StepImportOptions, import_step};

fn fixture(name: &str) -> Vec<u8> {
    let path = format!(
        "{}/../popkern-io/tests/fixtures/nonmanifold/{name}.step",
        env!("CARGO_MANIFEST_DIR")
    );
    std::fs::read(&path).unwrap_or_else(|e| panic!("{path}: {e}"))
}

/// The bodies the file imports, and the import's errors, under the
/// default import options (healing on: its vertex sewing must leave the
/// pieces' contact copies apart).
fn import(name: &str) -> (Vec<BRep>, Vec<String>) {
    let imported = import_step(&fixture(name), &StepImportOptions::default()).expect("import");
    let errors = imported
        .diagnostics
        .entries()
        .iter()
        .filter(|e| e.severity == Severity::Error)
        .map(|e| format!("{e:?}"))
        .collect();
    let bodies: Vec<BRep> = imported
        .document
        .products
        .iter()
        .flat_map(|p| &p.shapes)
        .flat_map(|s| s.solids.iter().cloned())
        .collect();
    (bodies, errors)
}

/// Whether some vertex of `a` lies on a vertex of `b`.
fn touch(a: &BRep, b: &BRep) -> bool {
    a.vertices.values().any(|v| {
        b.vertices
            .values()
            .any(|w| v.position.distance_to(&w.position) <= a.tolerance.linear())
    })
}

struct Row {
    name: &'static str,
    /// The solids OpenCASCADE reads back from the file, ascending.
    occ_volumes: &'static [f64],
}

const ROWS: &[Row] = &[
    Row {
        name: "nm_cubes_on_an_edge",
        occ_volumes: &[1.0, 1.0],
    },
    Row {
        name: "nm_cubes_at_a_corner",
        occ_volumes: &[1.0, 1.0],
    },
    Row {
        name: "nm_three_blocks_on_an_edge",
        occ_volumes: &[500.0, 500.0, 500.0],
    },
    Row {
        name: "nm_cylinder_on_a_block",
        occ_volumes: &[18.84955592153876, 72.0],
    },
];

#[test]
fn non_manifold_shells_import_as_separate_touching_bodies() {
    for row in ROWS {
        let (bodies, errors) = import(row.name);
        assert!(errors.is_empty(), "{}: {errors:?}", row.name);
        assert_eq!(bodies.len(), row.occ_volumes.len(), "{}", row.name);
        for body in &bodies {
            let findings = body.validate_all();
            assert!(findings.is_empty(), "{}: {findings:?}", row.name);
            assert!(
                bodies
                    .iter()
                    .any(|other| !std::ptr::eq(other, body) && touch(body, other)),
                "{}: each body touches another",
                row.name
            );
        }
        let mut volumes: Vec<f64> = bodies
            .iter()
            .map(|p| popkern_measure::body_volume(p, &p.tolerance).expect("volume"))
            .collect();
        volumes.sort_by(f64::total_cmp);
        for (ours, occ) in volumes.iter().zip(row.occ_volumes) {
            let gap = (ours - occ).abs() / occ;
            assert!(
                gap < 1e-6,
                "{}: volume {ours} vs OCC {occ} ({gap:.2e})",
                row.name
            );
        }
        for body in &bodies {
            let mesh = popkern_mesh::tessellate(body, &Default::default()).expect("mesh");
            assert_eq!(popkern_mesh::welded_open_edges(&mesh), 0, "{}", row.name);
        }
    }
}

/// Imported touching bodies are operands like any: a block bridging the
/// two cubes on an edge unions them into one solid region, which still
/// touches itself where the cubes meet above and below the block (a
/// general body; the default policy declines it).
#[test]
fn a_boolean_bridges_an_imported_contact() {
    let (bodies, _) = import("nm_cubes_on_an_edge");
    let tol = Tolerance::default();
    let square = [
        Point3::new(0.5, 0.5, 0.2),
        Point3::new(1.5, 0.5, 0.2),
        Point3::new(1.5, 1.5, 0.2),
        Point3::new(0.5, 1.5, 0.2),
    ];
    let bridge = BRep::prism(&square, Vec3::new(0.0, 0.0, 1.0), 0.6, tol).expect("bridge");
    let operands = [&bodies[0], &bodies[1], &bridge];
    assert!(nary_boolean(&operands, &SetExpr::Any, tol, &BoolOptions::default()).is_err());
    let union = nary_boolean(
        &operands,
        &SetExpr::Any,
        tol,
        &BoolOptions {
            non_manifold: NonManifoldPolicy::General,
        },
    )
    .expect("union");
    assert_eq!(union.solid_regions().len(), 1);
    let volume = popkern_measure::body_volume(&union, &tol).expect("volume");
    assert!((volume - 2.3).abs() < 1e-9, "{volume}");
    assert!(union.validate_all().is_empty());
}
