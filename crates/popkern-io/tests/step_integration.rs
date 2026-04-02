//! End-to-end integration tests for the STEP import/export pipeline.
//!
//! These tests exercise the full round trip: parse STEP text, build the entity
//! graph, map topology + product structure, heal, export, and re-import.

use popkern_io::{StepExportOptions, StepImportOptions, export_step, import_step};

// ---------------------------------------------------------------------------
// Test 1: Import the hand-written cube fixture
// ---------------------------------------------------------------------------

#[test]
fn import_cube_step_file() {
    let data = include_bytes!("fixtures/cube.step");
    let opts = StepImportOptions {
        healing: popkern_io::step::healing::HealingLevel::None,
        ..Default::default()
    };
    let result = import_step(data, &opts).unwrap();

    // Should have at least 1 product
    assert!(
        !result.document.products.is_empty(),
        "expected at least 1 product"
    );

    // Product name should match the fixture
    assert_eq!(result.document.products[0].name, "UnitCube");

    // Should have at least 1 solid
    let total_solids: usize = result
        .document
        .products
        .iter()
        .flat_map(|p| &p.shapes)
        .map(|s| s.solids.len())
        .sum();
    assert!(
        total_solids >= 1,
        "expected at least 1 solid, got {total_solids}"
    );

    // The cube solid should have 8 vertices, 12 edges, 6 faces
    let brep = &result.document.products[0].shapes[0].solids[0];
    assert_eq!(brep.vertices.len(), 8, "cube must have 8 vertices");
    assert_eq!(brep.edges.len(), 12, "cube must have 12 edges");
    assert_eq!(brep.faces.len(), 6, "cube must have 6 faces");
    assert_eq!(brep.shells.len(), 1, "cube must have 1 shell");
    assert_eq!(brep.solids.len(), 1, "cube must have 1 solid");
}

// ---------------------------------------------------------------------------
// Test 2: Import the empty fixture
// ---------------------------------------------------------------------------

#[test]
fn import_empty_step_file() {
    let data = include_bytes!("fixtures/empty.step");
    let result = import_step(data, &StepImportOptions::default()).unwrap();
    assert!(
        result.document.products.is_empty(),
        "empty file should produce no products"
    );
}

// ---------------------------------------------------------------------------
// Test 3: Import file with broken references — lenient parsing
// ---------------------------------------------------------------------------

#[test]
fn import_broken_refs_is_lenient() {
    let data = include_bytes!("fixtures/broken_refs.step");
    let result = import_step(data, &StepImportOptions::default()).unwrap();
    // Should have diagnostics about the broken ref (#999 does not exist)
    assert!(
        result.diagnostics.has_errors(),
        "broken refs should produce diagnostic errors"
    );
    // Import itself should not fail — errors are reported, not fatal
}

// ---------------------------------------------------------------------------
// Test 4: Round-trip — import cube → export → re-import
// ---------------------------------------------------------------------------

#[test]
fn round_trip_cube() {
    let data = include_bytes!("fixtures/cube.step");
    let import_opts = StepImportOptions {
        healing: popkern_io::step::healing::HealingLevel::None,
        ..Default::default()
    };
    let result1 = import_step(data, &import_opts).unwrap();

    // Export
    let exported = export_step(&result1.document, &StepExportOptions::default()).unwrap();

    // Verify exported data is valid Part 21
    let text = String::from_utf8_lossy(&exported);
    assert!(
        text.starts_with("ISO-10303-21;"),
        "exported data must start with Part 21 header"
    );
    assert!(
        text.contains("END-ISO-10303-21;"),
        "exported data must contain Part 21 terminator"
    );

    // Re-import
    let result2 = import_step(&exported, &import_opts).unwrap();

    // Compare topology counts if both have solids
    if let (Some(s1), Some(s2)) = (first_solid(&result1), first_solid(&result2)) {
        assert_eq!(
            s1.vertices.len(),
            s2.vertices.len(),
            "vertex count mismatch after round-trip"
        );
        assert_eq!(
            s1.edges.len(),
            s2.edges.len(),
            "edge count mismatch after round-trip"
        );
        assert_eq!(
            s1.faces.len(),
            s2.faces.len(),
            "face count mismatch after round-trip"
        );
    }
}

/// Helper: extract the first BRep solid from a STEP import result.
fn first_solid(
    result: &popkern_io::StepImportResult,
) -> Option<&popkern_core::topology::brep::BRep> {
    result
        .document
        .products
        .first()
        .and_then(|p| p.shapes.first())
        .and_then(|s| s.solids.first())
}

// ---------------------------------------------------------------------------
// Test 5: Export a programmatically-built cube and verify STEP content
// ---------------------------------------------------------------------------

#[test]
fn export_programmatic_cube() {
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;
    use popkern_io::{
        OpaqueBag, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument, StepProductNode,
        UnitContext,
    };
    use std::sync::Arc;

    let brep = BRep::unit_cube(Tolerance::default()).unwrap();

    let doc = StepDocument {
        products: vec![StepProductNode {
            name: "TestCube".into(),
            id: "CUBE-001".into(),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: Default::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };

    let exported = export_step(&doc, &StepExportOptions::default()).unwrap();

    // Should be valid Part 21
    let text = String::from_utf8_lossy(&exported);
    assert!(text.contains("ISO-10303-21;"));
    assert!(text.contains("CARTESIAN_POINT"));
    assert!(text.contains("MANIFOLD_SOLID_BREP"));
    assert!(text.contains("PRODUCT("));
    assert!(text.contains("END-ISO-10303-21;"));

    // Should be re-importable with product structure intact
    let import_opts = StepImportOptions {
        healing: popkern_io::step::healing::HealingLevel::None,
        ..Default::default()
    };
    let result =
        import_step(&exported, &import_opts).expect("re-import of exported STEP must not fail");
    assert!(
        !result.document.products.is_empty(),
        "re-imported document should have products"
    );
    assert!(
        !result.document.products[0].shapes.is_empty(),
        "re-imported product should have shapes"
    );
}

// ---------------------------------------------------------------------------
// Test 6: Schema detection — the cube fixture declares AUTOMOTIVE_DESIGN
// ---------------------------------------------------------------------------

#[test]
fn detects_schema_from_header() {
    let data = include_bytes!("fixtures/cube.step");
    let result = import_step(data, &StepImportOptions::default()).unwrap();
    // Import should succeed for an AUTOMOTIVE_DESIGN schema file
    assert!(
        !result.document.products.is_empty(),
        "AUTOMOTIVE_DESIGN schema should be accepted"
    );
}

// ---------------------------------------------------------------------------
// Test 7: Product filter — only import matching products
// ---------------------------------------------------------------------------

#[test]
fn product_filter_excludes_unmatched() {
    let data = include_bytes!("fixtures/cube.step");
    let opts = StepImportOptions {
        product_filter: Some(vec!["NonExistentProduct".into()]),
        healing: popkern_io::step::healing::HealingLevel::None,
        ..Default::default()
    };
    let result = import_step(data, &opts).unwrap();

    // Product structure is still loaded but shapes are cleared for non-matching names
    let total_shapes: usize = result
        .document
        .products
        .iter()
        .flat_map(|p| &p.shapes)
        .map(|s| s.solids.len())
        .sum();
    assert_eq!(
        total_shapes, 0,
        "filtered import should have 0 solids for non-matching filter"
    );
}

// ---------------------------------------------------------------------------
// Test 8: Export of empty document produces valid Part 21
// ---------------------------------------------------------------------------

#[test]
fn export_empty_document() {
    use popkern_io::{OpaqueBag, StepDocument, UnitContext};

    let doc = StepDocument {
        products: vec![],
        roots: vec![],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };

    let exported = export_step(&doc, &StepExportOptions::default()).unwrap();
    let text = String::from_utf8_lossy(&exported);
    assert!(text.contains("ISO-10303-21;"));
    assert!(text.contains("END-ISO-10303-21;"));

    // Re-import of empty export should succeed
    let result = import_step(&exported, &StepImportOptions::default()).unwrap();
    assert!(result.document.products.is_empty());
}

// ---------------------------------------------------------------------------
// Test 9: Programmatic cube → export → import → verify topology and positions
// ---------------------------------------------------------------------------

/// Build a unit cube from scratch, export to STEP, re-import, and verify that
/// the topology counts and vertex positions survive the round-trip.
#[test]
fn round_trip_programmatic_cube_topology_and_positions() {
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;
    use popkern_io::{
        OpaqueBag, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument, StepProductNode,
        UnitContext,
    };
    use std::sync::Arc;

    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).unwrap();

    // Collect vertex positions from the original BRep.
    let mut original_positions: Vec<[f64; 3]> = brep
        .vertices
        .values()
        .map(|v| [v.position.x, v.position.y, v.position.z])
        .collect();
    original_positions.sort_by(|a, b| {
        a[0].partial_cmp(&b[0])
            .unwrap()
            .then(a[1].partial_cmp(&b[1]).unwrap())
            .then(a[2].partial_cmp(&b[2]).unwrap())
    });

    let original_face_count = brep.faces.len();
    let original_edge_count = brep.edges.len();
    let original_vertex_count = brep.vertices.len();

    // Build a StepDocument wrapping the BRep.
    let doc = StepDocument {
        products: vec![StepProductNode {
            name: "RoundTripCube".into(),
            id: "RT-001".into(),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: Default::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };

    // Export to STEP bytes.
    let exported = export_step(&doc, &StepExportOptions::default()).expect("export should succeed");

    // Re-import with no healing to get raw topology.
    let import_opts = StepImportOptions {
        healing: popkern_io::step::healing::HealingLevel::None,
        ..Default::default()
    };
    let result = import_step(&exported, &import_opts).expect("re-import should succeed");

    // Verify product structure.
    assert!(
        !result.document.products.is_empty(),
        "re-imported document should have products"
    );
    assert_eq!(result.document.products[0].name, "RoundTripCube");

    // Extract the first solid.
    let reimported_brep = first_solid(&result).expect("should have at least one solid");

    // Verify topology counts.
    assert_eq!(
        reimported_brep.faces.len(),
        original_face_count,
        "face count mismatch: expected {original_face_count}, got {}",
        reimported_brep.faces.len()
    );
    assert_eq!(
        reimported_brep.edges.len(),
        original_edge_count,
        "edge count mismatch: expected {original_edge_count}, got {}",
        reimported_brep.edges.len()
    );
    assert_eq!(
        reimported_brep.vertices.len(),
        original_vertex_count,
        "vertex count mismatch: expected {original_vertex_count}, got {}",
        reimported_brep.vertices.len()
    );

    // Verify vertex positions are within tolerance.
    let mut reimported_positions: Vec<[f64; 3]> = reimported_brep
        .vertices
        .values()
        .map(|v| [v.position.x, v.position.y, v.position.z])
        .collect();
    reimported_positions.sort_by(|a, b| {
        a[0].partial_cmp(&b[0])
            .unwrap()
            .then(a[1].partial_cmp(&b[1]).unwrap())
            .then(a[2].partial_cmp(&b[2]).unwrap())
    });

    let position_tol = 1e-6;
    assert_eq!(
        original_positions.len(),
        reimported_positions.len(),
        "sorted position list lengths must match"
    );
    for (i, (orig, reimp)) in original_positions
        .iter()
        .zip(reimported_positions.iter())
        .enumerate()
    {
        let dist = ((orig[0] - reimp[0]).powi(2)
            + (orig[1] - reimp[1]).powi(2)
            + (orig[2] - reimp[2]).powi(2))
        .sqrt();
        assert!(
            dist < position_tol,
            "vertex {i} position mismatch: original={orig:?}, reimported={reimp:?}, dist={dist}"
        );
    }
}

// ---------------------------------------------------------------------------
// Test 10: Double round-trip — export → import → export → import
// ---------------------------------------------------------------------------

/// Verifies that a second round-trip through STEP preserves topology.
#[test]
fn double_round_trip_preserves_topology() {
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;
    use popkern_io::{
        OpaqueBag, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument, StepProductNode,
        UnitContext,
    };
    use std::sync::Arc;

    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let original_vertices = brep.vertices.len();
    let original_edges = brep.edges.len();
    let original_faces = brep.faces.len();

    let doc = StepDocument {
        products: vec![StepProductNode {
            name: "DoubleCube".into(),
            id: "D-001".into(),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: Default::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };

    let import_opts = StepImportOptions {
        healing: popkern_io::step::healing::HealingLevel::None,
        ..Default::default()
    };
    let export_opts = StepExportOptions::default();

    // First round-trip.
    let step1 = export_step(&doc, &export_opts).expect("export 1");
    let result1 = import_step(&step1, &import_opts).expect("import 1");

    // Second round-trip.
    let step2 = export_step(&result1.document, &export_opts).expect("export 2");
    let result2 = import_step(&step2, &import_opts).expect("import 2");

    let solid = first_solid(&result2).expect("solid after double round-trip");
    assert_eq!(
        solid.vertices.len(),
        original_vertices,
        "vertex count after double round-trip"
    );
    assert_eq!(
        solid.edges.len(),
        original_edges,
        "edge count after double round-trip"
    );
    assert_eq!(
        solid.faces.len(),
        original_faces,
        "face count after double round-trip"
    );
}
