use super::*;
use popkern_core::tolerance::Tolerance;
use std::sync::Arc;

use crate::step::document::{
    OpaqueBag, RepresentationContext, ShapeRep, ShapeRepKind, StepProductNode, UnitContext,
};
use crate::step::export::{StepExportOptions, export_step};

#[test]
fn export_unit_cube_entity_counts() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut ctx = ExportContext::new();
    let solid_id = ctx.export_brep(&brep).unwrap();

    // Count entity types
    let count =
        |name: &str| -> usize { ctx.entities.iter().filter(|e| e.type_name == name).count() };

    // A cube has 8 vertices, 12 edges, 6 faces
    assert_eq!(count("VERTEX_POINT"), 8, "Expected 8 VERTEX_POINTs");
    assert_eq!(count("EDGE_CURVE"), 12, "Expected 12 EDGE_CURVEs");
    assert_eq!(count("ADVANCED_FACE"), 6, "Expected 6 ADVANCED_FACEs");
    assert_eq!(count("CLOSED_SHELL"), 1, "Expected 1 CLOSED_SHELL");
    assert_eq!(
        count("MANIFOLD_SOLID_BREP"),
        1,
        "Expected 1 MANIFOLD_SOLID_BREP"
    );

    // Each face has one outer loop
    assert_eq!(count("FACE_OUTER_BOUND"), 6, "Expected 6 FACE_OUTER_BOUNDs");
    assert_eq!(count("EDGE_LOOP"), 6, "Expected 6 EDGE_LOOPs");

    // Each edge appears in two faces -> 24 oriented edges
    assert_eq!(count("ORIENTED_EDGE"), 24, "Expected 24 ORIENTED_EDGEs");

    // The solid ID should refer to a real entity
    let solid = ctx.entities.iter().find(|e| e.id == solid_id);
    assert!(solid.is_some());
    assert_eq!(solid.unwrap().type_name, "MANIFOLD_SOLID_BREP");
}

#[test]
fn export_cube_surfaces_are_planes() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut ctx = ExportContext::new();
    ctx.export_brep(&brep).unwrap();

    let plane_count = ctx
        .entities
        .iter()
        .filter(|e| e.type_name == "PLANE")
        .count();
    // A cube has 6 planar faces; some planes may be shared if
    // the deduplication catches them, but typically 6 distinct planes
    assert!(
        plane_count >= 3 && plane_count <= 6,
        "Expected 3-6 PLANEs for a cube, got {plane_count}"
    );
}

#[test]
fn export_brep_vertex_points_reference_cartesian_points() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut ctx = ExportContext::new();
    ctx.export_brep(&brep).unwrap();

    for entity in &ctx.entities {
        if entity.type_name == "VERTEX_POINT" {
            // VERTEX_POINT('', #ref_to_cartesian_point)
            assert_eq!(entity.params.len(), 2);
            assert!(matches!(entity.params[1], Parameter::EntityRef(_)));
        }
    }
}

#[test]
fn export_brep_with_inner_loop() {
    use popkern_core::constructors::builder::BRepBuilder;
    use popkern_core::geometry::Point3;

    let mut b = BRepBuilder::new(Tolerance::default());

    // Outer quad: v0-v1-v2-v3
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(4.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(4.0, 4.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 4.0, 0.0));

    // Inner quad (hole): v4-v5-v6-v7
    let v4 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v5 = b.add_vertex(Point3::new(3.0, 1.0, 0.0));
    let v6 = b.add_vertex(Point3::new(3.0, 3.0, 0.0));
    let v7 = b.add_vertex(Point3::new(1.0, 3.0, 0.0));

    // Outer loop edges
    let (_, he01, _) = b.add_edge(v0, v1);
    let (_, he12, _) = b.add_edge(v1, v2);
    let (_, he23, _) = b.add_edge(v2, v3);
    let (_, he30, _) = b.add_edge(v3, v0);

    // Inner loop edges (CW from outside = hole)
    let (_, he47, _) = b.add_edge(v4, v7);
    let (_, he76, _) = b.add_edge(v7, v6);
    let (_, he65, _) = b.add_edge(v6, v5);
    let (_, he54, _) = b.add_edge(v5, v4);

    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let face = b.add_face(shell);
    b.add_loop(&[he01, he12, he23, he30], face).unwrap();
    b.add_inner_loop(&[he47, he76, he65, he54], face).unwrap();

    let brep = b.build_unchecked();
    let mut ctx = ExportContext::new();
    let solid_id = ctx.export_brep(&brep).unwrap();

    let count =
        |name: &str| -> usize { ctx.entities.iter().filter(|e| e.type_name == name).count() };

    // Should have 1 FACE_OUTER_BOUND and 1 FACE_BOUND
    assert_eq!(count("FACE_OUTER_BOUND"), 1);
    assert_eq!(count("FACE_BOUND"), 1);
    // 2 edge loops (outer + inner)
    assert_eq!(count("EDGE_LOOP"), 2);
    assert_eq!(count("ADVANCED_FACE"), 1);

    let solid = ctx.entities.iter().find(|e| e.id == solid_id);
    assert!(solid.is_some());
}

#[test]
fn export_brep_edges_without_curves_synthesize_lines() {
    // unit_cube doesn't set curves, so the exporter should synthesize them
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut ctx = ExportContext::new();
    ctx.export_brep(&brep).unwrap();

    let line_count = ctx
        .entities
        .iter()
        .filter(|e| e.type_name == "LINE")
        .count();
    // All 12 edges should have synthesized LINE curves
    assert!(
        line_count >= 12,
        "Expected at least 12 LINEs, got {line_count}"
    );
}

#[test]
fn export_oriented_edge_has_correct_structure() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut ctx = ExportContext::new();
    ctx.export_brep(&brep).unwrap();

    for entity in &ctx.entities {
        if entity.type_name == "ORIENTED_EDGE" {
            // ORIENTED_EDGE('', *, *, #edge_curve, .T./.F.)
            assert_eq!(entity.params.len(), 5);
            assert_eq!(entity.params[0], Parameter::String("".into()));
            assert_eq!(entity.params[1], Parameter::Derived);
            assert_eq!(entity.params[2], Parameter::Derived);
            assert!(matches!(entity.params[3], Parameter::EntityRef(_)));
            assert!(
                entity.params[4] == Parameter::Enum("T".into())
                    || entity.params[4] == Parameter::Enum("F".into())
            );
        }
    }
}

#[test]
fn export_cube_has_surface_curves_and_pcurves() {
    // unit_cube has pcurves on all half-edges
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut ctx = ExportContext::new();
    ctx.export_brep(&brep).unwrap();

    let count =
        |name: &str| -> usize { ctx.entities.iter().filter(|e| e.type_name == name).count() };

    // Each of the 12 edges should have a SURFACE_CURVE wrapper
    assert_eq!(
        count("SURFACE_CURVE"),
        12,
        "Expected 12 SURFACE_CURVEs for cube"
    );

    // Each edge touches 2 faces -> 24 PCURVEs total
    assert_eq!(count("PCURVE"), 24, "Expected 24 PCURVEs for cube");

    // Each PCURVE needs a DEFINITIONAL_REPRESENTATION
    assert_eq!(
        count("DEFINITIONAL_REPRESENTATION"),
        24,
        "Expected 24 DEFINITIONAL_REPRESENTATIONs"
    );

    // Only one shared parametric context
    assert_eq!(
        count("PARAMETRIC_REPRESENTATION_CONTEXT"),
        1,
        "Expected 1 shared PARAMETRIC_REPRESENTATION_CONTEXT"
    );
}

#[test]
fn cylinder_export_has_seam_curve() {
    let tol = Tolerance::default();
    let cyl = BRep::cylinder(1.0, 2.0, tol).unwrap();
    let mut ctx = ExportContext::new();
    ctx.export_brep(&cyl).unwrap();

    let count =
        |name: &str| -> usize { ctx.entities.iter().filter(|e| e.type_name == name).count() };

    // A cylinder has a seam edge where both half-edges belong to the same
    // lateral face. That edge should use SEAM_CURVE, not SURFACE_CURVE.
    assert!(
        count("SEAM_CURVE") >= 1,
        "Expected at least 1 SEAM_CURVE for cylinder seam edge, got {}",
        count("SEAM_CURVE"),
    );
}

#[test]
fn cylinder_step_text_contains_seam_curve() {
    let tol = Tolerance::default();
    let cyl = BRep::cylinder(1.0, 2.0, tol).unwrap();
    let doc = crate::step::document::StepDocument {
        products: vec![StepProductNode {
            name: "Cyl".into(),
            id: "CYL001".into(),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![cyl],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: Default::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };

    let opts = StepExportOptions::default();
    let result = export_step(&doc, &opts).unwrap();
    let step_str = String::from_utf8(result).unwrap();
    assert!(
        step_str.contains("SEAM_CURVE"),
        "Cylinder STEP export should contain SEAM_CURVE for seam edge"
    );
}

#[test]
fn cube_has_no_seam_curves() {
    // A cube has no seam edges — all edges are shared between distinct faces.
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut ctx = ExportContext::new();
    ctx.export_brep(&brep).unwrap();

    let seam_count = ctx
        .entities
        .iter()
        .filter(|e| e.type_name == "SEAM_CURVE")
        .count();
    assert_eq!(
        seam_count, 0,
        "Cube should have 0 SEAM_CURVEs, got {seam_count}"
    );
}
