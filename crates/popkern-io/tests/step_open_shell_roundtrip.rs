use std::sync::Arc;

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{Point3, SurfaceGeometry, SurfaceOnFace};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, BrepExportQualityOptions, BrepTopologyOutputPolicy, BrepTrimRepairOptions, Edge,
    EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Shell, ShellHandle,
    Solid, Vertex, analyze_brep_export_quality,
};
use popkern_io::{
    HealingLevel, OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind,
    StepDocument, StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step,
    import_step,
};

#[test]
fn open_shell_sheet_body_round_trips_through_step() {
    let brep = BRep::open_box(Tolerance::default()).unwrap();
    let document = document_with_body(brep.clone());
    let exported = export_step(&document, &StepExportOptions::default()).unwrap();
    let text = String::from_utf8_lossy(&exported);
    assert!(text.contains("SHELL_BASED_SURFACE_MODEL"));
    assert!(text.contains("OPEN_SHELL"));

    let imported = import_step(&exported, &import_options()).unwrap();
    assert!(
        !imported.diagnostics.has_errors(),
        "STEP import diagnostics: {:?}",
        imported.diagnostics.entries()
    );
    let round_tripped = first_body(&imported).expect("expected one imported sheet body");

    assert_eq!(round_tripped.faces.len(), brep.faces.len());
    assert_eq!(round_tripped.edges.len(), brep.edges.len());
    assert_eq!(
        round_tripped.boundary_edges.len(),
        brep.boundary_edges.len()
    );
    round_tripped.validate().unwrap();
}

#[test]
fn collapsed_trim_degenerate_edge_curve_reimports_from_step() {
    let brep = collapsed_trim_sheet();
    let exported = export_step(&document_with_body(brep), &StepExportOptions::default()).unwrap();
    let text = String::from_utf8_lossy(&exported);
    assert!(text.contains("SHELL_BASED_SURFACE_MODEL"));
    assert!(text.contains("B_SPLINE_CURVE_WITH_KNOTS"));
    assert!(!text.contains("LINE("));

    let imported = import_step(&exported, &import_options()).unwrap();
    assert!(
        !imported.diagnostics.has_errors(),
        "STEP import diagnostics: {:?}",
        imported.diagnostics.entries()
    );
    let round_tripped = first_body(&imported).expect("expected one imported sheet body");
    assert_eq!(round_tripped.edges.len(), 1);
    assert_eq!(round_tripped.boundary_edges.len(), 1);

    let edge = round_tripped.edges.values().next().unwrap();
    assert!(
        matches!(
            edge.curve.as_ref().map(|curve| &curve.geometry),
            Some(CurveGeometry::Nurbs(_))
        ),
        "collapsed trim should reimport as a real degenerate NURBS edge curve"
    );

    let mut options = open_sheet_quality_options();
    options.require_edge_curves = true;
    options.validate_brep = false;
    let report = analyze_brep_export_quality(round_tripped, options);
    assert!(report.is_ready(), "export-quality report: {report:?}");
    assert_eq!(report.missing_edge_curves, 0);
    assert_eq!(report.edge_curve_blockers, 0);
}

fn document_with_body(brep: BRep) -> StepDocument {
    StepDocument {
        products: vec![StepProductNode {
            name: "sheet".into(),
            id: "sheet".into(),
            children: Vec::new(),
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::Other("SHAPE_REPRESENTATION".into()),
                solids: vec![brep],
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    }
}

fn import_options() -> StepImportOptions {
    StepImportOptions {
        healing: HealingLevel::None,
        ..Default::default()
    }
}

fn open_sheet_quality_options() -> BrepExportQualityOptions {
    BrepExportQualityOptions {
        topology_output: BrepTopologyOutputPolicy::OpenSheet,
        trim: BrepTrimRepairOptions {
            endpoint_tolerance: 1e-6,
            sampled_tolerance: 1e-6,
            samples: 5,
            validate_output: true,
        },
        require_edge_curves: false,
        allow_collapsed_trim_edge_curves: false,
        max_pcurve_deviation: None,
        max_edge_surface_deviation: None,
        validate_brep: true,
    }
}

fn first_body(result: &popkern_io::StepImportResult) -> Option<&BRep> {
    result
        .document
        .products
        .iter()
        .flat_map(|product| &product.shapes)
        .flat_map(|shape| &shape.solids)
        .next()
}

fn collapsed_trim_sheet() -> BRep {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: ShellHandle::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: FaceHandle::default(),
        solid,
    });
    brep.solids[solid].shell = shell;
    let vertex = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: HalfEdgeHandle::default(),
        name: None,
        tolerance: None,
    });
    let face = brep.faces.insert(Face {
        outer_loop: LoopHandle::default(),
        inner_loops: Vec::new(),
        shell,
        surface: Some(SurfaceOnFace {
            geometry: collapsed_edge_surface(),
            orientation: true,
        }),
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    let loop_ = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face,
    });
    let half_edge = brep.half_edges.insert(HalfEdge {
        origin: vertex,
        twin: HalfEdgeHandle::default(),
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge: EdgeHandle::default(),
        loop_,
        pcurve: Some(Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0)),
    });
    let edge = brep.edges.insert(Edge {
        half_edge,
        curve: None,
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: Default::default(),
    });
    brep.half_edges[half_edge].edge = edge;
    brep.half_edges[half_edge].next = half_edge;
    brep.half_edges[half_edge].prev = half_edge;
    brep.vertices[vertex].outgoing = half_edge;
    brep.loops[loop_].half_edge = half_edge;
    brep.faces[face].outer_loop = loop_;
    brep.shells[shell].face = face;
    brep.boundary_edges.insert(edge);
    brep
}

fn collapsed_edge_surface() -> SurfaceGeometry {
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0, 1.0, 1.0, 1.0],
    })
}
