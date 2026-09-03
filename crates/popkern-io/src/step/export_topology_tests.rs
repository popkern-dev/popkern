use super::*;
use crate::step::document::{
    OpaqueBag, RepresentationContext, ShapeRep, ShapeRepKind, StepProductNode, UnitContext,
};
use crate::step::export::{StepExportOptions, export_step};
use crate::step::export_topology_helpers::edge_half_edges;
use popkern_core::tolerance::Tolerance;
use std::sync::Arc;

#[test]
fn export_unit_cube_entity_counts() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut ctx = ExportContext::new();
    let solid_id = ctx.export_brep_items(&brep).unwrap()[0];

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
    ctx.export_brep_items(&brep).unwrap();

    let plane_count = ctx
        .entities
        .iter()
        .filter(|e| e.type_name == "PLANE")
        .count();
    // A cube has 6 planar faces; some planes may be shared if
    // the deduplication catches them, but typically 6 distinct planes
    assert!(
        (3..=6).contains(&plane_count),
        "Expected 3-6 PLANEs for a cube, got {plane_count}"
    );
}

#[test]
fn export_brep_vertex_points_reference_cartesian_points() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut ctx = ExportContext::new();
    ctx.export_brep_items(&brep).unwrap();

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
    let solid_id = ctx.export_brep_items(&brep).unwrap()[0];

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
    ctx.export_brep_items(&brep).unwrap();

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
    ctx.export_brep_items(&brep).unwrap();

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
    ctx.export_brep_items(&brep).unwrap();

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
        count("GEOMETRIC_REPRESENTATION_CONTEXT"),
        1,
        "Expected 1 shared 2D parametric context (GEOMETRIC_REPRESENTATION_CONTEXT complex)"
    );
}

#[test]
fn cylinder_export_has_seam_curve() {
    let tol = Tolerance::default();
    let cyl = BRep::cylinder(1.0, 2.0, tol).unwrap();
    let mut ctx = ExportContext::new();
    ctx.export_brep_items(&cyl).unwrap();

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

/// The seam's two pcurves must be BRANCH-SEPARATED (u values one
/// period apart — the UV trim contour cannot close otherwise) and
/// CO-PARAMETERIZED with the basis curve (both traced in the same
/// direction). Identical twins or opposed directions degrade
/// OpenCASCADE's trim of every closed-surface face (flagship boss:
/// the seam fillet sewed to its neighbours only along the chart
/// seam line).
#[test]
fn cylinder_seam_pcurves_are_branch_separated_and_co_parameterized() {
    use crate::step::parser::types::Parameter;
    fn walk_refs(p: &Parameter, out: &mut Vec<crate::step::parser::types::EntityId>) {
        match p {
            Parameter::EntityRef(id) => out.push(*id),
            Parameter::List(items) => {
                for item in items {
                    walk_refs(item, out);
                }
            }
            _ => {}
        }
    }
    fn refs_of(params: &[Parameter]) -> Vec<crate::step::parser::types::EntityId> {
        let mut out = Vec::new();
        for p in params {
            walk_refs(p, &mut out);
        }
        out
    }

    let tol = Tolerance::default();
    let cyl = BRep::cylinder(1.0, 2.0, tol).unwrap();
    let mut ctx = ExportContext::new();
    ctx.export_brep_items(&cyl).unwrap();

    let entity = |id: crate::step::parser::types::EntityId| {
        ctx.entities.iter().find(|e| e.id == id).unwrap()
    };
    let seam = ctx
        .entities
        .iter()
        .find(|e| e.type_name == "SEAM_CURVE")
        .expect("cylinder must export a SEAM_CURVE");
    let pcurve_ids: Vec<_> = refs_of(&seam.params)
        .into_iter()
        .filter(|&id| entity(id).type_name == "PCURVE")
        .collect();
    assert_eq!(pcurve_ids.len(), 2, "seam must carry exactly two pcurves");

    // PCURVE -> DEFINITIONAL_REPRESENTATION -> 2D curve -> CARTESIAN_POINTs.
    let cps_of = |pcurve_id: crate::step::parser::types::EntityId| -> Vec<(f64, f64)> {
        let dr = refs_of(&entity(pcurve_id).params)
            .into_iter()
            .find(|&id| entity(id).type_name == "DEFINITIONAL_REPRESENTATION")
            .expect("pcurve definitional representation");
        let curve = refs_of(&entity(dr).params)
            .into_iter()
            .find(|&id| entity(id).type_name.contains("CURVE"))
            .expect("pcurve basis 2D curve");
        refs_of(&entity(curve).params)
            .into_iter()
            .filter(|&id| entity(id).type_name == "CARTESIAN_POINT")
            .map(|id| {
                let vals: Vec<f64> = entity(id)
                    .params
                    .iter()
                    .flat_map(|p| match p {
                        Parameter::List(items) => items
                            .iter()
                            .filter_map(|v| match v {
                                Parameter::Real(r) => Some(*r),
                                Parameter::Int(n) => Some(*n as f64),
                                _ => None,
                            })
                            .collect::<Vec<f64>>(),
                        _ => Vec::new(),
                    })
                    .collect();
                (vals[0], vals[1])
            })
            .collect()
    };
    let a = cps_of(pcurve_ids[0]);
    let b = cps_of(pcurve_ids[1]);
    assert!(
        a.len() >= 2 && b.len() >= 2,
        "seam pcurves need control points"
    );
    let du = (a[0].0 - b[0].0).abs();
    assert!(
        (du - std::f64::consts::TAU).abs() < 1e-9,
        "seam pcurves must sit one period apart in u, got du = {du}"
    );
    let dir_a = (a.last().unwrap().1 - a[0].1).signum();
    let dir_b = (b.last().unwrap().1 - b[0].1).signum();
    assert_eq!(
        dir_a, dir_b,
        "seam pcurves must be co-parameterized with the basis curve"
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
    ctx.export_brep_items(&brep).unwrap();

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

#[cfg(toms_package)]
#[test]
fn open_boundary_brep_exports_as_surface_model() {
    use crate::step::document::StepDocument;
    use popkern_core::geometry::polyhedral_spline::{
        ControlNet, PolyhedralBrepExportOptions, build_polyhedral_spline_brep,
        construct_from_control_net,
    };

    let net = ControlNet::from_obj_str(include_str!(
        "../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T0.obj"
    ))
    .unwrap();
    let spline = construct_from_control_net(&net).unwrap().remove(0);
    let export = build_polyhedral_spline_brep(&spline, PolyhedralBrepExportOptions::default())
        .expect("PnS T0 BRep");
    assert!(export.brep.has_boundary_edges());

    let doc = StepDocument {
        products: vec![StepProductNode {
            name: "PnS T0".into(),
            id: "PNS-T0".into(),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::Other("SHAPE_REPRESENTATION".into()),
                solids: vec![export.brep],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: Default::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };

    let text =
        String::from_utf8(export_step(&doc, &StepExportOptions::default()).unwrap()).unwrap();
    assert!(text.contains("OPEN_SHELL"));
    assert!(text.contains("SHELL_BASED_SURFACE_MODEL"));
    assert!(!text.contains("MANIFOLD_SOLID_BREP"));
}

#[test]
fn collapsed_surface_trim_missing_curve_exports_degenerate_bspline() {
    let brep = collapsed_trim_sheet();
    let mut ctx = ExportContext::new();
    ctx.export_brep_items(&brep).unwrap();

    let count =
        |name: &str| -> usize { ctx.entities.iter().filter(|e| e.type_name == name).count() };

    assert_eq!(count("EDGE_CURVE"), 1);
    assert_eq!(count("B_SPLINE_CURVE_WITH_KNOTS"), 2);
    assert_eq!(
        count("LINE"),
        0,
        "collapsed trim must not export a zero-direction LINE"
    );
    assert_eq!(count("PCURVE"), 1);
}

#[test]
fn missing_edge_curve_synthesizes_from_valid_surface_trim() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = brep
        .edges
        .iter()
        .find_map(|(edge, _data)| {
            let has_surface_trim = brep
                .half_edges
                .iter()
                .any(|(_, half_edge)| half_edge.edge == edge && half_edge.pcurve.is_some());
            has_surface_trim.then_some(edge)
        })
        .expect("cube should expose edge surface trims");

    brep.edges[edge].curve = None;
    let curve = synthesize_missing_edge_curve_from_surface_trims(&brep, edge)
        .expect("valid surface trim should synthesize an export edge curve");
    brep.edges[edge].curve = Some(curve);

    let report = brep
        .measure_edge_pcurve_deviation(edge)
        .expect("synthesized edge should remain measurable against trims");
    assert!(
        report.max_deviation <= brep.tolerance.linear().max(1e-9),
        "synthesized edge curve must stay compatible with existing pcurves: {report:#?}"
    );
}

#[test]
fn missing_half_edge_pcurve_synthesizes_during_surface_curve_export() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = brep
        .edges
        .keys()
        .find(|&edge| {
            edge_half_edges(&brep, edge)
                .iter()
                .filter(|half_edge| brep.half_edges.get(**half_edge).is_some())
                .count()
                == 2
        })
        .expect("cube should have manifold edges");
    let half_edges = edge_half_edges(&brep, edge);
    brep.half_edges[half_edges[0]].pcurve = None;

    let mut ctx = ExportContext::new();
    ctx.export_brep_items(&brep).unwrap();

    let surface_curve = ctx
        .entities
        .iter()
        .find(|entity| entity.type_name == "SURFACE_CURVE")
        .expect("cube edge should export as SURFACE_CURVE");
    let pcurve_count = format!("{:?}", surface_curve.params[2])
        .matches("EntityRef")
        .count();
    assert_eq!(
        pcurve_count, 2,
        "missing face-use pcurve should be synthesized for STEP export: {surface_curve:?}"
    );
}

fn collapsed_trim_sheet() -> BRep {
    use popkern_core::geometry::pcurve::Pcurve;
    use popkern_core::geometry::{Point3, SurfaceGeometry, SurfaceOnFace, Vec3};
    use popkern_core::topology::{
        Edge, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Shell,
        ShellHandle, Solid, Vertex,
    };

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
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::zero(),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
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
