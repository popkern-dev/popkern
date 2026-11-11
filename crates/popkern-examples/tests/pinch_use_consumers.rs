//! A face touching itself at a vertex (a pinch use) through every
//! consumer (PK-1129, the consumer table of the tangency foundations
//! design): a 20 × 20 × 5 plate less a ball of radius 4 at (16, 10, 5),
//! whose dimple's rim touches the plate's top edge at (20, 10, 5). The
//! top face is one disc there, its boundary one loop through the vertex
//! twice. Each consumer carries it (validation, tessellation, measures,
//! STEP both ways, a boolean with it as an operand, an edge imprinted
//! from it into the wedge it enters) or declines with a
//! structured error (fillets and shells at the pinch), never serving a
//! body that fails validation.

use std::f64::consts::PI;
use std::sync::Arc;

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexNeighbourhood};
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step, import_step,
};

const PINCH: Point3 = Point3 {
    x: 20.0,
    y: 10.0,
    z: 5.0,
};

fn plate(lo: (f64, f64, f64), size: (f64, f64, f64), tol: Tolerance) -> BRep {
    let (x0, y0, z0) = lo;
    let (x1, y1, z1) = (x0 + size.0, y0 + size.1, z0 + size.2);
    let rect = |z: f64| {
        vec![
            Point3::new(x0, y0, z),
            Point3::new(x1, y0, z),
            Point3::new(x1, y1, z),
            Point3::new(x0, y1, z),
        ]
    };
    BRep::polygon_frustum(&rect(z0), &rect(z1), tol).expect("plate")
}

fn pinched(tol: Tolerance) -> BRep {
    let ball = BRep::sphere(4.0, tol)
        .and_then(|b| b.transform(&Transform3::translation(16.0, 10.0, 5.0)))
        .expect("ball");
    boolean(
        BoolOperation::Subtraction,
        plate((0.0, 0.0, 0.0), (20.0, 20.0, 5.0), tol),
        ball,
        tol,
        None,
    )
    .expect("the dimpled plate serves")
}

fn volume() -> f64 {
    2000.0 - 2.0 / 3.0 * PI * 64.0
}

fn area() -> f64 {
    2.0 * 400.0 + 4.0 * 100.0 - 16.0 * PI + 32.0 * PI
}

fn assert_measures(body: &BRep, tol: &Tolerance, what: &str) {
    let v = popkern_measure::body_volume(body, tol).expect("volume");
    let a = popkern_measure::body_surface_area(body, tol).expect("area");
    assert!((v - volume()).abs() < 1e-6 * volume(), "{what}: volume {v}");
    assert!((a - area()).abs() < 1e-6 * area(), "{what}: area {a}");
}

/// The body is valid, a manifold everywhere, with one pinch use: the top
/// face meets the vertex twice round it.
#[test]
fn validation_reads_the_pinch_as_a_manifold_vertex() {
    let tol = Tolerance::default();
    let body = pinched(tol);
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    let mut pinches = 0;
    for (v, n) in body.vertex_neighbourhoods() {
        assert!(n.is_manifold(), "{v:?} {n:?}");
        if let VertexNeighbourhood::Disc { pinch_uses } = &n
            && !pinch_uses.is_empty()
        {
            assert!(body.vertices[v].position.distance_to(&PINCH) < 1e-6);
            pinches += 1;
        }
    }
    assert_eq!(pinches, 1);
}

/// Tessellation meshes the face whole (its area is the face's) and the
/// measures are exact.
#[test]
fn tessellation_and_measures_carry_the_pinch() {
    let tol = Tolerance::default();
    let body = pinched(tol);
    assert_measures(&body, &tol, "served");
    let mesh = popkern_mesh::tessellate(&body, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let p = |i: usize| mesh.vertices[i].position;
    let mesh_area: f64 = mesh
        .triangles
        .iter()
        .map(|t| 0.5 * (p(t[1]) - p(t[0])).cross(&(p(t[2]) - p(t[0]))).length())
        .sum();
    assert!(
        (mesh_area - area()).abs() < 2e-3 * area(),
        "mesh area {mesh_area}"
    );
}

/// STEP carries the loop through the vertex twice both ways: the body
/// read back is valid with the same measures.
#[test]
fn step_round_trips_the_pinch() {
    let tol = Tolerance::default();
    let body = pinched(tol);
    let doc = StepDocument {
        products: vec![StepProductNode {
            name: "pinched".into(),
            id: "pinched".into(),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![body.clone()],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    let bytes = export_step(&doc, &StepExportOptions::default()).expect("export");
    let read = import_step(&bytes, &StepImportOptions::default()).expect("import");
    let solids: Vec<&BRep> = read
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter().flat_map(|s| s.solids.iter()))
        .collect();
    assert_eq!(solids.len(), 1);
    let back = solids[0];
    let errors = back.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    assert_eq!(back.faces.len(), body.faces.len());
    assert_measures(back, &tol, "read back");
}

/// A boolean with the pinched body as an operand, away from the pinch:
/// a block on the plate's far side unites with it.
#[test]
fn a_boolean_takes_the_pinched_body_as_an_operand() {
    let tol = Tolerance::default();
    let body = pinched(tol);
    let block = plate((-5.0, 5.0, 1.0), (10.0, 10.0, 3.0), tol);
    let united = boolean(BoolOperation::Union, body, block, tol, None).expect("union");
    let errors = united.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    // The block's half outside the plate adds 5 × 10 × 3.
    let v = popkern_measure::body_volume(&united, &tol).expect("volume");
    assert!((v - (volume() + 150.0)).abs() < 1e-6 * v, "{v}");
}

/// Fillets and chamfers of the edges at the pinch, and a shell, offset or
/// thicken of the body, do not handle a pinch use yet: each declines
/// with `CoreError::PinchUse` naming the face and the vertex, before
/// building anything (the shell built a skin whose Euler–Poincaré count
/// read 4 for 2 and declined only on validating it).
#[test]
fn fillets_chamfers_and_shells_decline_at_the_pinch() {
    use popkern_blend::{
        ChamferKind, ChamferOptions, ChamferSpec, FilletOptions, FilletRadius, FilletSpec,
        chamfer_edges, fillet_edges,
    };
    use popkern_core::error::CoreError;
    use popkern_shell::{ShellDirection, ShellError, ShellOptions, shell};
    let tol = Tolerance::default();
    let body = pinched(tol);
    let at_pinch = |e: &popkern_core::topology::Edge| {
        let he = &body.half_edges[e.half_edge];
        let mate = he.mate().and_then(|m| body.half_edges.get(m));
        let seam = mate.is_some_and(|m| body.loops[m.loop_].face == body.loops[he.loop_].face);
        let ends = [Some(he.origin), mate.map(|m| m.origin)];
        !seam
            && ends
                .iter()
                .flatten()
                .any(|&v| body.vertices[v].position.distance_to(&PINCH) < 1e-6)
    };
    let edges: Vec<_> = body
        .edges
        .iter()
        .filter(|(_, e)| at_pinch(e))
        .map(|(h, _)| h)
        .collect();
    assert!(!edges.is_empty());
    let pinch_named = |e: &CoreError| {
        matches!(e, CoreError::PinchUse { vertex, .. }
            if body.vertices[*vertex].position.distance_to(&PINCH) < 1e-6)
    };
    let fillets: Vec<FilletSpec> = edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.5),
        })
        .collect();
    match fillet_edges(&body, &fillets, &FilletOptions::default(), tol, None) {
        Err(popkern_blend::BlendError::Topology(e)) => assert!(pinch_named(&e), "{e}"),
        other => panic!("fillet: {:?}", other.map(|r| r.brep.faces.len())),
    }
    let chamfers: Vec<ChamferSpec> = edges
        .iter()
        .map(|&edge| ChamferSpec {
            edge,
            kind: ChamferKind::Equal(0.5),
        })
        .collect();
    match chamfer_edges(&body, &chamfers, &ChamferOptions::default(), tol, None) {
        Err(popkern_blend::BlendError::Topology(e)) => assert!(pinch_named(&e), "{e}"),
        other => panic!("chamfer: {:?}", other.map(|r| r.brep.faces.len())),
    }
    let shelled = shell(
        &body,
        0.5,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol,
        None,
    );
    match shelled {
        Err(ShellError::Core(e)) => assert!(pinch_named(&e), "{e}"),
        other => panic!("shell: {:?}", other.map(|r| r.brep.faces.len())),
    }
}

/// A quadratic arc on the plate's top from `a` through control `b` to `c`.
fn arc(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> popkern_core::geometry::curve::CurveOnEdge {
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
    use popkern_core::geometry::nurbs_curve::NurbsCurve;
    CurveOnEdge {
        geometry: CurveGeometry::Nurbs(NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: [a, b, c]
                .iter()
                .map(|&(x, y)| Point3::new(x, y, 5.0))
                .collect(),
            weights: vec![1.0; 3],
        }),
        t_start: 0.0,
        t_end: 1.0,
    }
}

/// An edge imprinted from the pinch into the top face enters one of
/// the two cusps between the dimple's rim and the plate's edge, both
/// leaving the vertex along the edge's line: an arc into either cusp,
/// less curved than the rim (0.1 against 0.25), to the far corner on
/// its side splits the face, the piece it cuts off its control
/// triangle's third (100 / 3), and is removed again; leaving into the
/// dimple, it enters no wedge and declines with the body unchanged.
#[test]
fn an_edge_is_imprinted_from_the_pinch_into_the_wedge_it_enters() {
    let tol = Tolerance::default();
    let body = pinched(tol);
    let vertex_at = |b: &BRep, p: Point3| {
        b.vertices
            .iter()
            .find(|(_, v)| v.position.distance_to(&p) < 1e-9)
            .map(|(h, _)| h)
            .expect("vertex")
    };
    let top = body
        .faces
        .keys()
        .find(|&f| {
            body.face_point_trim_containment(f, &Point3::new(2.0, 2.0, 5.0))
                .is_some_and(|c| c.distance < 1e-9)
        })
        .expect("top");
    let top_area = 400.0 - 16.0 * PI;
    let pinch = vertex_at(&body, PINCH);
    for (corner, control) in [((0.0, 20.0), (20.0, 20.0)), ((0.0, 0.0), (20.0, 0.0))] {
        let mut cut = body.clone();
        let far = vertex_at(&cut, Point3::new(corner.0, corner.1, 5.0));
        let made = cut
            .imprint_edge(top, pinch, far, arc((20.0, 10.0), control, corner))
            .unwrap_or_else(|e| panic!("into the cusp towards {corner:?}: {e}"));
        let errors = cut.validate_all();
        assert!(errors.is_empty(), "{corner:?}: {errors:?}");
        let new = made.new_face.expect("the face splits");
        let area = |f| popkern_measure::face_area(&cut, f, &tol).expect("area");
        let (a_new, a_top) = (area(new), area(top));
        let (small, large) = (a_new.min(a_top), a_new.max(a_top));
        assert!(
            (small - 100.0 / 3.0).abs() < 1e-9 * top_area,
            "{corner:?}: {small}"
        );
        assert!(
            (large - (top_area - 100.0 / 3.0)).abs() < 1e-9 * top_area,
            "{corner:?}: {large}"
        );
        assert_measures(&cut, &tol, "imprinted");
        cut.kill_imprinted_edge(made.edge).expect("the inverse");
        assert!(cut.validate_all().is_empty(), "{corner:?} removed");
        assert_eq!(cut.faces.len(), body.faces.len());
        assert_measures(&cut, &tol, "removed");
    }
    // Straight into the dimple: no wedge holds it.
    let mut into = body.clone();
    let corner = vertex_at(&into, Point3::new(0.0, 0.0, 5.0));
    let line = arc((20.0, 10.0), (10.0, 5.0), (0.0, 0.0));
    let declined = into.imprint_edge(top, pinch, corner, line);
    assert!(
        matches!(
            declined,
            Err(popkern_core::euler::EulerError::NotImprintable { .. })
        ),
        "{declined:?}"
    );
    assert_eq!(into.edges.len(), body.edges.len());
}
