use super::*;
use crate::euler::DiagnosticHandle;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::gregory::{GregoryPatch, RectangularGregory};
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, VertexHandle};

#[test]
fn live_three_face_tangency_names_unsupported_support_gap() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let face = incident_face(&cube, vertex);
    cube.faces[face].surface.as_mut().expect("surface").geometry = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.25,
    };
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(
        report.strategy.family,
        CornerApexStrategyFamily::ThreeFaceBlend
    );
    assert_eq!(
        report.topology.surface_compatibility,
        SurfaceCompatibility::MissingPrimitive("three-face-tangency-torus-support")
    );
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive("three-face-tangency-torus-support")
    );
    assert_ne!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive(
            "three-face-tangency-nonplanar-support"
        )
    );
}

#[test]
fn live_three_face_tangency_names_nurbs_support_gap() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let face = incident_face(&cube, vertex);
    cube.faces[face].surface.as_mut().expect("surface").geometry =
        SurfaceGeometry::Nurbs(planar_nurbs_support());
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(
        report.topology.surface_compatibility,
        SurfaceCompatibility::MissingPrimitive("three-face-tangency-nurbs-support")
    );
    assert_ne!(
        report.topology.surface_compatibility,
        SurfaceCompatibility::Full
    );
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive("three-face-tangency-nurbs-support")
    );
    assert_ne!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive(
            "three-face-tangency-nonplanar-support"
        )
    );
}

#[test]
fn live_three_face_tangency_names_gregory_support_as_conditional() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let face = incident_face(&cube, vertex);
    cube.faces[face].surface.as_mut().expect("surface").geometry =
        SurfaceGeometry::Gregory(GregoryPatch::Rectangular(Box::new(flat_gregory_support())));
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(
        report.topology.surface_compatibility,
        SurfaceCompatibility::Conditional("nested-blend-support")
    );
    assert_ne!(
        report.topology.surface_compatibility,
        SurfaceCompatibility::Full
    );
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive(
            "three-face-tangency-gregory-support"
        )
    );
    assert_ne!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive(
            "three-face-tangency-nonplanar-support"
        )
    );
}

#[test]
fn live_three_face_tangency_accepts_analytic_sphere_support() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let face = incident_face(&cube, vertex);
    let center = match &cube.faces[face].surface.as_ref().expect("surface").geometry {
        SurfaceGeometry::Plane { normal, .. } => {
            cube.vertices[vertex].position + normal.normalized().expect("normal") * 1.0
        }
        _ => panic!("fixture face is planar"),
    };
    cube.faces[face].surface.as_mut().expect("surface").geometry = SurfaceGeometry::Sphere {
        center,
        radius: 1.0,
    };
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));
    assert_eq!(
        report.strategy.family,
        CornerApexStrategyFamily::ThreeFaceBlend
    );
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::Ready
    );
}

#[test]
fn live_three_face_tangency_accepts_analytic_cylinder_support() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let face = incident_face(&cube, vertex);
    let (origin, axis) = match &cube.faces[face].surface.as_ref().expect("surface").geometry {
        SurfaceGeometry::Plane { normal, .. } => {
            let normal = normal.normalized().expect("normal");
            (
                cube.vertices[vertex].position + normal,
                perpendicular(normal),
            )
        }
        _ => panic!("fixture face is planar"),
    };
    cube.faces[face].surface.as_mut().expect("surface").geometry = SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius: 1.0,
    };
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(
        report.strategy.family,
        CornerApexStrategyFamily::ThreeFaceBlend
    );
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::Ready
    );
}

#[test]
fn live_three_face_tangency_accepts_analytic_cone_support() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let face = incident_face(&cube, vertex);
    let normal = match &cube.faces[face].surface.as_ref().expect("surface").geometry {
        SurfaceGeometry::Plane { normal, .. } => normal.normalized().expect("normal"),
        _ => panic!("fixture face is planar"),
    };
    let position = cube.vertices[vertex].position;
    cube.faces[face].surface.as_mut().expect("surface").geometry =
        cone_through_vertex_with_normal(position, normal);
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        BlendKind::FilletConstantRadius,
        1_000_000,
    ));

    assert_eq!(
        report.strategy.family,
        CornerApexStrategyFamily::ThreeFaceBlend
    );
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::Ready
    );
}

#[test]
fn live_three_face_tangency_accepts_two_cylinder_one_plane_supports() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let faces = incident_faces(&cube, vertex);
    cube.faces[faces[0]]
        .surface
        .as_mut()
        .expect("surface")
        .geometry = SurfaceGeometry::Cylinder {
        origin: Point3::new(1.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    cube.faces[faces[1]]
        .surface
        .as_mut()
        .expect("surface")
        .geometry = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 1.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    cube.faces[faces[2]]
        .surface
        .as_mut()
        .expect("surface")
        .geometry = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(
        report.strategy.family,
        CornerApexStrategyFamily::ThreeFaceBlend
    );
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::Ready
    );
}

#[test]
fn live_mutual_trim_names_nonlinear_spine_gap() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);
    cube.edges[edges[0]].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        t_start: 0.0,
        t_end: std::f64::consts::FRAC_PI_2,
    });

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges[..2],
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(report.strategy.family, CornerApexStrategyFamily::MutualTrim);
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive("nonlinear-spine-mutual-trim-probe")
    );
}

fn request(
    vertex: VertexHandle,
    edges: &[EdgeHandle],
    kind: BlendKind,
    radius_nanometers: i64,
) -> CornerApexDiagnosticRequest {
    CornerApexDiagnosticRequest {
        scope: LocalCornerScope {
            vertex: DiagnosticHandle::Vertex(vertex),
        },
        selected_edges: edges
            .iter()
            .map(|edge| LocalBlendEdgeSelection {
                edge: DiagnosticHandle::Edge(*edge),
                blend_kind: kind,
                radius_nanometers: Some(radius_nanometers),
            })
            .collect(),
        selection_intent: "unit-test",
        operation_context: "live-corner-apex-diagnostics",
    }
}

fn vertex_at(brep: &BRep, x: f64, y: f64, z: f64) -> VertexHandle {
    brep.vertices
        .iter()
        .find(|(_, vertex)| {
            vertex.position.x == x && vertex.position.y == y && vertex.position.z == z
        })
        .map(|(handle, _)| handle)
        .expect("vertex")
}

fn incident_edges(brep: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    let mut edges = std::collections::BTreeSet::new();
    for (_, half_edge) in &brep.half_edges {
        if half_edge.origin == vertex
            || brep
                .half_edges
                .get(half_edge.twin)
                .is_some_and(|twin| twin.origin == vertex)
        {
            edges.insert(half_edge.edge);
        }
    }
    edges.into_iter().collect()
}

fn incident_face(brep: &BRep, vertex: VertexHandle) -> crate::topology::FaceHandle {
    incident_faces(brep, vertex)
        .into_iter()
        .next()
        .expect("incident face")
}

fn incident_faces(brep: &BRep, vertex: VertexHandle) -> Vec<crate::topology::FaceHandle> {
    brep.faces
        .iter()
        .filter(|(_, face)| loop_vertices(brep, face.outer_loop).contains(&vertex))
        .map(|(handle, _)| handle)
        .collect()
}

fn loop_vertices(brep: &BRep, loop_: crate::topology::LoopHandle) -> Vec<VertexHandle> {
    let mut vertices = Vec::new();
    let start = brep.loops[loop_].half_edge;
    let mut cursor = start;
    loop {
        let half_edge = &brep.half_edges[cursor];
        vertices.push(half_edge.origin);
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    vertices
}

fn perpendicular(normal: Vec3) -> Vec3 {
    let reference = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    normal
        .cross(&reference)
        .normalized()
        .expect("perpendicular")
}

fn cone_through_vertex_with_normal(vertex: Point3, normal: Vec3) -> SurfaceGeometry {
    let half_angle = std::f64::consts::FRAC_PI_4;
    let tangent = perpendicular(normal);
    let axis = (tangent * half_angle.cos() - normal * half_angle.sin())
        .normalized()
        .expect("cone axis");
    let radial_direction = (normal + axis * half_angle.sin())
        .normalized()
        .expect("cone radial direction");
    let apex = vertex - radial_direction * half_angle.tan() - axis;
    SurfaceGeometry::Cone {
        apex,
        axis,
        half_angle,
    }
}

fn planar_nurbs_support() -> NurbsSurface {
    NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0, 1.0, 1.0, 1.0],
    }
}

fn flat_gregory_support() -> RectangularGregory {
    let mut boundary_cps = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for (u_index, row) in boundary_cps.iter_mut().enumerate() {
        for (v_index, point) in row.iter_mut().enumerate() {
            *point = Point3::new(u_index as f64 / 3.0, v_index as f64 / 3.0, 0.0);
        }
    }
    RectangularGregory {
        boundary_cps,
        inner_cps: [(Point3::new(0.5, 0.5, 0.0), Point3::new(0.5, 0.5, 0.0)); 4],
    }
}
