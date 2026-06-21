use std::collections::HashMap;

use super::readiness::{PatchEdgeRole, validate_brep_readiness};
use super::types::{PatchEdge, PolyhedralSpline, PolyhedralSplineError, PolyhedralSplineFamily};
use crate::geometry::Point3;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use crate::tolerance::Tolerance;
use crate::topology::{
    BRep, Edge, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Shell, ShellHandle,
    Solid, TopologyRole, Vertex, VertexHandle,
};
use slotmap::Key;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PolyhedralBrepExportOptions {
    pub tolerance: Tolerance,
    pub samples: usize,
    pub c0_tolerance: f64,
    pub operation: OperationId,
}

impl Default for PolyhedralBrepExportOptions {
    fn default() -> Self {
        let tolerance = Tolerance::default();
        Self {
            tolerance,
            samples: 11,
            c0_tolerance: tolerance.linear(),
            operation: OperationId(0),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolyhedralBrepFaceProvenance {
    pub patch_index: usize,
    pub face: FaceHandle,
}

#[derive(Debug, Clone)]
pub struct PolyhedralBrepExport {
    pub brep: BRep,
    pub faces: Vec<PolyhedralBrepFaceProvenance>,
    pub boundary_edge_count: usize,
    pub seam_edge_count: usize,
    pub collapsed_singularity_count: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct VertexKey(u64, u64, u64);

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct EdgeKey(VertexHandle, VertexHandle);

pub fn build_polyhedral_spline_brep(
    spline: &PolyhedralSpline,
    options: PolyhedralBrepExportOptions,
) -> Result<PolyhedralBrepExport, PolyhedralSplineError> {
    let readiness = validate_brep_readiness(spline, options.samples, options.c0_tolerance)?;
    let mut brep = BRep::new(options.tolerance);
    let solid = brep.solids.insert(Solid {
        shell: ShellHandle::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: FaceHandle::default(),
        solid,
    });
    brep.solids[solid].shell = shell;

    let mut vertices = HashMap::new();
    let mut edges = HashMap::new();
    let mut faces = Vec::with_capacity(spline.patches.len());
    let mut collapsed_singularity_count = 0usize;

    for (patch_index, patch) in spline.patches.iter().enumerate() {
        let face = insert_patch_face(&mut brep, shell, options.operation, patch_index);
        if brep.shells[shell].face == FaceHandle::default() {
            brep.shells[shell].face = face;
        }
        brep.faces[face].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(patch.to_bezier_nurbs()),
            orientation: true,
        });
        let mut loop_half_edges = Vec::with_capacity(4);
        for boundary in patch_boundary() {
            let start = patch.evaluate(boundary.uv_start.0, boundary.uv_start.1);
            let end = patch.evaluate(boundary.uv_end.0, boundary.uv_end.1);
            let curve = boundary_curve_on_edge(patch, boundary);
            let start_vertex = intern_vertex(&mut brep, &mut vertices, start, options.c0_tolerance);
            let end_vertex = intern_vertex(&mut brep, &mut vertices, end, options.c0_tolerance);
            let edge_role = readiness
                .edge_uses
                .iter()
                .find(|edge_use| {
                    edge_use.patch_index == patch_index && edge_use.edge == boundary.edge
                })
                .map(|edge_use| edge_use.role)
                .ok_or_else(|| export_error(spline, "missing edge role in readiness report"))?;
            if matches!(edge_role, PatchEdgeRole::CollapsedSingularity { .. }) {
                collapsed_singularity_count += 1;
            }
            let half_edge = intern_half_edge(
                &mut brep,
                &mut edges,
                EdgeBuild {
                    start_vertex,
                    end_vertex,
                    curve,
                    role: edge_role,
                    family: spline.family,
                },
            )?;
            brep.half_edges[half_edge].pcurve =
                Some(Pcurve::linear(boundary.uv_start, boundary.uv_end, 0.0, 1.0));
            loop_half_edges.push(half_edge);
        }
        wire_loop(&mut brep, face, &loop_half_edges)?;
        faces.push(PolyhedralBrepFaceProvenance { patch_index, face });
    }

    brep.validate()
        .map_err(|err| export_error(spline, format!("post-build validation failed: {err}")))?;
    Ok(PolyhedralBrepExport {
        boundary_edge_count: brep.boundary_edges.len(),
        seam_edge_count: edges.len().saturating_sub(brep.boundary_edges.len()),
        collapsed_singularity_count,
        brep,
        faces,
    })
}

#[derive(Debug, Clone, Copy)]
struct BoundaryEdge {
    edge: PatchEdge,
    uv_start: (f64, f64),
    uv_end: (f64, f64),
}

struct EdgeBuild {
    start_vertex: VertexHandle,
    end_vertex: VertexHandle,
    curve: Option<CurveOnEdge>,
    role: PatchEdgeRole,
    family: PolyhedralSplineFamily,
}

fn patch_boundary() -> [BoundaryEdge; 4] {
    [
        BoundaryEdge {
            edge: PatchEdge::VMin,
            uv_start: (0.0, 0.0),
            uv_end: (1.0, 0.0),
        },
        BoundaryEdge {
            edge: PatchEdge::UMax,
            uv_start: (1.0, 0.0),
            uv_end: (1.0, 1.0),
        },
        BoundaryEdge {
            edge: PatchEdge::VMax,
            uv_start: (1.0, 1.0),
            uv_end: (0.0, 1.0),
        },
        BoundaryEdge {
            edge: PatchEdge::UMin,
            uv_start: (0.0, 1.0),
            uv_end: (0.0, 0.0),
        },
    ]
}

fn insert_patch_face(
    brep: &mut BRep,
    shell: ShellHandle,
    operation: OperationId,
    patch_index: usize,
) -> FaceHandle {
    let name = brep.name_pool.intern(
        EntityName {
            operation,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: patch_index as u32,
            },
            ancestors: Vec::new(),
        },
        None,
    );
    brep.faces.insert(Face {
        outer_loop: LoopHandle::default(),
        inner_loops: Vec::new(),
        shell,
        surface: None,
        name: Some(name),
        tolerance: None,
        role: TopologyRole::Production,
    })
}

fn intern_vertex(
    brep: &mut BRep,
    vertices: &mut HashMap<VertexKey, VertexHandle>,
    point: Point3,
    tolerance: f64,
) -> VertexHandle {
    for (&_, &handle) in vertices.iter() {
        if brep.vertices[handle].position.distance_to(&point) <= tolerance {
            return handle;
        }
    }
    let key = vertex_key(point, tolerance);
    *vertices.entry(key).or_insert_with(|| {
        brep.vertices.insert(Vertex {
            position: point,
            outgoing: HalfEdgeHandle::default(),
            name: None,
            tolerance: None,
        })
    })
}

fn intern_half_edge(
    brep: &mut BRep,
    edges: &mut HashMap<EdgeKey, HalfEdgeHandle>,
    build: EdgeBuild,
) -> Result<HalfEdgeHandle, PolyhedralSplineError> {
    let key = edge_key(build.start_vertex, build.end_vertex);
    let collapsed = matches!(build.role, PatchEdgeRole::CollapsedSingularity { .. });
    if !collapsed && let Some(&other) = edges.get(&key) {
        let edge = brep.half_edges[other].edge;
        if !matches!(
            build.role,
            PatchEdgeRole::Seam { .. } | PatchEdgeRole::CollapsedSingularity { .. }
        ) {
            return Err(PolyhedralSplineError::BrepExport {
                family: build.family,
                reason: "boundary edge unexpectedly reused".to_string(),
            });
        }
        let half_edge = brep.half_edges.insert(HalfEdge {
            origin: build.start_vertex,
            twin: other,
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge,
            loop_: LoopHandle::default(),
            pcurve: None,
        });
        brep.half_edges[other].twin = half_edge;
        return Ok(half_edge);
    }

    let edge = brep.insert_edge_tracked(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: build.curve,
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: TopologyRole::Production,
    });
    let half_edge = brep.half_edges.insert(HalfEdge {
        origin: build.start_vertex,
        twin: HalfEdgeHandle::default(),
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge,
        loop_: LoopHandle::default(),
        pcurve: None,
    });
    brep.edges[edge].half_edge = half_edge;
    if brep.vertices[build.start_vertex].outgoing == HalfEdgeHandle::default() {
        brep.vertices[build.start_vertex].outgoing = half_edge;
    }
    if matches!(
        build.role,
        PatchEdgeRole::Boundary | PatchEdgeRole::CollapsedSingularity { .. }
    ) {
        brep.boundary_edges.insert(edge);
    }
    if !collapsed {
        edges.insert(key, half_edge);
    }
    Ok(half_edge)
}

fn boundary_curve_on_edge(
    patch: &super::types::PolyhedralPatch,
    boundary: BoundaryEdge,
) -> Option<CurveOnEdge> {
    let (degree, control_points) = match boundary.edge {
        PatchEdge::VMin => (
            patch.degree_u,
            (0..=patch.degree_u)
                .map(|u| patch.cp(u, 0))
                .collect::<Vec<_>>(),
        ),
        PatchEdge::UMax => (
            patch.degree_v,
            (0..=patch.degree_v)
                .map(|v| patch.cp(patch.degree_u, v))
                .collect::<Vec<_>>(),
        ),
        PatchEdge::VMax => (
            patch.degree_u,
            (0..=patch.degree_u)
                .rev()
                .map(|u| patch.cp(u, patch.degree_v))
                .collect::<Vec<_>>(),
        ),
        PatchEdge::UMin => (
            patch.degree_v,
            (0..=patch.degree_v)
                .rev()
                .map(|v| patch.cp(0, v))
                .collect::<Vec<_>>(),
        ),
    };
    if control_points
        .first()
        .zip(control_points.last())
        .is_some_and(|(start, end)| start.distance_to(end) <= 1.0e-12)
    {
        return None;
    }
    Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(NurbsCurve {
            degree,
            knots: bezier_curve_knots(degree),
            weights: vec![1.0; control_points.len()],
            control_points,
        }),
        t_start: 0.0,
        t_end: 1.0,
    })
}

fn bezier_curve_knots(degree: usize) -> Vec<f64> {
    let mut knots = vec![0.0; degree + 1];
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    knots
}

fn wire_loop(
    brep: &mut BRep,
    face: FaceHandle,
    half_edges: &[HalfEdgeHandle],
) -> Result<(), PolyhedralSplineError> {
    let loop_handle = brep.loops.insert(Loop {
        half_edge: half_edges[0],
        face,
    });
    let len = half_edges.len();
    for (index, &half_edge) in half_edges.iter().enumerate() {
        brep.half_edges[half_edge].loop_ = loop_handle;
        brep.half_edges[half_edge].next = half_edges[(index + 1) % len];
        brep.half_edges[half_edge].prev = half_edges[(index + len - 1) % len];
    }
    brep.faces[face].outer_loop = loop_handle;
    Ok(())
}

fn edge_key(a: VertexHandle, b: VertexHandle) -> EdgeKey {
    if a.data().as_ffi() <= b.data().as_ffi() {
        EdgeKey(a, b)
    } else {
        EdgeKey(b, a)
    }
}

fn vertex_key(point: Point3, tolerance: f64) -> VertexKey {
    let scale = 1.0 / tolerance.max(1.0e-12);
    VertexKey(
        (point.x * scale).round().to_bits(),
        (point.y * scale).round().to_bits(),
        (point.z * scale).round().to_bits(),
    )
}

fn export_error(spline: &PolyhedralSpline, reason: impl Into<String>) -> PolyhedralSplineError {
    PolyhedralSplineError::BrepExport {
        family: spline.family,
        reason: reason.into(),
    }
}
