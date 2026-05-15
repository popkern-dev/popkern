//! Final construction audit for bounded-surface boundary evidence.

use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::edge_triple::{
    BoundarySide, EdgeAnchors, EdgeHint, EdgeTripleOptions, SurfaceSide, compute_edge_geometry,
};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{
    BRep, FaceHandle, HalfEdgeHandle, SurfaceBoundaryHintDecline, SurfaceBoundaryHintDeclineReason,
};

pub(crate) fn audit_gregory_analytic_boundary_evidence(brep: &mut BRep) {
    let mut half_edges = Vec::new();
    for (_, face) in &brep.faces {
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        if !matches!(surface.geometry, SurfaceGeometry::Gregory(_)) {
            continue;
        }
        let start = brep.loops[face.outer_loop].half_edge;
        let mut cursor = start;
        for _ in 0..=brep.half_edges.len() {
            half_edges.push(cursor);
            cursor = brep.half_edges[cursor].next;
            if cursor == start {
                break;
            }
        }
    }
    for half_edge in half_edges {
        audit_half_edge(brep, half_edge);
    }
}

fn audit_half_edge(brep: &mut BRep, half_edge: HalfEdgeHandle) {
    if brep.surface_boundary_hint(half_edge).is_some()
        || brep.surface_boundary_decline(half_edge).is_some()
    {
        return;
    }
    if !has_analytic_neighbor(brep, half_edge) {
        return;
    }
    if try_author_exact_boundary_owned_triple(brep, half_edge) {
        return;
    }
    let Some(pcurve) = brep.half_edges[half_edge].pcurve.clone() else {
        brep.set_surface_boundary_decline(
            half_edge,
            SurfaceBoundaryHintDecline {
                reason: SurfaceBoundaryHintDeclineReason::MissingPcurve,
                best_boundary_deviation: f64::INFINITY,
                best_boundary: None,
                endpoint_gap: 0.0,
            },
        );
        return;
    };
    let classification = classify_rectangular_boundary_pcurve(&pcurve);
    if let Some(boundary) = classification.boundary {
        brep.set_surface_boundary_hint(half_edge, boundary);
    } else {
        brep.set_surface_boundary_decline(
            half_edge,
            SurfaceBoundaryHintDecline {
                reason: SurfaceBoundaryHintDeclineReason::NotOnCanonicalBoundary,
                best_boundary_deviation: classification.best_deviation,
                best_boundary: classification.boundary,
                endpoint_gap: 0.0,
            },
        );
    }
}

fn try_author_exact_boundary_owned_triple(brep: &mut BRep, half_edge: HalfEdgeHandle) -> bool {
    let Some((cap_face, other_face, twin)) = adjacent_faces(brep, half_edge) else {
        return false;
    };
    let Some(cap_surface) = brep
        .faces
        .get(cap_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| surface.geometry.clone())
    else {
        return false;
    };
    if !matches!(cap_surface, SurfaceGeometry::Gregory(_)) {
        return false;
    }
    let Some(other_surface) = brep
        .faces
        .get(other_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| surface.geometry.clone())
    else {
        return false;
    };
    if !is_supported_analytic_surface(&other_surface) {
        return false;
    }
    let Some(start) = brep
        .half_edges
        .get(half_edge)
        .and_then(|data| brep.vertices.get(data.origin).map(|vertex| vertex.position))
    else {
        return false;
    };
    let Some(end) = brep
        .half_edges
        .get(twin)
        .and_then(|data| brep.vertices.get(data.origin).map(|vertex| vertex.position))
    else {
        return false;
    };
    let boundaries = [
        BoundarySide::UMin,
        BoundarySide::UMax,
        BoundarySide::VMin,
        BoundarySide::VMax,
    ];
    for boundary in boundaries {
        let hint = EdgeHint::SurfaceBoundary {
            side: SurfaceSide::A,
            boundary,
        };
        let Ok(triple) = compute_edge_geometry(
            &cap_surface,
            &other_surface,
            EdgeAnchors::Open { start, end },
            Some(hint),
            &brep.tolerance,
            &EdgeTripleOptions::default(),
        ) else {
            continue;
        };
        let Some(edge) = brep.half_edges.get(half_edge).map(|data| data.edge) else {
            return false;
        };
        brep.edges[edge].curve = Some(CurveOnEdge {
            geometry: triple.curve_3d,
            t_start: triple.t_start,
            t_end: triple.t_end,
        });
        brep.half_edges[half_edge].pcurve = Some(triple.pcurve_a);
        brep.half_edges[twin].pcurve = Some(triple.pcurve_b);
        return brep.set_surface_boundary_hint(half_edge, boundary);
    }
    false
}

fn adjacent_faces(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<(FaceHandle, FaceHandle, HalfEdgeHandle)> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let cap_face = brep.loops.get(half_edge_data.loop_)?.face;
    let twin = half_edge_data.twin;
    let other_face = brep.loops.get(brep.half_edges.get(twin)?.loop_)?.face;
    (cap_face != other_face).then_some((cap_face, other_face, twin))
}

fn has_analytic_neighbor(brep: &BRep, half_edge: HalfEdgeHandle) -> bool {
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return false;
    };
    let Some(twin) = brep.half_edges.get(half_edge_data.twin) else {
        return false;
    };
    let Some(twin_face) = brep.loops.get(twin.loop_).map(|loop_data| loop_data.face) else {
        return false;
    };
    brep.faces
        .get(twin_face)
        .and_then(|face| face.surface.as_ref())
        .is_some_and(|surface| is_supported_analytic_surface(&surface.geometry))
}

fn is_supported_analytic_surface(surface: &SurfaceGeometry) -> bool {
    matches!(
        surface,
        SurfaceGeometry::Plane { .. }
            | SurfaceGeometry::Cylinder { .. }
            | SurfaceGeometry::Sphere { .. }
            | SurfaceGeometry::Torus { .. }
    )
}

struct BoundaryPcurveClassification {
    boundary: Option<BoundarySide>,
    best_deviation: f64,
}

fn classify_rectangular_boundary_pcurve(pcurve: &Pcurve) -> BoundaryPcurveClassification {
    let (t0, t1) = pcurve.parameter_range();
    let mut max_u_min = 0.0_f64;
    let mut max_u_max = 0.0_f64;
    let mut max_v_min = 0.0_f64;
    let mut max_v_max = 0.0_f64;
    for i in 0..9 {
        let s = i as f64 / 8.0;
        let uv = pcurve.point_at(t0 + s * (t1 - t0));
        max_u_min = max_u_min.max(uv.x.abs());
        max_u_max = max_u_max.max((1.0 - uv.x).abs());
        max_v_min = max_v_min.max(uv.y.abs());
        max_v_max = max_v_max.max((1.0 - uv.y).abs());
    }
    let candidates = [
        (BoundarySide::UMin, max_u_min),
        (BoundarySide::UMax, max_u_max),
        (BoundarySide::VMin, max_v_min),
        (BoundarySide::VMax, max_v_max),
    ];
    let Some(best) = candidates.into_iter().min_by(|a, b| a.1.total_cmp(&b.1)) else {
        return BoundaryPcurveClassification {
            boundary: None,
            best_deviation: f64::INFINITY,
        };
    };
    BoundaryPcurveClassification {
        boundary: (best.1 <= 1.0e-6).then_some(best.0),
        best_deviation: best.1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::curve::CurveOnEdge;
    use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregory};
    use popkern_core::geometry::{Point3, SurfaceOnFace, Vec3};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::{
        Edge, Face, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Shell, ShellHandle, Solid, Vertex,
    };

    #[test]
    fn audit_authors_exact_triple_for_canonical_gregory_boundary() {
        let mut brep = gregory_plane_boundary_fixture();
        let edge = brep.edges.keys().next().unwrap();
        let half_edge = brep.edges[edge].half_edge;
        brep.half_edges[half_edge].pcurve =
            Some(Pcurve::linear((0.25, 0.0), (0.25, 1.0), 0.0, 1.0));

        audit_gregory_analytic_boundary_evidence(&mut brep);

        assert_eq!(
            brep.surface_boundary_hint(half_edge),
            Some(BoundarySide::UMin)
        );
        assert!(brep.surface_boundary_decline(half_edge).is_none());
        let pcurve = brep.half_edges[half_edge].pcurve.as_ref().unwrap();
        let (t0, t1) = pcurve.parameter_range();
        assert_eq!(pcurve.point_at(t0).x, 0.0);
        assert_eq!(pcurve.point_at(t1).x, 0.0);
    }

    fn gregory_plane_boundary_fixture() -> BRep {
        let mut brep = BRep::new(Tolerance::default());
        let solid = brep.solids.insert(Solid {
            shell: ShellHandle::default(),
        });
        let shell = brep.shells.insert(Shell {
            face: Default::default(),
            solid,
        });
        brep.solids[solid].shell = shell;
        let start = Point3::new(0.0, 0.0, 0.0);
        let end = Point3::new(0.0, 1.0, 0.0);
        let va = brep.vertices.insert(vertex(start));
        let vb = brep.vertices.insert(vertex(end));
        let face_a = brep.faces.insert(face(shell, gregory_surface(start, end)));
        let face_b = brep.faces.insert(face(shell, plane_surface()));
        let loop_a = brep.loops.insert(Loop {
            half_edge: HalfEdgeHandle::default(),
            face: face_a,
        });
        let loop_b = brep.loops.insert(Loop {
            half_edge: HalfEdgeHandle::default(),
            face: face_b,
        });
        let he_a = brep.half_edges.insert(half_edge(va, loop_a));
        let he_b = brep.half_edges.insert(half_edge(vb, loop_b));
        let edge = brep.edges.insert(Edge {
            half_edge: he_a,
            curve: CurveOnEdge::line_segment(
                start + Vec3::new(0.0, 0.0, 0.1),
                end + Vec3::new(0.0, 0.0, 0.1),
                1e-12,
            ),
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });
        brep.half_edges[he_a].edge = edge;
        brep.half_edges[he_b].edge = edge;
        brep.half_edges[he_a].twin = he_b;
        brep.half_edges[he_b].twin = he_a;
        brep.half_edges[he_a].next = he_a;
        brep.half_edges[he_a].prev = he_a;
        brep.half_edges[he_b].next = he_b;
        brep.half_edges[he_b].prev = he_b;
        brep.loops[loop_a].half_edge = he_a;
        brep.loops[loop_b].half_edge = he_b;
        brep.faces[face_a].outer_loop = loop_a;
        brep.faces[face_b].outer_loop = loop_b;
        brep.vertices[va].outgoing = he_a;
        brep.vertices[vb].outgoing = he_b;
        brep.shells[shell].face = face_a;
        brep
    }

    fn vertex(position: Point3) -> Vertex {
        Vertex {
            position,
            outgoing: HalfEdgeHandle::default(),
            name: None,
            tolerance: None,
        }
    }

    fn half_edge(origin: popkern_core::topology::VertexHandle, loop_: LoopHandle) -> HalfEdge {
        HalfEdge {
            origin,
            twin: HalfEdgeHandle::default(),
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge: Default::default(),
            loop_,
            pcurve: None,
        }
    }

    fn face(shell: ShellHandle, surface: SurfaceOnFace) -> Face {
        Face {
            outer_loop: LoopHandle::default(),
            inner_loops: Vec::new(),
            shell,
            surface: Some(surface),
            name: None,
            tolerance: None,
            role: Default::default(),
        }
    }

    fn gregory_surface(start: Point3, end: Point3) -> SurfaceOnFace {
        let lerp = |s: f64| start + (end - start) * s;
        let offset = Vec3::new(1.0, 0.0, 0.0);
        let u0 = [lerp(0.0), lerp(1.0 / 3.0), lerp(2.0 / 3.0), lerp(1.0)];
        let row = |scale: f64| {
            [
                u0[0] + offset * scale,
                u0[1] + offset * scale,
                u0[2] + offset * scale,
                u0[3] + offset * scale,
            ]
        };
        SurfaceOnFace {
            geometry: SurfaceGeometry::Gregory(GregoryPatch::Rectangular(Box::new(
                RectangularGregory {
                    boundary_cps: [u0, row(1.0 / 3.0), row(2.0 / 3.0), row(1.0)],
                    inner_cps: [(row(1.0 / 3.0)[1], row(1.0 / 3.0)[1]); 4],
                },
            ))),
            orientation: true,
        }
    }

    fn plane_surface() -> SurfaceOnFace {
        SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
            orientation: true,
        }
    }
}
