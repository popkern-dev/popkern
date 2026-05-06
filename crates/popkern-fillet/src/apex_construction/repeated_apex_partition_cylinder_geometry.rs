use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveOnEdge, Point3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::fillet_surface::FilletSurfaceResult;

use super::super::error::ApexBuildError;
use super::super::repeated_apex_loop_plan::RepeatedApexLoopPlan;

pub(super) fn author_missing_on_cylinder_line_trim_for_endpoint(
    brep: &mut BRep,
    face: FaceHandle,
    point: Point3,
) -> Result<Option<EdgeHandle>, ApexBuildError> {
    let Some(half_edge) = containing_on_cylinder_line_chord_half_edge(brep, face, point) else {
        return Ok(None);
    };
    let edge = brep
        .half_edges
        .get(half_edge)
        .map(|half_edge_data| half_edge_data.edge)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author missing cylinder trim curve: \
                 half-edge {half_edge:?} is stale"
            ),
        })?;
    let start = brep
        .half_edges
        .get(half_edge)
        .and_then(|half_edge_data| brep.vertices.get(half_edge_data.origin))
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author missing cylinder trim curve: \
                 half-edge {half_edge:?} has stale start vertex"
            ),
        })?;
    let end_vertex =
        half_edge_end_vertex(brep, half_edge).ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author missing cylinder trim curve: \
                 half-edge {half_edge:?} has stale twin"
            ),
        })?;
    let end = brep
        .vertices
        .get(end_vertex)
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author missing cylinder trim curve: \
                 end vertex {end_vertex:?} is stale"
            ),
        })?;
    let curve =
        CurveOnEdge::line_segment(start, end, brep.tolerance.linear()).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot author missing cylinder trim curve: \
                 half-edge {half_edge:?} has degenerate endpoints"
                ),
            }
        })?;
    brep.edges[edge].curve = Some(curve);
    Ok(Some(edge))
}

fn containing_on_cylinder_line_chord_half_edge(
    brep: &BRep,
    face: FaceHandle,
    point: Point3,
) -> Option<HalfEdgeHandle> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        let start_point = brep.vertices.get(half_edge.origin)?.position;
        let end_vertex = half_edge_end_vertex(brep, cursor)?;
        let end_point = brep.vertices.get(end_vertex)?.position;
        if point_lies_on_cylinder_face(brep, face, start_point)
            && point_lies_on_cylinder_face(brep, face, end_point)
            && point_lies_on_segment(point, start_point, end_point, brep.tolerance.linear())
        {
            return Some(cursor);
        }
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    None
}

fn point_lies_on_segment(point: Point3, start: Point3, end: Point3, tolerance: f64) -> bool {
    let edge = end - start;
    let len2 = edge.dot(&edge);
    if len2 <= tolerance.max(1.0e-12).powi(2) {
        return false;
    }
    let t = (point - start).dot(&edge) / len2;
    if t <= tolerance || t >= 1.0 - tolerance {
        return false;
    }
    let projected = start + edge * t;
    point.distance_to(&projected) <= tolerance.max(edge.length() * 1.0e-8)
}

pub(super) fn point_is_not_edge_endpoint(brep: &BRep, edge: EdgeHandle, point: Point3) -> bool {
    let Some(edge_data) = brep.edges.get(edge) else {
        return false;
    };
    let Some(half_edge) = brep.half_edges.get(edge_data.half_edge) else {
        return false;
    };
    let Some(twin) = brep.half_edges.get(half_edge.twin) else {
        return false;
    };
    let start = brep.vertices[half_edge.origin].position;
    let end = brep.vertices[twin.origin].position;
    let tol = brep.tolerance.linear();
    point.distance_to(&start) > tol && point.distance_to(&end) > tol
}

pub(super) fn point_lies_on_cylinder_face(brep: &BRep, face: FaceHandle, point: Point3) -> bool {
    let Some(surface) = brep
        .faces
        .get(face)
        .and_then(|face_data| face_data.surface.as_ref())
    else {
        return false;
    };
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = surface.geometry
    else {
        return false;
    };
    let axis_len2 = axis.dot(&axis);
    if axis_len2 <= 1.0e-24 {
        return false;
    }
    let rel = point - origin;
    let axial = rel.dot(&axis) / axis_len2;
    let center = origin + axis * axial;
    let radial = point.distance_to(&center);
    (radial - radius).abs() <= brep.tolerance.linear().max(radius.abs() * 1.0e-6)
}

pub(super) fn cylinder_radial_vector(
    point: Point3,
    origin: Point3,
    axis: popkern_core::geometry::Vec3,
) -> Option<popkern_core::geometry::Vec3> {
    let axis_len2 = axis.dot(&axis);
    if axis_len2 <= 1.0e-24 {
        return None;
    }
    let rel = point - origin;
    let axial = rel.dot(&axis) / axis_len2;
    let center = origin + axis * axial;
    Some(point - center)
}

pub(super) fn target_cylinder_face_for_fillet(
    brep: &BRep,
    fillet_index: usize,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<FaceHandle, ApexBuildError> {
    let target = fillet_surfaces
        .get(fillet_index)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!("missing fillet surface for FilletEnd[{fillet_index}]"),
        })?;
    brep.faces
        .iter()
        .filter_map(|(face_handle, face)| {
            let surface = face.surface.as_ref()?;
            cylinder_surfaces_match(&surface.geometry, &target.surface).then_some(face_handle)
        })
        .next()
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] has no live \
                 cylinder face matching its fillet surface"
            ),
        })
}

fn cylinder_surfaces_match(a: &SurfaceGeometry, b: &SurfaceGeometry) -> bool {
    let (
        SurfaceGeometry::Cylinder {
            origin: a_origin,
            axis: a_axis,
            radius: a_radius,
        },
        SurfaceGeometry::Cylinder {
            origin: b_origin,
            axis: b_axis,
            radius: b_radius,
        },
    ) = (a, b)
    else {
        return false;
    };
    if (a_radius - b_radius).abs() > a_radius.abs().max(b_radius.abs()) * 1.0e-5 + 1.0e-9 {
        return false;
    }
    let Some(a_axis_unit) = a_axis.normalized() else {
        return false;
    };
    let Some(b_axis_unit) = b_axis.normalized() else {
        return false;
    };
    if a_axis_unit.dot(&b_axis_unit).abs() < 1.0 - 1.0e-5 {
        return false;
    }
    let delta = *a_origin - *b_origin;
    let axial = delta.dot(&b_axis_unit);
    let perpendicular = delta - b_axis_unit * axial;
    perpendicular.length() <= 1.0e-5
}

pub(super) fn half_edge_end_vertex(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<VertexHandle> {
    let twin = brep.half_edges.get(half_edge)?.twin;
    brep.half_edges.get(twin).map(|twin_data| twin_data.origin)
}

pub(super) fn arc_half_edges_between(
    loop_plan: &RepeatedApexLoopPlan,
    start_pos: usize,
    end_pos: usize,
) -> Vec<HalfEdgeHandle> {
    let mut out = Vec::new();
    let n = loop_plan.half_edges.len();
    if n == 0 {
        return out;
    }
    let mut pos = start_pos % n;
    loop {
        out.push(loop_plan.half_edges[pos]);
        pos = (pos + 1) % n;
        if pos == end_pos {
            break;
        }
        if out.len() > n {
            break;
        }
    }
    out
}
pub(super) struct CylinderTrimSurface {
    origin: Point3,
    axis: popkern_core::geometry::Vec3,
}

impl CylinderTrimSurface {
    pub(super) fn from_face(
        brep: &BRep,
        face: FaceHandle,
        fillet_index: usize,
    ) -> Result<Self, ApexBuildError> {
        let Some(surface) = brep
            .faces
            .get(face)
            .and_then(|face_data| face_data.surface.as_ref())
        else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] target cylinder \
                     face {face:?} has no surface"
                ),
            });
        };
        let SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius: _,
        } = surface.geometry
        else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] target face \
                     {face:?} is not a cylinder"
                ),
            });
        };
        Ok(Self { origin, axis })
    }

    pub(super) fn validate_generator_segment(
        &self,
        brep: &BRep,
        fillet_index: usize,
        start: VertexHandle,
        end: VertexHandle,
    ) -> Result<(), ApexBuildError> {
        let start_point = brep.vertices[start].position;
        let end_point = brep.vertices[end].position;
        let start_radial =
            cylinder_radial_vector(start_point, self.origin, self.axis).ok_or_else(|| {
                ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] cannot classify \
                     generator trim {start:?}->{end:?}: degenerate cylinder axis"
                    ),
                }
            })?;
        let end_radial =
            cylinder_radial_vector(end_point, self.origin, self.axis).ok_or_else(|| {
                ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] cannot classify \
                     generator trim {start:?}->{end:?}: degenerate cylinder axis"
                    ),
                }
            })?;
        let Some(start_unit) = start_radial.normalized() else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] generator trim \
                     {start:?}->{end:?} has degenerate start radial vector"
                ),
            });
        };
        let Some(end_unit) = end_radial.normalized() else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] generator trim \
                     {start:?}->{end:?} has degenerate end radial vector"
                ),
            });
        };
        if start_unit.dot(&end_unit) < 1.0 - 1.0e-5 {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] generator trim \
                     {start:?}->{end:?} does not keep a constant cylinder radial direction"
                ),
            });
        }
        Ok(())
    }
}
