//! Transactional planned-boundary seam authoring for three-quad cap faces.
//!
//! The sphere-cap rung installs analytic cap surfaces after the topology has
//! already been split. At that point the outer cap seams must be authored as a
//! single edge triple: one 3-D curve plus both pcurves. Generic reprojection is
//! deliberately downstream of this module; it may preserve already-authored
//! analytic pcurves, but it does not get to invent the seam branch.

use super::three_quad_face_matching::{face_surface_matches, half_edge_endpoints_match_curve};
use crate::blend::boundary::{BlendBoundary, BoundarySegment};
use crate::fillet_surface::FilletSurfaceResult;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::fitting;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};
use popkern_math::PortableMath;

#[cfg(test)]
#[path = "three_quad_planned_seam_authoring_tests.rs"]
mod tests;

const MATCH_TOL: f64 = 1.0e-5;
const SAMPLE_COUNT: usize = 17;

/// Author the two planned outer-boundary seams owned by `subface_index` when
/// that sub-face has been replaced by an analytic sphere cap.
///
/// Postconditions for every committed seam:
/// - the edge's 3-D curve evaluates on the cap surface and neighbour surface;
/// - both half-edge pcurves evaluate to that same 3-D curve;
/// - the stored curve is oriented from the cap half-edge origin to its
///   destination;
/// - failure leaves the edge curve and both pcurves unchanged.
pub(super) fn author_planned_sphere_cap_boundary_seams(
    brep: &mut BRep,
    sub_face: FaceHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    subface_index: usize,
) -> Result<usize, String> {
    let Some(SurfaceGeometry::Sphere { .. }) = brep
        .faces
        .get(sub_face)
        .and_then(|f| f.surface.as_ref())
        .map(|sof| &sof.geometry)
    else {
        return Ok(0);
    };

    let mut authored = 0usize;
    for segment_index in [2 * subface_index, 2 * subface_index + 1] {
        let segment = boundary
            .segments
            .get(segment_index)
            .ok_or_else(|| format!("planned seam segment {segment_index} missing"))?;
        let he = find_cap_half_edge_for_segment(brep, sub_face, segment, fillet_surfaces)
            .ok_or_else(|| {
                format!(
                    "planned seam segment {segment_index} has no matching cap half-edge on {sub_face:?}"
                )
            })?;
        author_one_planned_seam(brep, he, sub_face, segment).map_err(|err| {
            format!(
                "planned seam segment {segment_index} on {sub_face:?} failed postconditions: {err}"
            )
        })?;
        authored += 1;
    }
    Ok(authored)
}

fn find_cap_half_edge_for_segment(
    brep: &BRep,
    sub_face: FaceHandle,
    segment: &BoundarySegment,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<HalfEdgeHandle> {
    let face_data = brep.faces.get(sub_face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..(brep.half_edges.len() + 1) {
        let he = brep.half_edges.get(cursor)?;
        let twin_face = brep.loops.get(brep.half_edges.get(he.twin)?.loop_)?.face;
        let matches = match segment {
            BoundarySegment::FaceEdge { face, curve } => {
                twin_face == *face
                    || half_edge_endpoints_match_curve(brep, cursor, curve, MATCH_TOL)
            }
            BoundarySegment::FilletEnd {
                fillet_index,
                curve,
            } => {
                half_edge_endpoints_match_curve(brep, cursor, curve, MATCH_TOL)
                    || fillet_surfaces.get(*fillet_index).is_some_and(|fillet| {
                        face_surface_matches(brep, twin_face, &fillet.surface)
                    })
            }
        };
        if matches {
            return Some(cursor);
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    None
}

fn author_one_planned_seam(
    brep: &mut BRep,
    cap_he: HalfEdgeHandle,
    cap_face: FaceHandle,
    segment: &BoundarySegment,
) -> Result<(), String> {
    let twin = brep
        .half_edges
        .get(cap_he)
        .map(|he| he.twin)
        .ok_or_else(|| format!("cap half-edge {cap_he:?} missing"))?;
    let edge = brep
        .half_edges
        .get(cap_he)
        .map(|he| he.edge)
        .ok_or_else(|| format!("cap half-edge {cap_he:?} missing edge"))?;
    let other_face = brep
        .loops
        .get(
            brep.half_edges
                .get(twin)
                .ok_or_else(|| format!("twin half-edge {twin:?} missing"))?
                .loop_,
        )
        .map(|loop_data| loop_data.face)
        .ok_or_else(|| format!("twin loop for {twin:?} missing"))?;
    if other_face == cap_face {
        return Err("cap seam twin remains on same face".to_string());
    }

    let cap_surface = brep
        .faces
        .get(cap_face)
        .and_then(|f| f.surface.as_ref())
        .map(|sof| sof.geometry.clone())
        .ok_or_else(|| format!("cap face {cap_face:?} has no surface"))?;
    let other_surface = brep
        .faces
        .get(other_face)
        .and_then(|f| f.surface.as_ref())
        .map(|sof| sof.geometry.clone())
        .ok_or_else(|| format!("neighbour face {other_face:?} has no surface"))?;

    let start = brep
        .vertices
        .get(brep.half_edges[cap_he].origin)
        .map(|v| v.position)
        .ok_or_else(|| format!("cap half-edge {cap_he:?} origin missing"))?;
    let end = brep
        .vertices
        .get(brep.half_edges[twin].origin)
        .map(|v| v.position)
        .ok_or_else(|| format!("twin half-edge {twin:?} origin missing"))?;

    let authored = build_planned_edge_triple(segment, &cap_surface, &other_surface, start, end)?;

    validate_edge_triple(&authored, &cap_surface, &other_surface, start, end)?;

    let old_edge_curve = brep.edges.get(edge).and_then(|e| e.curve.clone());
    let old_cap_pcurve = brep.half_edges.get(cap_he).and_then(|he| he.pcurve.clone());
    let old_other_pcurve = brep.half_edges.get(twin).and_then(|he| he.pcurve.clone());

    brep.edges[edge].curve = Some(authored.curve);
    brep.half_edges[cap_he].pcurve = Some(authored.cap_pcurve);
    brep.half_edges[twin].pcurve = Some(authored.other_pcurve);

    if let Err(err) = validate_committed_edge(brep, cap_he, &cap_surface, &other_surface) {
        brep.edges[edge].curve = old_edge_curve;
        brep.half_edges[cap_he].pcurve = old_cap_pcurve;
        brep.half_edges[twin].pcurve = old_other_pcurve;
        return Err(err);
    }

    Ok(())
}

struct AuthoredSeam {
    curve: CurveOnEdge,
    cap_pcurve: Pcurve,
    other_pcurve: Pcurve,
}

fn build_planned_edge_triple(
    segment: &BoundarySegment,
    cap_surface: &SurfaceGeometry,
    other_surface: &SurfaceGeometry,
    start: Point3,
    end: Point3,
) -> Result<AuthoredSeam, String> {
    if let (
        SurfaceGeometry::Sphere { center, radius, .. },
        SurfaceGeometry::Plane { origin, normal, .. },
    ) = (cap_surface, other_surface)
        && matches!(segment, BoundarySegment::FaceEdge { .. })
    {
        return build_sphere_plane_intersection_triple(
            *center,
            *radius,
            *origin,
            *normal,
            cap_surface,
            other_surface,
            start,
            end,
        );
    }

    let curve = oriented_nurbs_curve_from_segment(segment, start, end)
        .ok_or_else(|| "planned segment could not be reparameterized to half-edge".to_string())?;
    let curve_on_edge = CurveOnEdge {
        geometry: CurveGeometry::Nurbs(curve.clone()),
        t_start: 0.0,
        t_end: 1.0,
    };
    let cap_pcurve = sample_pcurve_on_surface(&curve, cap_surface, 0.0, 1.0)
        .ok_or_else(|| "planned segment did not project to cap surface".to_string())?;
    let other_pcurve = sample_pcurve_on_surface(&curve, other_surface, 0.0, 1.0)
        .ok_or_else(|| "planned segment did not project to neighbour surface".to_string())?;
    Ok(AuthoredSeam {
        curve: curve_on_edge,
        cap_pcurve,
        other_pcurve,
    })
}

#[allow(clippy::too_many_arguments)]
fn build_sphere_plane_intersection_triple(
    sphere_center: Point3,
    sphere_radius: f64,
    plane_origin: Point3,
    plane_normal: popkern_core::geometry::Vec3,
    cap_surface: &SurfaceGeometry,
    other_surface: &SurfaceGeometry,
    start: Point3,
    end: Point3,
) -> Result<AuthoredSeam, String> {
    let normal = plane_normal
        .normalized()
        .ok_or_else(|| "plane normal is degenerate".to_string())?;
    let signed_distance = (sphere_center - plane_origin).dot(&normal);
    let radius_sq = sphere_radius * sphere_radius - signed_distance * signed_distance;
    if radius_sq <= 0.0 {
        return Err(format!(
            "sphere-plane seam is empty or tangent: sphere radius {sphere_radius:.6e}, \
             plane offset {signed_distance:.6e}"
        ));
    }
    let circle_radius = radius_sq.sqrt();
    let circle_center = sphere_center - normal * signed_distance;
    let start_vec = start - circle_center;
    let end_vec = end - circle_center;
    let start_radial = start_vec.length();
    let end_radial = end_vec.length();
    let radial_tol = circle_radius.max(1.0) * 1.0e-4;
    if (start_radial - circle_radius).abs() > radial_tol
        || (end_radial - circle_radius).abs() > radial_tol
    {
        return Err(format!(
            "sphere-plane endpoints are not on intersection circle: start radial {:.6e}, \
             end radial {:.6e}, expected {circle_radius:.6e}",
            start_radial, end_radial
        ));
    }
    let u_dir = start_vec
        .normalized()
        .ok_or_else(|| "sphere-plane seam start is at circle centre".to_string())?;
    let v_dir = normal
        .cross(&u_dir)
        .normalized()
        .ok_or_else(|| "sphere-plane seam basis is degenerate".to_string())?;
    let end_unit = end_vec
        .normalized()
        .ok_or_else(|| "sphere-plane seam end is at circle centre".to_string())?;
    let mut end_angle = end_unit.dot(&v_dir).patan2(end_unit.dot(&u_dir));
    while end_angle > std::f64::consts::PI {
        end_angle -= std::f64::consts::TAU;
    }
    while end_angle < -std::f64::consts::PI {
        end_angle += std::f64::consts::TAU;
    }

    let mut points = Vec::with_capacity(SAMPLE_COUNT);
    let mut params = Vec::with_capacity(SAMPLE_COUNT);
    for i in 0..SAMPLE_COUNT {
        let s = i as f64 / (SAMPLE_COUNT - 1) as f64;
        let angle = s * end_angle;
        points.push(
            circle_center
                + u_dir * (circle_radius * angle.pcos())
                + v_dir * (circle_radius * angle.psin()),
        );
        params.push(s);
    }
    let curve = fitting::interpolate_curve(&points, &params, 5)
        .map_err(|e| format!("sphere-plane seam interpolation failed: {e}"))?;
    let cap_pcurve = sample_pcurve_on_surface(&curve, cap_surface, 0.0, 1.0)
        .ok_or_else(|| "sphere-plane seam did not project to sphere cap".to_string())?;
    let other_pcurve = sample_pcurve_on_surface(&curve, other_surface, 0.0, 1.0)
        .ok_or_else(|| "sphere-plane seam did not project to support plane".to_string())?;
    Ok(AuthoredSeam {
        curve: CurveOnEdge {
            geometry: CurveGeometry::Nurbs(curve),
            t_start: 0.0,
            t_end: 1.0,
        },
        cap_pcurve,
        other_pcurve,
    })
}

fn oriented_nurbs_curve_from_segment(
    segment: &BoundarySegment,
    start: Point3,
    end: Point3,
) -> Option<NurbsCurve> {
    let source = match segment {
        BoundarySegment::FilletEnd { curve, .. } | BoundarySegment::FaceEdge { curve, .. } => curve,
    };
    let (t0, t1) = source.parameter_range();
    let p0 = source.point_at(t0);
    let p1 = source.point_at(t1);
    let forward = start.distance_to(&p0) + end.distance_to(&p1)
        <= start.distance_to(&p1) + end.distance_to(&p0);

    let mut points = Vec::with_capacity(SAMPLE_COUNT);
    let mut params = Vec::with_capacity(SAMPLE_COUNT);
    for i in 0..SAMPLE_COUNT {
        let s = i as f64 / (SAMPLE_COUNT - 1) as f64;
        let source_s = if forward { s } else { 1.0 - s };
        let t = t0 + source_s * (t1 - t0);
        points.push(source.point_at(t));
        params.push(s);
    }
    fitting::interpolate_curve(&points, &params, 5).ok()
}

fn sample_pcurve_on_surface(
    curve: &NurbsCurve,
    surface: &SurfaceGeometry,
    t_start: f64,
    t_end: f64,
) -> Option<Pcurve> {
    let mut points = Vec::with_capacity(SAMPLE_COUNT);
    let mut params = Vec::with_capacity(SAMPLE_COUNT);
    for i in 0..SAMPLE_COUNT {
        let s = i as f64 / (SAMPLE_COUNT - 1) as f64;
        let t = t_start + s * (t_end - t_start);
        let point = curve.point_at(t);
        let (u, v) = project_point_to_uv(surface, &point)?;
        points.push(Point3::new(u, v, 0.0));
        params.push(t);
    }
    fitting::interpolate_curve(&points, &params, 5)
        .ok()
        .map(Pcurve::from_nurbs)
}

fn validate_edge_triple(
    seam: &AuthoredSeam,
    cap_surface: &SurfaceGeometry,
    other_surface: &SurfaceGeometry,
    start: Point3,
    end: Point3,
) -> Result<(), String> {
    let p_start = seam.curve.geometry.point_at(seam.curve.t_start);
    let p_end = seam.curve.geometry.point_at(seam.curve.t_end);
    let endpoint_tol = 1.0e-5;
    if p_start.distance_to(&start) > endpoint_tol || p_end.distance_to(&end) > endpoint_tol {
        return Err(format!(
            "edge curve endpoints do not match half-edge vertices: start {:.3e}, end {:.3e}",
            p_start.distance_to(&start),
            p_end.distance_to(&end)
        ));
    }
    validate_sampled_agreement(
        &seam.curve,
        &seam.cap_pcurve,
        &seam.other_pcurve,
        cap_surface,
        other_surface,
    )
}

fn validate_committed_edge(
    brep: &BRep,
    cap_he: HalfEdgeHandle,
    cap_surface: &SurfaceGeometry,
    other_surface: &SurfaceGeometry,
) -> Result<(), String> {
    let twin = brep.half_edges[cap_he].twin;
    let edge = brep.half_edges[cap_he].edge;
    let curve = brep.edges[edge]
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge:?} lost 3-D curve during seam commit"))?;
    let cap_pc = brep.half_edges[cap_he]
        .pcurve
        .as_ref()
        .ok_or_else(|| format!("cap half-edge {cap_he:?} lost pcurve during seam commit"))?;
    let other_pc = brep.half_edges[twin]
        .pcurve
        .as_ref()
        .ok_or_else(|| format!("twin half-edge {twin:?} lost pcurve during seam commit"))?;
    validate_sampled_agreement(curve, cap_pc, other_pc, cap_surface, other_surface)
}

fn validate_sampled_agreement(
    curve: &CurveOnEdge,
    cap_pcurve: &Pcurve,
    other_pcurve: &Pcurve,
    cap_surface: &SurfaceGeometry,
    other_surface: &SurfaceGeometry,
) -> Result<(), String> {
    let tolerance = 2.0e-4;
    let mut worst_cap = 0.0_f64;
    let mut worst_other = 0.0_f64;
    let mut worst_between = 0.0_f64;
    for i in 0..SAMPLE_COUNT {
        let s = i as f64 / (SAMPLE_COUNT - 1) as f64;
        let t = curve.t_start + s * (curve.t_end - curve.t_start);
        let p_curve = curve.geometry.point_at(t);
        let cap_uv = cap_pcurve.point_at(t);
        let other_uv = other_pcurve.point_at(t);
        let p_cap = cap_surface.point_at(cap_uv.x, cap_uv.y);
        let p_other = other_surface.point_at(other_uv.x, other_uv.y);
        worst_cap = worst_cap.max(p_curve.distance_to(&p_cap));
        worst_other = worst_other.max(p_curve.distance_to(&p_other));
        worst_between = worst_between.max(p_cap.distance_to(&p_other));
    }
    if worst_cap > tolerance || worst_other > tolerance || worst_between > tolerance {
        return Err(format!(
            "sample agreement exceeded tolerance: curve-cap {worst_cap:.3e}, \
             curve-neighbour {worst_other:.3e}, cap-neighbour {worst_between:.3e}"
        ));
    }
    Ok(())
}
