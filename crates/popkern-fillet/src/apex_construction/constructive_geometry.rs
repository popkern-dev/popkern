//! Geometry and surface attachment helpers for constructive apex building.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::fillet_surface::FilletSurfaceResult;

use super::error::ApexBuildError;
use super::plan::ApexPlan;

pub(super) fn planned_fillet_index_for_sector(
    plan: &ApexPlan,
    sector_index: usize,
) -> Result<usize, ApexBuildError> {
    let edge_index = 2 * sector_index;
    let edge = plan
        .target_topology
        .apex_edges
        .get(edge_index)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "FilletEnd sector {sector_index} references missing apex_edge[{edge_index}]"
            ),
        })?;
    match edge.role {
        super::target_topology::TargetEdgeRole::CapVsCyl { fillet_idx } => Ok(fillet_idx),
        ref other => Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "FilletEnd sector {sector_index} expected CapVsCyl apex_edge[{edge_index}], got {other:?}"
            ),
        }),
    }
}

pub(super) fn attach_fillet_end_curve_to_edge(
    brep: &mut BRep,
    edge: EdgeHandle,
    planned_curve: &NurbsCurve,
    fillet_index: usize,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<(), ApexBuildError> {
    let p = planned_curve.degree;
    let n = planned_curve.control_points.len() - 1;
    let t_start = planned_curve.knots[p];
    let t_end = planned_curve.knots[n + 1];
    let p_start = planned_curve.point_at(t_start);
    let p_end = planned_curve.point_at(t_end);

    let he =
        brep.edges
            .get(edge)
            .map(|e| e.half_edge)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!("attach_fillet_end_curve_to_edge: edge {edge:?} missing"),
            })?;
    let twin = brep.half_edges[he].twin;
    let he_origin = brep.vertices[brep.half_edges[he].origin].position;
    let twin_origin = brep.vertices[brep.half_edges[twin].origin].position;
    let surface = fillet_surfaces
        .get(fillet_index)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!("missing fillet surface for FilletEnd[{fillet_index}]"),
        })?
        .surface
        .clone();
    if let Some((geometry, arc_start, arc_end)) =
        cylinder_arc_through_edge_endpoints(&surface, he_origin, twin_origin)
    {
        if std::env::var_os("POPKERN_A2_TRACE").is_some() {
            eprintln!(
                "[a2-trace §14.9/fillet-end-arc] constructing cylinder-owned FilletEnd[{fillet_index}] \
                 edge curve on {edge:?}"
            );
        }
        brep.edges[edge].curve = Some(CurveOnEdge {
            geometry,
            t_start: arc_start,
            t_end: arc_end,
        });
        return Ok(());
    }

    let forward_cost = he_origin.distance_to(&p_start) + twin_origin.distance_to(&p_end);
    let reversed_cost = he_origin.distance_to(&p_end) + twin_origin.distance_to(&p_start);
    if forward_cost.min(reversed_cost) <= 1.0e-5 {
        return attach_curve_to_edge_aligned(brep, edge, planned_curve);
    }

    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "FilletEnd[{fillet_index}] planned curve endpoint mismatch and no cylinder-owned \
             analytic arc is available for edge {edge:?}"
        ),
    })
}

pub(super) fn cylinder_arc_through_edge_endpoints(
    surface: &SurfaceGeometry,
    start: Point3,
    end: Point3,
) -> Option<(CurveGeometry, f64, f64)> {
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = surface
    else {
        return None;
    };
    let axis_len2 = axis.dot(axis);
    if axis_len2 <= 1.0e-24 || *radius <= 0.0 {
        return None;
    }
    let v0 = (start - *origin).dot(axis) / axis_len2;
    let v1 = (end - *origin).dot(axis) / axis_len2;
    if (v0 - v1).abs() > 1.0e-5 {
        return None;
    }
    let v_const = 0.5 * (v0 + v1);
    let center = *origin + *axis * v_const;
    let (perp1, perp2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(*axis);
    let start_rel = start - center;
    let end_rel = end - center;
    let radial0 = (start_rel.dot(&perp1).powi(2) + start_rel.dot(&perp2).powi(2)).sqrt();
    let radial1 = (end_rel.dot(&perp1).powi(2) + end_rel.dot(&perp2).powi(2)).sqrt();
    if (radial0 - *radius).abs() > radius * 0.05 || (radial1 - *radius).abs() > radius * 0.05 {
        return None;
    }
    let theta_start = start_rel.dot(&perp2).atan2(start_rel.dot(&perp1));
    let mut theta_end = end_rel.dot(&perp2).atan2(end_rel.dot(&perp1));
    while theta_end - theta_start > std::f64::consts::PI {
        theta_end -= std::f64::consts::TAU;
    }
    while theta_end - theta_start < -std::f64::consts::PI {
        theta_end += std::f64::consts::TAU;
    }
    Some((
        CurveGeometry::Ellipse {
            center,
            u_basis: perp1 * *radius,
            v_basis: perp2 * *radius,
            start_angle: theta_start,
            end_angle: theta_end,
        },
        theta_start,
        theta_end,
    ))
}

pub(super) fn attach_curve_to_edge_aligned(
    brep: &mut BRep,
    edge: EdgeHandle,
    curve: &NurbsCurve,
) -> Result<(), ApexBuildError> {
    let p = curve.degree;
    let n = curve.control_points.len() - 1;
    let t_start = curve.knots[p];
    let t_end = curve.knots[n + 1];
    let p_start = curve.point_at(t_start);
    let p_end = curve.point_at(t_end);

    let he =
        brep.edges
            .get(edge)
            .map(|e| e.half_edge)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!("attach_curve_to_edge_aligned: edge {edge:?} missing"),
            })?;
    let twin = brep.half_edges[he].twin;
    let he_origin = brep.vertices[brep.half_edges[he].origin].position;
    let twin_origin = brep.vertices[brep.half_edges[twin].origin].position;

    let forward_cost = he_origin.distance_to(&p_start) + twin_origin.distance_to(&p_end);
    let reversed_cost = he_origin.distance_to(&p_end) + twin_origin.distance_to(&p_start);
    let edge_curve = if forward_cost.min(reversed_cost) > 1.0e-5 {
        if std::env::var_os("POPKERN_A2_TRACE").is_some() {
            eprintln!(
                "[a2-trace §14.9/edge-curve] replacing mismatched planned edge curve on {edge:?}: \
                 endpoint_cost={:.3e}, edge=({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}), \
                 planned=({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                forward_cost.min(reversed_cost),
                he_origin.x,
                he_origin.y,
                he_origin.z,
                twin_origin.x,
                twin_origin.y,
                twin_origin.z,
                p_start.x,
                p_start.y,
                p_start.z,
                p_end.x,
                p_end.y,
                p_end.z
            );
        }
        CurveOnEdge::line_segment(he_origin, twin_origin, 1.0e-9).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "attach_curve_to_edge_aligned: edge {edge:?} has degenerate endpoints"
                ),
            }
        })?
    } else if forward_cost <= reversed_cost {
        canonical_curve_on_edge(curve, t_start, t_end, 1.0e-9)
    } else {
        canonical_curve_on_edge(curve, t_end, t_start, 1.0e-9)
    };

    brep.edges[edge].curve = Some(edge_curve);
    Ok(())
}

fn canonical_curve_on_edge(
    curve: &NurbsCurve,
    t_start: f64,
    t_end: f64,
    tolerance: f64,
) -> CurveOnEdge {
    let geometry = CurveGeometry::Nurbs(curve.clone());
    let geometry = geometry
        .canonical_line_in_range((t_start, t_end), tolerance)
        .unwrap_or(geometry);
    CurveOnEdge {
        geometry,
        t_start,
        t_end,
    }
}

pub(super) fn assign_fillet_strip_surface(
    brep: &mut BRep,
    strip_face: FaceHandle,
    fillet_index: usize,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<(), ApexBuildError> {
    let surface = fillet_surfaces
        .get(fillet_index)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!("missing fillet surface for FilletEnd[{fillet_index}]"),
        })?
        .surface
        .clone();
    brep.faces[strip_face].surface = Some(SurfaceOnFace {
        geometry: surface.clone(),
        orientation: true,
    });

    let Some(face) = brep.faces.get(strip_face) else {
        return Ok(());
    };
    let outer_loop = face.outer_loop;
    let Some(loop_data) = brep.loops.get(outer_loop) else {
        return Ok(());
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    let mut iters = 0usize;
    let max_iters = brep.half_edges.len() + 1;
    while brep.half_edges.contains_key(cursor) {
        let edge = brep.half_edges[cursor].edge;
        if brep
            .edges
            .get(edge)
            .and_then(|e| e.curve.as_ref())
            .is_none()
        {
            let start = brep.vertices[brep.half_edges[cursor].origin].position;
            let twin = brep.half_edges[cursor].twin;
            let end = brep.vertices[brep.half_edges[twin].origin].position;
            if let Some((geometry, t_start, t_end)) =
                cylinder_arc_through_edge_endpoints(&surface, start, end)
            {
                brep.edges[edge].curve = Some(CurveOnEdge {
                    geometry,
                    t_start,
                    t_end,
                });
            }
        }
        if let Some(curve) = brep.edges.get(edge).and_then(|e| e.curve.as_ref()) {
            let p0 = curve.geometry.point_at(curve.t_start);
            let p1 = curve.geometry.point_at(curve.t_end);
            let uv0 = crate::construct::helpers::project_to_uv(&surface, &p0);
            let uv1 = crate::construct::helpers::project_to_uv(&surface, &p1);
            brep.half_edges[cursor].pcurve =
                Some(Pcurve::linear(uv0, uv1, curve.t_start, curve.t_end));
        }
        cursor = brep.half_edges[cursor].next;
        iters += 1;
        if cursor == start || iters >= max_iters {
            break;
        }
    }
    Ok(())
}
