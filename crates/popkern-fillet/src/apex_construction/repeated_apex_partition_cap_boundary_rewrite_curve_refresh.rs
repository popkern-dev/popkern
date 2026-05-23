//! Curve refresh for cap-boundary rewrites that reauthor shared FilletEnd seams.

use super::*;

pub(super) fn refresh_reauthored_fillet_end_curves_for_plan(
    required: bool,
    brep: &mut BRep,
    topology_plans: &[CylinderTrimRewriteTopologyPlan],
) -> Result<(), ApexBuildError> {
    if !required {
        return Ok(());
    }
    for plan in topology_plans {
        let FilletEndSeamRewrite::ReauthorCapBoundaryAndCylinderSide {
            edge,
            current_cylinder_half_edge,
            ..
        } = plan.seam_rewrite
        else {
            continue;
        };
        let start = brep
            .half_edges
            .get(current_cylinder_half_edge)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot refresh reauthored FilletEnd[{}] \
                     curve: stale cap-side half-edge {current_cylinder_half_edge:?}",
                    plan.fillet_index
                ),
            })?;
        let end = half_edge_end_vertex(brep, current_cylinder_half_edge).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot refresh reauthored FilletEnd[{}] \
                     curve: half-edge {current_cylinder_half_edge:?} has stale twin",
                    plan.fillet_index
                ),
            }
        })?;
        let curve = exact_cylinder_arc_trim_curve_between(
            brep,
            plan.cylinder_face,
            start,
            end,
            "reauthored shared FilletEnd seam",
            plan.fillet_index,
        )
        .or_else(|_| {
            exact_cylinder_trim_curve_between(
                brep,
                plan.cylinder_face,
                start,
                end,
                "reauthored shared FilletEnd seam",
                plan.fillet_index,
            )
        })?;
        brep.edges[edge].curve = Some(curve);
    }
    Ok(())
}

pub(super) fn author_materialized_shared_fillet_end_trims_for_plan(
    brep: &mut BRep,
    topology_plans: &[CylinderTrimRewriteTopologyPlan],
) -> Result<(), ApexBuildError> {
    for plan in topology_plans {
        let (edge, half_edge_a, half_edge_b, label, defer_cap_side) = match plan.seam_rewrite {
            FilletEndSeamRewrite::ReuseCurrentTwinPair {
                edge,
                cap_half_edge,
                current_cylinder_half_edge,
            } => (
                edge,
                cap_half_edge,
                current_cylinder_half_edge,
                "reused shared FilletEnd seam",
                false,
            ),
            FilletEndSeamRewrite::ReauthorCapBoundaryAndCylinderSide {
                edge,
                cap_half_edge,
                current_cylinder_half_edge,
            } => (
                edge,
                current_cylinder_half_edge,
                cap_half_edge,
                "reauthored shared FilletEnd seam",
                true,
            ),
        };
        let Some(edge_half_edge) = brep.edges.get(edge).map(|edge_data| edge_data.half_edge) else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot author shared FilletEnd[{}] trim: \
                     stale edge {edge:?}",
                    plan.fillet_index
                ),
            });
        };
        let start = brep
            .half_edges
            .get(edge_half_edge)
            .map(|he| he.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot author shared FilletEnd[{}] trim: \
                     stale canonical half-edge {edge_half_edge:?}",
                    plan.fillet_index
                ),
            })?;
        let end = half_edge_end_vertex(brep, edge_half_edge).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot author shared FilletEnd[{}] trim: \
                     canonical half-edge {edge_half_edge:?} has stale twin",
                    plan.fillet_index
                ),
            }
        })?;
        let curve = exact_cylinder_arc_trim_curve_between(
            brep,
            plan.cylinder_face,
            start,
            end,
            label,
            plan.fillet_index,
        )
        .or_else(|_| {
            exact_cylinder_trim_curve_between(
                brep,
                plan.cylinder_face,
                start,
                end,
                label,
                plan.fillet_index,
            )
        })?;
        author_materialized_shared_fillet_end_trim(
            brep,
            plan,
            edge,
            half_edge_a,
            half_edge_b,
            curve,
            defer_cap_side,
        )?;
    }
    Ok(())
}

fn author_materialized_shared_fillet_end_trim(
    brep: &mut BRep,
    plan: &CylinderTrimRewriteTopologyPlan,
    edge: EdgeHandle,
    half_edge_a: HalfEdgeHandle,
    half_edge_b: HalfEdgeHandle,
    curve: CurveOnEdge,
    defer_cap_side: bool,
) -> Result<(), ApexBuildError> {
    let tolerance = brep.tolerance.linear().max(1.0e-6);
    let surface_a = if defer_cap_side {
        None
    } else {
        surface_for_half_edge(brep, half_edge_a).cloned()
    };
    let surface_b = surface_for_half_edge(brep, half_edge_b).cloned();
    match (surface_a, surface_b, defer_cap_side) {
        (Some(surface_a), Some(surface_b), false) => {
            let pcurve_a = projected_shared_trim_pcurve(brep, half_edge_a, &curve, &surface_a, 33)
                .ok_or_else(|| pcurve_projection_error(plan, half_edge_a, "paired side"))?;
            let pcurve_b = projected_shared_trim_pcurve(brep, half_edge_b, &curve, &surface_b, 33)
                .ok_or_else(|| pcurve_projection_error(plan, half_edge_b, "paired side"))?;
            commit_shared_trim_edge_transactionally(
                brep,
                SharedTrimEdgeInput {
                    edge,
                    half_edge_a,
                    half_edge_b,
                    curve,
                    pcurve_a,
                    pcurve_b,
                    tolerance,
                    validate_brep: false,
                    provenance: "repeated-apex materialized shared FilletEnd trim",
                },
            )
            .map_err(|failure| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot author shared FilletEnd[{}] trim: \
                     shared-trim validation declined {:?} with curve/surface {:.3e}, side \
                     agreement {:.3e}, endpoint {:.3e}",
                    plan.fillet_index,
                    failure.reason,
                    failure.validation.max_curve_surface_deviation,
                    failure.validation.max_side_agreement_deviation,
                    failure.validation.max_endpoint_gap
                ),
            })?;
        }
        (surface_a, surface_b, _) => {
            brep.edges[edge].curve = Some(curve.clone());
            clear_stale_reauthored_pcurve(brep, half_edge_a);
            clear_stale_reauthored_pcurve(brep, half_edge_b);
            let authored_a = author_single_surface_trim_pcurve(
                brep,
                plan,
                half_edge_a,
                &curve,
                surface_a.as_ref(),
                tolerance,
            )?;
            let authored_b = author_single_surface_trim_pcurve(
                brep,
                plan,
                half_edge_b,
                &curve,
                surface_b.as_ref(),
                tolerance,
            )?;
            if !authored_a && !authored_b {
                validate_curve_on_planned_cylinder(brep, plan, &curve, tolerance)?;
            }
        }
    }
    Ok(())
}

fn clear_stale_reauthored_pcurve(brep: &mut BRep, half_edge: HalfEdgeHandle) {
    if let Some(half_edge_data) = brep.half_edges.get_mut(half_edge) {
        half_edge_data.pcurve = None;
    }
}

fn author_single_surface_trim_pcurve(
    brep: &mut BRep,
    plan: &CylinderTrimRewriteTopologyPlan,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
    surface: Option<&SurfaceGeometry>,
    tolerance: f64,
) -> Result<bool, ApexBuildError> {
    let Some(surface) = surface else {
        return Ok(false);
    };
    let pcurve = projected_shared_trim_pcurve(brep, half_edge, curve, surface, 33)
        .ok_or_else(|| pcurve_projection_error(plan, half_edge, "single surfaced side"))?;
    let deviation = single_surface_pcurve_deviation(brep, half_edge, curve, &pcurve, surface)
        .ok_or_else(|| pcurve_projection_error(plan, half_edge, "single surfaced side"))?;
    if deviation > tolerance {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author shared FilletEnd[{}] trim: \
                 single-surface pcurve validation for half-edge {half_edge:?} exceeded tolerance \
                 with deviation {:.3e} > {:.3e}",
                plan.fillet_index, deviation, tolerance
            ),
        });
    }
    brep.half_edges[half_edge].pcurve = Some(pcurve);
    Ok(true)
}

fn surface_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<&SurfaceGeometry> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let loop_data = brep.loops.get(half_edge_data.loop_)?;
    brep.faces
        .get(loop_data.face)?
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)
}

fn validate_curve_on_planned_cylinder(
    brep: &BRep,
    plan: &CylinderTrimRewriteTopologyPlan,
    curve: &CurveOnEdge,
    tolerance: f64,
) -> Result<(), ApexBuildError> {
    let Some(surface) = brep
        .faces
        .get(plan.cylinder_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author shared FilletEnd[{}] trim: \
                 planned cylinder face {:?} has no surface",
                plan.fillet_index, plan.cylinder_face
            ),
        });
    };
    let samples = 17usize;
    let mut max_deviation: f64 = 0.0;
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let curve_t = curve.t_start + (curve.t_end - curve.t_start) * fraction;
        let point = curve.geometry.point_at(curve_t);
        let Some((u, v)) = surface.closest_point(&point) else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot author shared FilletEnd[{}] trim: \
                     planned cylinder projection failed on face {:?}",
                    plan.fillet_index, plan.cylinder_face
                ),
            });
        };
        max_deviation = max_deviation.max(point.distance_to(&surface.point_at(u, v)));
    }
    if max_deviation > tolerance {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author shared FilletEnd[{}] trim: \
                 exact curve is not on planned cylinder face {:?}, deviation {:.3e} > {:.3e}",
                plan.fillet_index, plan.cylinder_face, max_deviation, tolerance
            ),
        });
    }
    Ok(())
}

fn pcurve_projection_error(
    plan: &CylinderTrimRewriteTopologyPlan,
    half_edge: HalfEdgeHandle,
    side_kind: &'static str,
) -> ApexBuildError {
    ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition cannot author shared FilletEnd[{}] trim: \
             {side_kind} half-edge {half_edge:?} pcurve projection failed",
            plan.fillet_index
        ),
    }
}

fn projected_shared_trim_pcurve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
    surface: &SurfaceGeometry,
    sample_count: usize,
) -> Option<Pcurve> {
    let forward = half_edge_curve_forward(brep, half_edge, curve)?;
    let samples = sample_count.max(2);
    let (period_u, period_v) = surface.periodicity();
    let mut prev_uv = None;
    let mut uv_points = Vec::with_capacity(samples);
    let mut params = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let curve_fraction = if forward { fraction } else { 1.0 - fraction };
        let t = curve.t_start + (curve.t_end - curve.t_start) * curve_fraction;
        let point = curve.geometry.point_at(t);
        let (mut u, mut v) = project_shared_trim_point(surface, &point, prev_uv)?;
        unwrap_periodic_projection(&mut u, &mut v, prev_uv, period_u, period_v);
        prev_uv = Some((u, v));
        uv_points.push(Point3::new(u, v, 0.0));
        params.push(fraction);
    }
    let degree = 5.min(uv_points.len().saturating_sub(1));
    popkern_core::geometry::fitting::interpolate_curve(&uv_points, &params, degree)
        .ok()
        .map(Pcurve::from_nurbs)
}

fn half_edge_curve_forward(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
) -> Option<bool> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let origin = brep.vertices.get(half_edge_data.origin)?.position;
    let curve_start = curve.geometry.point_at(curve.t_start);
    let curve_end = curve.geometry.point_at(curve.t_end);
    Some(origin.distance_to(&curve_start) <= origin.distance_to(&curve_end))
}

fn single_surface_pcurve_deviation(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
) -> Option<f64> {
    let forward = half_edge_curve_forward(brep, half_edge, curve)?;
    let (pc_start, pc_end) = pcurve.parameter_range();
    let samples = 17usize;
    let mut max_deviation: f64 = 0.0;
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let curve_fraction = if forward { fraction } else { 1.0 - fraction };
        let curve_t = curve.t_start + (curve.t_end - curve.t_start) * curve_fraction;
        let pcurve_t = pc_start + (pc_end - pc_start) * fraction;
        let uv = pcurve.point_at(pcurve_t);
        let curve_point = curve.geometry.point_at(curve_t);
        let surface_point = surface.point_at(uv.x, uv.y);
        max_deviation = max_deviation.max(curve_point.distance_to(&surface_point));
    }
    Some(max_deviation)
}

fn project_shared_trim_point(
    surface: &SurfaceGeometry,
    point: &Point3,
    prev_uv: Option<(f64, f64)>,
) -> Option<(f64, f64)> {
    let exact = if let Some((prev_u, prev_v)) = prev_uv {
        popkern_core::geometry::uv_projection::project_point_to_uv_near(
            surface, point, prev_u, prev_v,
        )
    } else {
        popkern_core::geometry::uv_projection::project_point_to_uv(surface, point)
    };
    exact.or_else(|| surface.closest_point(point))
}

fn unwrap_periodic_projection(
    u: &mut f64,
    v: &mut f64,
    prev_uv: Option<(f64, f64)>,
    period_u: Option<f64>,
    period_v: Option<f64>,
) {
    let Some((prev_u, prev_v)) = prev_uv else {
        return;
    };
    if let Some(period) = period_u {
        while *u - prev_u > period * 0.5 {
            *u -= period;
        }
        while prev_u - *u > period * 0.5 {
            *u += period;
        }
    }
    if let Some(period) = period_v {
        while *v - prev_v > period * 0.5 {
            *v -= period;
        }
        while prev_v - *v > period * 0.5 {
            *v += period;
        }
    }
}
