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
        let (edge, half_edge_a, half_edge_b, label) = match plan.seam_rewrite {
            FilletEndSeamRewrite::ReuseCurrentTwinPair {
                edge,
                cap_half_edge,
                current_cylinder_half_edge,
            } => (
                edge,
                cap_half_edge,
                current_cylinder_half_edge,
                "reused shared FilletEnd seam",
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
        let mut pcurve_source = brep.clone();
        pcurve_source.edges[edge].curve = Some(curve.clone());
        let pcurve_a =
            projected_shared_trim_pcurve(&pcurve_source, half_edge_a, 33).ok_or_else(|| {
                ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition cannot author shared FilletEnd[{}] trim: \
                         half-edge {half_edge_a:?} pcurve projection failed",
                        plan.fillet_index
                    ),
                }
            })?;
        let pcurve_b =
            projected_shared_trim_pcurve(&pcurve_source, half_edge_b, 33).ok_or_else(|| {
                ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition cannot author shared FilletEnd[{}] trim: \
                         half-edge {half_edge_b:?} pcurve projection failed",
                        plan.fillet_index
                    ),
                }
            })?;
        commit_shared_trim_edge_transactionally(
            brep,
            SharedTrimEdgeInput {
                edge,
                half_edge_a,
                half_edge_b,
                curve,
                pcurve_a,
                pcurve_b,
                tolerance: brep.tolerance.linear().max(1.0e-6),
                validate_brep: false,
                provenance: "repeated-apex materialized shared FilletEnd trim",
            },
        )
        .map_err(|failure| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author shared FilletEnd[{}] trim: \
                 shared-trim validation declined {:?} with curve/surface {:.3e}, side agreement \
                 {:.3e}, endpoint {:.3e}",
                plan.fillet_index,
                failure.reason,
                failure.validation.max_curve_surface_deviation,
                failure.validation.max_side_agreement_deviation,
                failure.validation.max_endpoint_gap
            ),
        })?;
    }
    Ok(())
}

fn projected_shared_trim_pcurve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    sample_count: usize,
) -> Option<Pcurve> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let surface = &brep
        .faces
        .get(brep.loops.get(half_edge_data.loop_)?.face)?
        .surface
        .as_ref()?
        .geometry;
    let curve = brep.edges.get(half_edge_data.edge)?.curve.as_ref()?;
    let forward = half_edge_curve_forward(brep, half_edge)?;
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
        let (mut u, mut v) = if let Some((prev_u, prev_v)) = prev_uv {
            popkern_core::geometry::uv_projection::project_point_to_uv_near(
                surface, &point, prev_u, prev_v,
            )?
        } else {
            popkern_core::geometry::uv_projection::project_point_to_uv(surface, &point)?
        };
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

fn half_edge_curve_forward(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<bool> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let curve = brep.edges.get(half_edge_data.edge)?.curve.as_ref()?;
    let origin = brep.vertices.get(half_edge_data.origin)?.position;
    let curve_start = curve.geometry.point_at(curve.t_start);
    let curve_end = curve.geometry.point_at(curve.t_end);
    Some(origin.distance_to(&curve_start) <= origin.distance_to(&curve_end))
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
