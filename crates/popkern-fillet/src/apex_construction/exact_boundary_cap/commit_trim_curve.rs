//! Canonical shared-trim curve selection and pcurve-trace deviation
//! measures for the exact-boundary commit.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{CurveEval, CurveOnEdge};
use popkern_core::healing::deviation::edge_surface_deviation;
use popkern_core::topology::BRep;

use super::{
    SharedTrimCurveSelection, curve_from_pcurve_in_edge_direction, edge_endpoint_gap,
    orient_curve_to_edge_vertices,
};

pub(crate) fn canonical_shared_trim_curve(
    brep: &BRep,
    selection: SharedTrimCurveSelection<'_>,
) -> Option<CurveOnEdge> {
    let SharedTrimCurveSelection {
        edge,
        cap_half_edge,
        cap_pcurve,
        cap_surface,
        support_half_edge,
        support_pcurve,
        support_surface,
    } = selection;
    let cap_curve =
        curve_from_pcurve_in_edge_direction(brep, edge, cap_half_edge, cap_pcurve, cap_surface)
            .and_then(|curve| orient_curve_to_edge_vertices(brep, edge, curve));
    let support_curve = curve_from_pcurve_in_edge_direction(
        brep,
        edge,
        support_half_edge,
        support_pcurve,
        support_surface,
    )
    .and_then(|curve| orient_curve_to_edge_vertices(brep, edge, curve));
    match (cap_curve, support_curve) {
        (Some(cap_curve), Some(support_curve)) => {
            let cap_endpoint_gap = edge_endpoint_gap(brep, edge, &cap_curve)?;
            let support_endpoint_gap = edge_endpoint_gap(brep, edge, &support_curve)?;
            // The cap curve traces the committed cap boundary row. When it
            // already lies on the support surface within tolerance, it is
            // the seam authority: preferring the support pcurve's trace
            // here replaced exact FilletEnd arcs with a fit of a ~1e-5
            // stale cylinder pcurve, and that displacement became the
            // marginal cap/cylinder tessellation seams (PK-124).
            let cap_on_support = curve_surface_deviation(&cap_curve, support_surface);
            if cap_endpoint_gap <= support_endpoint_gap.max(1.0e-6) && cap_on_support <= 1.0e-6 {
                return Some(cap_curve);
            }
            if matches!(support_surface, SurfaceGeometry::Cylinder { .. })
                && support_endpoint_gap <= cap_endpoint_gap.max(1.0e-6)
            {
                return Some(support_curve);
            }
            if cap_endpoint_gap <= support_endpoint_gap {
                Some(cap_curve)
            } else {
                Some(support_curve)
            }
        }
        (Some(curve), None) | (None, Some(curve)) => Some(curve),
        (None, None) => None,
    }
}

/// Worst sampled distance from an edge curve to a surface.
fn curve_surface_deviation(curve: &CurveOnEdge, surface: &SurfaceGeometry) -> f64 {
    edge_surface_deviation(&curve.geometry, (curve.t_start, curve.t_end), surface, 17).max_deviation
}

/// Worst sampled distance between a pcurve's surface trace and an edge
/// curve, matched by the better of forward/reverse fraction alignment.
pub(crate) fn pcurve_trace_curve_deviation(
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    curve: &CurveOnEdge,
) -> f64 {
    let (pc_start, pc_end) = pcurve.parameter_range();
    let samples = 17usize;
    let mut forward = 0.0_f64;
    let mut reverse = 0.0_f64;
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let uv = pcurve.point_at(pc_start + (pc_end - pc_start) * fraction);
        let trace = surface.point_at(uv.x, uv.y);
        let forward_t = curve.t_start + (curve.t_end - curve.t_start) * fraction;
        let reverse_t = curve.t_end + (curve.t_start - curve.t_end) * fraction;
        forward = forward.max(trace.distance_to(&curve.geometry.point_at(forward_t)));
        reverse = reverse.max(trace.distance_to(&curve.geometry.point_at(reverse_t)));
    }
    forward.min(reverse)
}
