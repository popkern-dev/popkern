//! Edge-curve reauthoring for KP20 cap/support boundary solves.

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::fitting::fit_curve;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::topology::{BRep, HalfEdgeHandle};

pub(super) fn reauthor_from_samples(
    brep: &mut BRep,
    cap_half_edge: HalfEdgeHandle,
    samples: &[Point3],
) -> Result<(), CoreError> {
    let Some(half_edge) = brep.half_edges.get(cap_half_edge) else {
        return Err(CoreError::InvalidHandle(
            "cap/support boundary solve half-edge invalid".into(),
        ));
    };
    let edge = half_edge.edge;
    let mut edge_samples = samples.to_vec();
    if brep.edges[edge].half_edge != cap_half_edge {
        edge_samples.reverse();
    }
    let curve = fit_nurbs_curve(&edge_samples)
        .map(|geometry| CurveOnEdge {
            geometry: CurveGeometry::Nurbs(geometry),
            t_start: 0.0,
            t_end: 1.0,
        })
        .or_else(|| line_curve_from_endpoints(brep, &edge_samples))
        .ok_or_else(|| {
            CoreError::PreconditionFailed(
                "cap/support boundary solve found degenerate edge samples".into(),
            )
        })?;
    brep.edges[edge].curve = Some(curve);
    Ok(())
}

fn fit_nurbs_curve(points: &[Point3]) -> Option<NurbsCurve> {
    if points.len() < 10 {
        return None;
    }
    let params = (0..points.len())
        .map(|index| index as f64 / (points.len() - 1) as f64)
        .collect::<Vec<_>>();
    fit_curve(points, &params, 9, 10).ok()
}

fn line_curve_from_endpoints(brep: &BRep, points: &[Point3]) -> Option<CurveOnEdge> {
    let start = *points.first()?;
    let end = *points.last()?;
    CurveOnEdge::line_segment(start, end, brep.tolerance.linear().max(1.0e-12))
}
