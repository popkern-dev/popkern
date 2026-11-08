//! Whether a residual fragment fits its plane carrier: every side on
//! the carrier and every vertex already in production topology.

use super::*;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

fn point_plane_distance(carrier: &SurfaceGeometry, point: popkern_core::geometry::Point3) -> f64 {
    let Some((u, v)) = carrier.closest_point(&point) else {
        return f64::INFINITY;
    };
    carrier.point_at(u, v).distance_to(&point)
}

/// Every fragment side fits the plane carrier and every fragment vertex
/// already participates in production topology.
pub(super) fn fragment_absorbability(
    brep: &BRep,
    fragment: FaceHandle,
    carrier: &SurfaceGeometry,
    production_vertices: &std::collections::BTreeSet<VertexHandle>,
) -> Result<(), String> {
    let fragment_data = brep
        .faces
        .get(fragment)
        .ok_or_else(|| "stale fragment".to_owned())?;
    let cycle = collect_loop_half_edges(brep, fragment_data.outer_loop)
        .map_err(|_| "fragment loop does not close".to_owned())?;
    let fit_tolerance = brep.tolerance.linear().max(1.0e-9) * 10.0;
    for side in cycle {
        let data = brep
            .half_edges
            .get(side)
            .ok_or_else(|| format!("stale side {side:?}"))?;
        if !production_vertices.contains(&data.origin) {
            return Err(format!(
                "vertex {:?} is scaffold-interior (not on any production loop)",
                data.origin
            ));
        }
        let origin = brep
            .vertices
            .get(data.origin)
            .ok_or_else(|| format!("stale vertex on {side:?}"))?;
        let vertex_gap = point_plane_distance(carrier, origin.position);
        if vertex_gap > fit_tolerance {
            return Err(format!(
                "vertex {:?} off the carrier by {vertex_gap:.3e}",
                data.origin
            ));
        }
        if let Some(curve) = brep
            .edges
            .get(data.edge)
            .and_then(|edge| edge.curve.as_ref())
        {
            let midpoint = curve.geometry.point_at(0.5 * (curve.t_start + curve.t_end));
            let midpoint_gap = point_plane_distance(carrier, midpoint);
            if midpoint_gap > fit_tolerance {
                return Err(format!(
                    "edge {:?} midpoint off the carrier by {midpoint_gap:.3e}",
                    data.edge
                ));
            }
        }
    }
    Ok(())
}
