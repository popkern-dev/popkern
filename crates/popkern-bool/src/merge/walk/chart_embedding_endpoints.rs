//! A chart run's ends read from the nodes they meet (split from
//! `chart_embedding.rs` under the 500-line rule).

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::{
    project_point_to_uv_multisheet, project_point_to_uv_near,
};

use crate::merge::graph::IntersectionGraph;

/// A chart point.
type Uv = (f64, f64);
/// A run's samples and its endpoint-refined samples.
type Run = (Vec<Uv>, Option<Uv>, Option<Uv>);

/// The two-point run of half-edge `h` with no curve geometry at all:
/// its node positions projected (the second on the branch nearest the
/// first, the same nearest-branch snap as the projection source; the
/// `_near` helper ignores hints on analytic surfaces, so the first is
/// projected sheet-correct).
pub(crate) fn endpoint_projection(
    graph: &IntersectionGraph,
    h: usize,
    surface: &SurfaceGeometry,
) -> Option<Run> {
    let a = graph.nodes[graph.half_edges[h].origin].position;
    let b = graph.nodes[graph.half_edges[graph.half_edges[h].twin].origin].position;
    let ua = project_point_to_uv_multisheet(surface, &a)?;
    let (mut ub, mut vb) = project_point_to_uv_near(surface, &b, ua.0, ua.1)
        .or_else(|| project_point_to_uv_multisheet(surface, &b))?;
    let (pp_u, pp_v) = surface.parameter_periodicity();
    if let Some(pp) = pp_u {
        ub += ((ua.0 - ub) / pp).round() * pp;
    }
    if let Some(pv) = pp_v {
        vb += ((ua.1 - vb) / pv).round() * pv;
    }
    Some((vec![ua, (ub, vb)], None, None))
}

/// Each end of `uv` (a run between nodes at `ends`) given its node's
/// own azimuth, the projection of the node onto `surface` on the period
/// image nearest the end, where the end maps within `band` of its node;
/// then `snap_pole_end_azimuths` for ends on a pole. Near a chart
/// pole a fitted record's u at its end carries its 3D noise over the
/// distance from the axis: a rod 3e-5 beside a ball's pole met the
/// ball's seam at a node whose two record ends read u = -2.4e-5 and
/// beyond, the seam's border window (a millionth of the period)
/// refused the seam, and the face fell to the 3D walk, declining
/// (PK-1176). The node is the meeting's authority: every run there
/// takes its azimuth.
pub(crate) fn snap_end_azimuths(
    surface: &SurfaceGeometry,
    uv: &mut [(f64, f64)],
    ends: (Point3, Point3),
    refined: (Option<Uv>, Option<Uv>),
) {
    let n = uv.len();
    if let (Some(period), true) = (surface.parameter_periodicity().0, n >= 2) {
        for (end, node) in [(0, ends.0), (n - 1, ends.1)] {
            let (u, v) = uv[end];
            if super::chart_poles::row_collapses(surface, u, v) {
                continue;
            }
            let band = 1e-4 * surface.partial_u(u, v).length().max(1e-12) * period;
            if surface.point_at(u, v).distance_to(&node) > band.max(1e-9) {
                continue;
            }
            if let Some((nu, _)) = project_point_to_uv_near(surface, &node, u, v) {
                uv[end].0 = nu + ((u - nu) / period).round() * period;
            }
        }
    }
    super::chart_poles::snap_pole_end_azimuths(surface, uv, refined.0, refined.1);
}
