//! Pole-node detection for chart subdivision (extracted from
//! `chart_subdivision.rs` under the file-size ratchet). A
//! chart-degenerate endpoint (sphere pole, cone apex) is ONE
//! geometric point whose UV image is the whole u-line at the pole v —
//! every incident half-edge lands there at ITS OWN meridian azimuth.
//! The subdivision stays valid: the pole is a single branch point,
//! and the rotational successor there orders departures by azimuth
//! (the endpoint u), not by local UV direction. The rotation sign
//! folds σ with the pole side (rotation sense reverses between the
//! two poles).

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::chart_subdivision::FaceChartSubdivision;
use crate::merge::graph::IntersectionGraph;

/// Record pole nodes (with rotation side) on `sub` for sphere/cone
/// charts.
pub(crate) fn detect_pole_nodes(
    graph: &IntersectionGraph,
    surface: &SurfaceGeometry,
    sub: &mut FaceChartSubdivision,
    u_period: Option<f64>,
) {
    if u_period.is_some()
        && matches!(
            surface,
            SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Cone { .. }
        )
    {
        let (mut v_lo, mut v_hi) = (f64::MAX, f64::MIN);
        for emb in sub.embeddings.values() {
            for &(_, v) in &emb.uv {
                v_lo = v_lo.min(v);
                v_hi = v_hi.max(v);
            }
        }
        let v_mid = 0.5 * (v_lo + v_hi);
        let hes: Vec<usize> = sub.embeddings.keys().copied().collect();
        for h in hes {
            let emb = &sub.embeddings[&h];
            let (n_start, n_end) = (
                graph.half_edges[h].origin,
                graph.half_edges[graph.half_edges[h].twin].origin,
            );
            for (node, &(u, v)) in [(n_start, &emb.uv[0]), (n_end, &emb.uv[emb.uv.len() - 1])] {
                if sub.poles.contains_key(&node) {
                    continue;
                }
                let p0 = surface.point_at(u, v);
                let p1 = surface.point_at(u + 1e-4, v);
                let scale =
                    1.0 + p0.distance_to(&popkern_core::geometry::Point3::new(0.0, 0.0, 0.0));
                if popkern_core::trace::trace_enabled("pole") {
                    eprintln!(
                        "  [pole?] node {node} uv=({u:.4},{v:.4}) p0=({:.3},{:.3},{:.3}) du-dist={:.3e}",
                        p0.x,
                        p0.y,
                        p0.z,
                        p0.distance_to(&p1)
                    );
                }
                if p0.distance_to(&p1) < 1e-11 * scale {
                    // Pole side decides the rotation sense; σ is not
                    // yet calibrated here, so store the side and fold
                    // σ in at the end of the build.
                    let side = if v >= v_mid { 1.0 } else { -1.0 };
                    sub.poles.insert(node, side);
                }
            }
        }
    }
}
