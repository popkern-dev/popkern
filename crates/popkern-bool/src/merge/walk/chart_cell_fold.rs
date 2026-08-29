//! Whole-period fold of chart embeddings into a sub-period face's
//! chart cell (PK-304/305 — extracted from `chart_subdivision.rs`
//! under the file-size ratchet). See the call site for the measured
//! rationale: the anchor is the run's GRAPH-NODE endpoints projected
//! to the surface and folded into the cell, not the run midpoint.

use popkern_core::geometry::surface::SurfaceGeometry;

use super::chart_subdivision::FaceChartSubdivision;
use crate::clip::chart_cell::ChartCell;
use crate::merge::graph::{FaceRef, IntersectionGraph};

/// Fold every embedding of `sub` into `cell` by whole-period shifts
/// per run, anchored to the run's graph-node endpoints.
pub(crate) fn fold_embeddings_into_cell(
    graph: &IntersectionGraph,
    surface: &SurfaceGeometry,
    sub: &mut FaceChartSubdivision,
    cell: &ChartCell,
    pu: &Option<f64>,
    pv: &Option<f64>,
    face: &FaceRef,
) {
    // Whole-period fold per RUN (a sub-period face's runs never
    // straddle a seam, so one shift per run is exact). The anchor
    // is the run's GRAPH-NODE endpoints projected to the surface
    // and folded into the cell — NOT the run midpoint: a run
    // hugging the cell border (the seam-adjacent quarter's own
    // boundary at u ≡ 0 ≡ 2π) rounds a midpoint fold to either
    // border arbitrarily, and one wrong-side run flips the
    // boundary shoelace so σ mis-calibrates and the role assembly
    // turns region-separating records into holes (measured:
    // A-face 5v1 walked σ = −1, its whole partition served as ONE
    // Inside mega-fragment, and the union discarded the face).
    let fold_to_cell = |x: f64, range: (f64, f64), p: &Option<f64>| match p {
        Some(p) if p.is_finite() && *p > 0.0 => {
            x + ((0.5 * (range.0 + range.1) - x) / p).round() * p
        }
        _ => x,
    };
    let node_uv = |n: usize| -> Option<(f64, f64)> {
        let pos = graph.nodes[n].position;
        let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(surface, &pos)?;
        Some((
            fold_to_cell(u, cell.u_range, pu),
            fold_to_cell(v, cell.v_range, pv),
        ))
    };
    let keys: Vec<usize> = sub.embeddings.keys().copied().collect();
    for h in keys {
        let (n0, n1) = (
            graph.half_edges[h].origin,
            graph.half_edges[graph.half_edges[h].twin].origin,
        );
        let emb = sub.embeddings.get_mut(&h).unwrap_or_else(|| unreachable!());
        let start = emb.uv[0];
        let end = emb.uv[emb.uv.len() - 1];
        let anchor = |p: &Option<f64>,
                      run_a: f64,
                      run_b: f64,
                      node_a: Option<f64>,
                      node_b: Option<f64>,
                      range: (f64, f64)|
         -> f64 {
            let Some(period) = p else { return 0.0 };
            match (node_a, node_b) {
                (Some(na), Some(nb)) => {
                    // Whole-period shift aligning run endpoints to
                    // the node anchors.
                    (((na - run_a) + (nb - run_b)) / (2.0 * period)).round() * period
                }
                _ => {
                    let mid = 0.5 * (run_a + run_b);
                    ((0.5 * (range.0 + range.1) - mid) / period).round() * period
                }
            }
        };
        let (a0, a1) = (node_uv(n0), node_uv(n1));
        let du = anchor(
            pu,
            start.0,
            end.0,
            a0.map(|q| q.0),
            a1.map(|q| q.0),
            cell.u_range,
        );
        let dv = anchor(
            pv,
            start.1,
            end.1,
            a0.map(|q| q.1),
            a1.map(|q| q.1),
            cell.v_range,
        );
        if du != 0.0 || dv != 0.0 {
            for q in emb.uv.iter_mut() {
                q.0 += du;
                q.1 += dv;
            }
        }
        if popkern_core::trace::trace_enabled("fold-scan") {
            let e2 = &sub.embeddings[&h];
            eprintln!(
                "      [fold] {face:?} he {h} src {:?} ({:.3},{:.3})->({:.3},{:.3}) shift ({du:.3},{dv:.3})",
                graph.half_edges[h].source,
                e2.uv[0].0,
                e2.uv[0].1,
                e2.uv[e2.uv.len() - 1].0,
                e2.uv[e2.uv.len() - 1].1
            );
        }
    }
}
