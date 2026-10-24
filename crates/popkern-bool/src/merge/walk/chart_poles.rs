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

/// Whether the chart can carry a collapsed row (see `detect_pole_nodes`).
fn pole_shaped(surface: &SurfaceGeometry) -> bool {
    matches!(
        surface,
        SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Cone { .. } | SurfaceGeometry::Nurbs(_)
    )
}

/// Whether the chart's u-row through `(u, v)` collapses to one point
/// (a sphere pole, cone apex or a revolved profile on its axis).
pub(crate) fn row_collapses(surface: &SurfaceGeometry, u: f64, v: f64) -> bool {
    let p0 = surface.point_at(u, v);
    let p1 = surface.point_at(u + 1e-4, v);
    let scale = 1.0 + p0.distance_to(&popkern_core::geometry::Point3::new(0.0, 0.0, 0.0));
    p0.distance_to(&p1) < 1e-11 * scale
}

/// Give a run's pole endpoints their arrival azimuth (PK-877).
///
/// A pole sample's u is arbitrary (the whole row is one point), and
/// the pcurve or projection source fills it with whatever its pole
/// branch returns: the drilled dome's seam piece from the pole to the
/// rim ran at u = 0 but sampled u = −π at the pole, the seam-border
/// alignment read that −π as the seam's azimuth, declared a misfit,
/// and the sphere walked by 3D-Newell, which wired the closed rim
/// backwards (the hole wall and the sphere zone served inside-out).
/// The successor rule at a pole orders departures by this endpoint u
/// as each run's own meridian azimuth, so the endpoint takes the
/// azimuth the run approaches it along: the endpoint-refined
/// δ-sample's u when it is off the pole, else the neighbouring
/// sample's. A run lying wholly on the pole row keeps its samples.
///
/// The δ-sample is projected on its own and can alias a period off
/// the run (PK-691): the upright drilled dome's seam piece sampled
/// u = 2π along its meridian while its δ-sample read u = 0, so the
/// pole end took 0, the run swept the whole chart diagonally, both
/// seam halves sat on the u = 2π border copy, and the seam-border
/// rebalance balanced the rim node by moving the closed rim ring's
/// start a period, against its own embedding. The δ azimuth is folded
/// onto the neighbouring off-pole sample's period branch.
pub(crate) fn snap_pole_end_azimuths(
    surface: &SurfaceGeometry,
    uv: &mut [(f64, f64)],
    refined_first: Option<(f64, f64)>,
    refined_last: Option<(f64, f64)>,
) {
    let n = uv.len();
    if n < 2 || !pole_shaped(surface) {
        return;
    }
    let off_pole = |q: &(f64, f64)| !row_collapses(surface, q.0, q.1);
    let period = surface.parameter_periodicity().0;
    for (end, inner, refined) in [(0, 1, refined_first), (n - 1, n - 2, refined_last)] {
        if off_pole(&uv[end]) {
            continue;
        }
        let neighbour = Some(uv[inner]).filter(off_pole);
        let toward = match (refined.filter(off_pole), neighbour, period) {
            (Some((u, v)), Some((u_in, _)), Some(p)) if p > 0.0 => {
                Some((u + ((u_in - u) / p).round() * p, v))
            }
            (Some(r), _, _) => Some(r),
            (None, n, _) => n,
        };
        if let Some((u, _)) = toward {
            uv[end].0 = u;
        }
    }
}

/// Record pole nodes (with rotation side) on `sub` for sphere/cone
/// charts.
pub(crate) fn detect_pole_nodes(
    graph: &IntersectionGraph,
    surface: &SurfaceGeometry,
    sub: &mut FaceChartSubdivision,
    u_period: Option<f64>,
) {
    if u_period.is_some() && pole_shaped(surface) {
        // NURBS included since PK-456: a u-closed surface of
        // revolution whose profile touches the axis collapses a whole
        // v-border row to one point — the same pole shape as the
        // analytic sphere/cone, and the du-collapse test below is the
        // authority either way (the enum match is only a cheap
        // filter). Without it, a pole-ended seam clusters into u
        // branches and the face declines on degree imbalance
        // (measured: every x/y-revolved sphere-dimple boolean).
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
                if popkern_core::trace::trace_enabled("pole") {
                    let p0 = surface.point_at(u, v);
                    eprintln!(
                        "  [pole?] node {node} uv=({u:.4},{v:.4}) p0=({:.3},{:.3},{:.3}) du-dist={:.3e}",
                        p0.x,
                        p0.y,
                        p0.z,
                        p0.distance_to(&surface.point_at(u + 1e-4, v))
                    );
                }
                if row_collapses(surface, u, v) {
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

#[cfg(test)]
mod tests {
    use super::snap_pole_end_azimuths;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::geometry::{Point3, Vec3};
    use std::f64::consts::TAU;

    fn sphere() -> SurfaceGeometry {
        SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
            axis: Vec3::new(0.0, 0.0, 1.0),
        }
    }

    #[test]
    fn a_pole_end_takes_its_runs_azimuth_on_the_runs_period_branch() {
        // The upright drilled dome's seam piece (PK-691): samples on
        // the u = 2π meridian, the δ-sample projected to u = 0.
        let mut uv = vec![(0.0, 0.0), (TAU, 0.018), (TAU, 0.036)];
        snap_pole_end_azimuths(&sphere(), &mut uv, Some((0.0, 0.001)), None);
        assert_eq!(uv[0].0, TAU);
        // The arrival end too, and a δ-sample already on the branch
        // keeps its own azimuth.
        let mut uv = vec![(1.0, 0.036), (1.0, 0.018), (-3.0, 0.0)];
        snap_pole_end_azimuths(&sphere(), &mut uv, None, Some((1.01, 0.001)));
        assert_eq!(uv[2].0, 1.01);
    }

    #[test]
    fn a_pole_end_without_a_delta_sample_takes_its_neighbours_azimuth() {
        let mut uv = vec![(0.0, 0.0), (TAU, 0.018), (TAU, 0.036)];
        snap_pole_end_azimuths(&sphere(), &mut uv, None, None);
        assert_eq!(uv[0].0, TAU);
    }
}
