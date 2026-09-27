//! Direction authority for CLOSED record-ring chart embeddings.
//!
//! Split from `chart_embedding` (file-size cap). The embedding builder
//! orients every run by its graph nodes — whichever run end sits at
//! the half-edge's origin node starts the run. A closed intersection
//! ring (plane × cylinder at a boss root, a bore rim) starts and ends
//! at ONE node, so that test is blind and the run keeps whatever
//! direction its source pcurve happened to have. Record pairs store
//! `curve_uv2` on the twin in the RECORD's direction (not the twin's
//! traversal), and when the pcurve's parameter window does not match
//! the 3D curve's (a closed record rotated to start at an edge-face
//! hit) the full-range pcurve is sampled forward: both sides of the
//! ring then embed the SAME 3D direction. The kept fragments on the
//! two faces claim one half-edge between them, the assembler rewires
//! the later face's loop reversed, and that face ships inside-out
//! (measured: a boss unioned onto a plate came back with its wall and
//! top cap facing into the solid; every flag-trusting consumer — the
//! rolling ball for the root fillet, which then carved a groove into
//! the boss — inherited the inversion).
//!
//! The graph's standing convention supplies the missing authority:
//! the half-edge of a pair that holds the 3D curve traverses it in
//! increasing parameter, its twin in decreasing parameter.

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use crate::merge::graph::{EdgeSource, IntersectionGraph};

/// Coarse parameter samples used to locate a run point on the curve.
const CURVE_SAMPLES: usize = 64;

/// Whether the closed run `uv` (in `raw`'s embedding, on `surface`)
/// walks against the direction `raw` must traverse its record curve.
///
/// Applies only to intersection-record pairs (`EdgeSource::Shared`)
/// whose ends meet at one node; everything else returns `false` and
/// keeps the node-endpoint orientation. Decides by a vote of the run's
/// 3D chord directions against the curve tangent at the nearest curve
/// sample; returns `false` unless the decisive votes disagree by a
/// clear majority (more than three to one, at least two), so
/// near-degenerate or chart-folded runs keep their existing direction.
pub(super) fn closed_run_against_curve(
    graph: &IntersectionGraph,
    raw: usize,
    curve: &CurveOnEdge,
    surface: &SurfaceGeometry,
    uv: &[(f64, f64)],
) -> bool {
    let twin = graph.half_edges[raw].twin;
    if graph.half_edges[raw].origin != graph.half_edges[twin].origin
        || !matches!(graph.half_edges[raw].source, EdgeSource::Shared { .. })
        || uv.len() < 3
    {
        return false;
    }
    let raw_holds_curve = graph.half_edges[raw].curve.is_some();
    let (t0, t1) = (curve.t_start, curve.t_end);
    let samples: Vec<(f64, popkern_core::geometry::Point3)> = (0..=CURVE_SAMPLES)
        .map(|k| {
            let t = t0 + (t1 - t0) * k as f64 / CURVE_SAMPLES as f64;
            (t, curve.geometry.point_at(t))
        })
        .collect();
    // Sample spacing sets the chord scale below which a vote is noise.
    let spacing = samples
        .windows(2)
        .map(|w| w[0].1.distance_to(&w[1].1))
        .fold(0.0_f64, f64::max);
    if spacing <= 0.0 {
        return false;
    }
    let (mut with, mut against) = (0usize, 0usize);
    for w in uv.windows(2) {
        let a = surface.point_at(w[0].0, w[0].1);
        let b = surface.point_at(w[1].0, w[1].1);
        let chord = b - a;
        let len = chord.length();
        if len <= spacing * 1.0e-3 {
            continue;
        }
        let mid = a + chord * 0.5;
        let Some(&(t, near)) = samples
            .iter()
            .min_by(|p, q| p.1.distance_to(&mid).total_cmp(&q.1.distance_to(&mid)))
        else {
            continue;
        };
        // A chord across a chart fold does not follow the curve.
        if near.distance_to(&mid) > len + spacing {
            continue;
        }
        let Some(tangent) = curve.geometry.derivative_at(t, 1).normalized() else {
            continue;
        };
        // Decisive only when the chord runs clearly along the tangent.
        let cos = chord.dot(&tangent) / len;
        if cos > 0.5 {
            with += 1;
        } else if cos < -0.5 {
            against += 1;
        }
    }
    // `raw` must walk WITH increasing t when it holds the curve.
    let (agree, disagree) = if raw_holds_curve {
        (with, against)
    } else {
        (against, with)
    };
    disagree >= 2 && disagree > 3 * agree
}
