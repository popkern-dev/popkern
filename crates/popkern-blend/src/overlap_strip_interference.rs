//! PK-140 mid-face strip interference detection.
//!
//! Split from `overlap.rs` (file-size cap). See
//! [`detect_strip_interference`].

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;

/// A mid-face strip interference between two vertex-disjoint fillets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StripInterference {
    /// Spec index of the first fillet.
    pub fillet_a: usize,
    /// Spec index of the second fillet.
    pub fillet_b: usize,
    /// Closest approach of the two fillet edges (sampled), the space
    /// both removal strips must share on the common support face.
    /// Drives the honest radius bound when a ridge-completed pair
    /// later fails to construct (PK-380).
    pub gap: f64,
}

/// Detect mid-face strip interference: two fillets on VERTEX-DISJOINT
/// edges sharing a support face whose removal strips overlap. Their
/// rails are typically PARALLEL (opposite edges of one face), so the
/// crossing test in [`detect_overlaps`] never fires — the overlap is
/// a region containment: one fillet's rail lies closer to the other
/// fillet's edge than that fillet's own rail does.
///
/// Corner pairs (edges sharing a vertex) are excluded — adjacent
/// fillets legitimately meet at vertex blends and their rails always
/// approach the shared corner.
pub fn detect_strip_interference(
    edge_geometries: &[crate::construct::EdgeGeometry],
    brep: &popkern_core::topology::BRep,
    tolerance: &Tolerance,
) -> Vec<StripInterference> {
    // A rail is sampled densely and each sample measured against the
    // other curve's closest point, the nearest sample polished along
    // the rail: 17 samples measured against 17 missed a boss root
    // fillet's foot circle running 0.01 into a plate edge round's band
    // (sampled 0.110 and 0.097 against true 0.07 inside the 0.08
    // width), and the strips were served overlapping with 52 open mesh
    // edges (PK-1189).
    let edge_window = |index: usize| -> Option<EdgeWindow<'_>> {
        let analysis = &edge_geometries[index].analysis;
        brep.edges
            .get(analysis.edge)
            .and_then(|edge| edge.curve.as_ref())
            .map(|curve| (&curve.geometry, curve.t_start, curve.t_end))
    };
    let mut interferences = Vec::new();
    for i in 0..edge_geometries.len() {
        for j in (i + 1)..edge_geometries.len() {
            let (a, b) = (&edge_geometries[i].analysis, &edge_geometries[j].analysis);
            // Corner pairs are legitimate vertex-blend territory.
            if a.vertex_start == b.vertex_start
                || a.vertex_start == b.vertex_end
                || a.vertex_end == b.vertex_start
                || a.vertex_end == b.vertex_end
            {
                continue;
            }
            // Shared support face pairs, with each fillet's trim on
            // that face.
            let shared: Vec<(&NurbsCurve, &NurbsCurve)> = [
                (a.face_1, &edge_geometries[i].trim_3d_1),
                (a.face_2, &edge_geometries[i].trim_3d_2),
            ]
            .into_iter()
            .flat_map(|(fa, trim_a)| {
                [
                    (b.face_1, &edge_geometries[j].trim_3d_1),
                    (b.face_2, &edge_geometries[j].trim_3d_2),
                ]
                .into_iter()
                .filter(move |(fb, _)| *fb == fa)
                .map(move |(_, trim_b)| (trim_a, trim_b))
            })
            .collect();
            if shared.is_empty() {
                continue;
            }
            let (Some(edge_i), Some(edge_j)) = (edge_window(i), edge_window(j)) else {
                continue;
            };
            let closed_i = a.vertex_start == a.vertex_end;
            let closed_j = b.vertex_start == b.vertex_end;
            let trace = popkern_core::trace::trace_enabled("interference");
            let hit = shared.iter().any(|(trim_a, trim_b)| {
                // Rails are authored with overhang beyond the edge's
                // swept extent (the padded offset-surface intersection
                // spans the whole face domain). The strip only exists
                // where its spine parameter lies within the edge span,
                // so overhang samples must not vote — an overhanging
                // rail crossing another fillet's band 35 mm past its
                // own edge's end is not an interference (measured:
                // the skewed-rib plate's junction rail extended to
                // x = ±67 on an 80-long plate and crossed both rim
                // strips, declining a part whose real strips are 7
                // apart — PK-185). A sample is in-band when its foot
                // on the edge is interior (past a thirty-second of the
                // edge from either end); a closed edge has no ends.
                if coarsely_apart(trim_a, trim_b, edge_i, edge_j) {
                    return false;
                }
                let rail_a = in_band_params(trim_a, edge_i, closed_i);
                let rail_b = in_band_params(trim_b, edge_j, closed_j);
                let width_a = rail_to_edge(trim_a, &rail_a, edge_i);
                let width_b = rail_to_edge(trim_b, &rail_b, edge_j);
                let b_into_a = rail_to_edge(trim_b, &rail_b, edge_i);
                let a_into_b = rail_to_edge(trim_a, &rail_a, edge_j);
                let slack = tolerance.linear().max(1.0e-9) * 10.0;
                let hit = b_into_a + slack < width_a || a_into_b + slack < width_b;
                if trace {
                    eprintln!(
                        "[idetect] pair ({i},{j}) hit={hit}: width_a={width_a:.4} width_b={width_b:.4} b_into_a={b_into_a:.4} a_into_b={a_into_b:.4}"
                    );
                }
                hit
            });
            if hit {
                interferences.push(StripInterference {
                    fillet_a: i,
                    fillet_b: j,
                    gap: edge_to_edge(edge_i, edge_j),
                });
            }
        }
    }
    interferences
}

type EdgeWindow<'a> = (&'a popkern_core::geometry::curve::CurveGeometry, f64, f64);

const RAIL_SAMPLES: usize = 64;

/// Samples per curve in the coarse pass.
const COARSE_SAMPLES: usize = 16;

/// Points of a curve over its window, and half the widest gap between
/// neighbours (how far its points may lie from a sample).
fn coarse_points(
    points: impl Fn(f64) -> popkern_core::geometry::Point3,
) -> (Vec<popkern_core::geometry::Point3>, f64) {
    let pts: Vec<popkern_core::geometry::Point3> = (0..=COARSE_SAMPLES)
        .map(|k| points(k as f64 / COARSE_SAMPLES as f64))
        .collect();
    let gap = pts
        .windows(2)
        .map(|w| w[0].distance_to(&w[1]))
        .fold(0.0_f64, f64::max);
    (pts, 0.5 * gap)
}

fn sampled_min(a: &[popkern_core::geometry::Point3], b: &[popkern_core::geometry::Point3]) -> f64 {
    a.iter()
        .flat_map(|p| b.iter().map(move |q| p.distance_to(q)))
        .fold(f64::INFINITY, f64::min)
}

/// Whether coarse samples prove neither rail can come closer to the other
/// strip's edge than that strip's own rail does: the sampled gap less both
/// curves' sampling reach still exceeds the sampled width (an upper bound
/// of the true width). Keeps the dense measure to pairs near enough to
/// interfere (a 12-edge cube's 66 pairs allocated 51 MB measured densely).
fn coarsely_apart(
    rail_a: &NurbsCurve,
    rail_b: &NurbsCurve,
    edge_i: EdgeWindow<'_>,
    edge_j: EdgeWindow<'_>,
) -> bool {
    let rail = |c: &NurbsCurve| {
        let (t0, t1) = c.parameter_range();
        coarse_points(|f| c.point_at(t0 + (t1 - t0) * f))
    };
    let edge = |(c, t0, t1): EdgeWindow<'_>| coarse_points(|f| c.point_at(t0 + (t1 - t0) * f));
    let ((ra, ha), (rb, hb)) = (rail(rail_a), rail(rail_b));
    let ((ei, hi), (ej, hj)) = (edge(edge_i), edge(edge_j));
    let width_a = sampled_min(&ra, &ei);
    let width_b = sampled_min(&rb, &ej);
    let b_into_a = sampled_min(&rb, &ei) - hb - hi;
    let a_into_b = sampled_min(&ra, &ej) - ha - hj;
    b_into_a > width_a && a_into_b > width_b
}

/// The foot of `p` on the edge's window: its parameter and distance.
fn edge_foot(edge: EdgeWindow<'_>, p: &popkern_core::geometry::Point3) -> (f64, f64) {
    let (curve, t0, t1) = edge;
    let t = popkern_core::geometry::curve::closest_point_in_window(curve, p, t0, t1);
    (t, curve.point_at(t).distance_to(p))
}

/// The rail's sample parameters whose feet fall inside the edge (all of
/// them on a closed edge); the whole rail when none does, so the pair
/// still measures.
fn in_band_params(rail: &NurbsCurve, edge: EdgeWindow<'_>, closed: bool) -> Vec<f64> {
    let (r0, r1) = rail.parameter_range();
    let all: Vec<f64> = (0..=RAIL_SAMPLES)
        .map(|k| r0 + (r1 - r0) * k as f64 / RAIL_SAMPLES as f64)
        .collect();
    if closed {
        return all;
    }
    let (_, t0, t1) = edge;
    let (lo, hi) = (t0.min(t1), t0.max(t1));
    let margin = (hi - lo) / 32.0;
    let kept: Vec<f64> = all
        .iter()
        .copied()
        .filter(|&t| {
            let (foot, _) = edge_foot(edge, &rail.point_at(t));
            foot > lo + margin && foot < hi - margin
        })
        .collect();
    if kept.is_empty() { all } else { kept }
}

/// The least distance from the rail's samples `params` to the edge,
/// the nearest sample polished along the rail between its neighbours.
fn rail_to_edge(rail: &NurbsCurve, params: &[f64], edge: EdgeWindow<'_>) -> f64 {
    let dist = |t: f64| edge_foot(edge, &rail.point_at(t)).1;
    let Some((k, mut best)) = params
        .iter()
        .map(|&t| dist(t))
        .enumerate()
        .min_by(|x, y| x.1.total_cmp(&y.1))
    else {
        return f64::INFINITY;
    };
    let mut lo = params[k.saturating_sub(1)];
    let mut hi = params[(k + 1).min(params.len() - 1)];
    for _ in 0..48 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if dist(m1) < dist(m2) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    best = best.min(dist(0.5 * (lo + hi)));
    best
}

/// The closest approach of two edges (sampled on the first, measured
/// to the second's closest point).
fn edge_to_edge(a: EdgeWindow<'_>, b: EdgeWindow<'_>) -> f64 {
    let (curve, t0, t1) = a;
    (0..=RAIL_SAMPLES)
        .map(|k| {
            let t = t0 + (t1 - t0) * k as f64 / RAIL_SAMPLES as f64;
            edge_foot(b, &curve.point_at(t)).1
        })
        .fold(f64::INFINITY, f64::min)
}

/// PK-140 thin-wall precheck: returns `Some((requested_radius,
/// available_extent))` when a PLANAR support face's full extent
/// perpendicular to the filleted edge is smaller than the requested
/// radius — the rail would exit the face across its far boundary.
/// Curved supports are skipped (their reach is not a vertex-extent
/// question); the max-vertex-reach extent only errs permissive.
pub(crate) fn planar_support_extent_overrun(
    brep: &popkern_core::topology::BRep,
    analysis: &crate::analysis::EdgeAnalysis,
    radius: &crate::FilletRadius,
    tolerance: &Tolerance,
) -> Option<(f64, f64)> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let requested = match radius {
        crate::FilletRadius::Constant(r) => *r,
        crate::FilletRadius::Variable { start, end } => start.max(*end),
        crate::FilletRadius::ByPoints(points) => {
            points.iter().map(|(_, r)| *r).fold(0.0_f64, f64::max)
        }
    };
    if !requested.is_finite() || requested <= 0.0 {
        return None;
    }
    let edge_a = brep.vertices.get(analysis.vertex_start)?.position;
    let edge_b = brep.vertices.get(analysis.vertex_end)?.position;
    let edge_dir = (edge_b - edge_a).normalized()?;
    let slack = tolerance.linear().max(1.0e-9) * 10.0;
    for (face, surface) in [
        (analysis.face_1, &analysis.surface_1),
        (analysis.face_2, &analysis.surface_2),
    ] {
        if !matches!(surface, SurfaceGeometry::Plane { .. }) {
            continue;
        }
        let mut extent = 0.0_f64;
        for p in brep.face_outer_loop_positions(face) {
            let d = p - edge_a;
            let perp = d - edge_dir * d.dot(&edge_dir);
            extent = extent.max(perp.length());
        }
        if extent > 0.0 && requested > extent + slack {
            return Some((requested, extent));
        }
    }
    None
}
