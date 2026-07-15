//! Arc/boundary crossing detection for the setback vertex surgery,
//! in the strip face's own UV domain.
//!
//! The legacy detector took, per boundary edge, the minimum 3-D
//! distance between the setback arc and the edge curve and accepted
//! approaches under a fixed 5e-3. That predicate cannot tell a
//! transversal trim crossing from a near-parallel approach: on
//! hexagon6_variable the arc hugs a pre-carve split row at 3e-3
//! along its whole length and grazes two mitre chords near one end,
//! yielding 1 or 3 "crossings" where the topology has exactly 2.
//!
//! A trim crossing is a UV-domain event, so this detector projects
//! both the arc and every boundary edge into the face surface's UV
//! chart and intersects the 2-D polylines. Distances are quality
//! margins, not gates. Periodic charts (cylinder `u` in `[0, 2π)`)
//! are unwrapped by continuity and mean-shifted onto the arc's
//! branch before intersection.

use super::*;

const ARC_SAMPLES: usize = 128;
const EDGE_SAMPLES: usize = 64;

/// Find where the arc crosses the face's trim boundary. Returns the
/// two crossings sorted by arc parameter, each carrying the EDGE's
/// own 3-D point at the crossing (the splitter requires exact point
/// inversion on the edge).
pub(super) fn arc_boundary_crossings(
    brep: &BRep,
    face: FaceHandle,
    cubic: &[Point3; 4],
) -> Result<[(f64, EdgeHandle, Point3); 2], String> {
    let surface = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|sof| sof.geometry.clone())
        .ok_or_else(|| format!("crossing face {face:?} has no surface"))?;
    let u_period = match &surface {
        popkern_core::geometry::SurfaceGeometry::Cylinder { .. }
        | popkern_core::geometry::SurfaceGeometry::Torus { .. }
        | popkern_core::geometry::SurfaceGeometry::Sphere { .. } => Some(std::f64::consts::TAU),
        _ => None,
    };

    let arc = bezier_cubic_curve(cubic);
    let arc_eval: &dyn CurveEval = &popkern_core::geometry::CurveGeometry::Nurbs(arc);
    let arc_uv = project_polyline(&surface, u_period, |s| arc_eval.point_at(s), ARC_SAMPLES)
        .ok_or_else(|| "setback arc does not project onto the strip chart".to_string())?;
    let arc_u_mean = mean_u(&arc_uv);

    let trace = std::env::var_os("PK_CARVE_TRACE").is_some();
    let mut hits: Vec<(f64, EdgeHandle, Point3, f64)> = Vec::new();
    for he in face_loop_half_edges(brep, face)? {
        let edge = brep.half_edges[he].edge;
        let Some(curve) = brep.edges[edge].curve.as_ref() else {
            continue;
        };
        let (e0, e1) = (curve.t_start, curve.t_end);
        if !(e0.is_finite() && e1.is_finite()) {
            continue;
        }
        let edge_eval: &dyn CurveEval = &curve.geometry;
        let Some(mut edge_uv) = project_polyline(
            &surface,
            u_period,
            |s| edge_eval.point_at(e0 + (e1 - e0) * s),
            EDGE_SAMPLES,
        ) else {
            continue;
        };
        // Bring the edge polyline onto the arc's period branch.
        if let Some(period) = u_period {
            let shift = ((arc_u_mean - mean_u(&edge_uv)) / period).round() * period;
            for p in &mut edge_uv {
                p.0 += shift;
            }
        }
        for (arc_s, edge_s) in polyline_intersections(&arc_uv, &edge_uv) {
            let edge_t = e0 + (e1 - e0) * edge_s;
            let point = edge_eval.point_at(edge_t);
            let margin = arc_eval.point_at(arc_s).distance_to(&point);
            if trace {
                eprintln!(
                    "    [crossings] edge {edge:?}: UV crossing at arc t {arc_s:.3} \
                     edge s {edge_s:.3} margin {margin:.3e}"
                );
            }
            hits.push((arc_s, edge, point, margin));
        }
    }
    if trace {
        let a0 = arc_eval.point_at(0.0);
        let a1 = arc_eval.point_at(1.0);
        eprintln!(
            "    [crossings] arc ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}): {} UV hits",
            a0.x,
            a0.y,
            a0.z,
            a1.x,
            a1.y,
            a1.z,
            hits.len()
        );
    }
    // Dedup near-identical arc parameters (a crossing at a shared
    // vertex between consecutive edges registers on both).
    hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    hits.dedup_by(|a, b| (a.0 - b.0).abs() < 1.0e-3);
    // Exact-corner endpoint contact: with uniform radii the trim
    // intersections are exact (measured 4e-8 corner gaps on
    // pentagon5_uniform), so the arc ENDS ON the strip boundary
    // instead of crossing it — there is no sign flip at that end for
    // the polyline test to find. Synthesize the missing crossing at
    // an arc endpoint that lies on a boundary edge within a band.
    if hits.len() < 2 {
        let chord = arc_eval.point_at(0.0).distance_to(&arc_eval.point_at(1.0));
        let band = (1.0e-3 * chord).max(1.0e-9);
        for t_missing in [0.0, 1.0] {
            if hits.iter().any(|h| (h.0 - t_missing).abs() < 1.0e-3) {
                continue;
            }
            let p_end = arc_eval.point_at(t_missing);
            let mut best: Option<(f64, EdgeHandle, Point3)> = None;
            for he in face_loop_half_edges(brep, face)? {
                let edge = brep.half_edges[he].edge;
                let Some(curve) = brep.edges[edge].curve.as_ref() else {
                    continue;
                };
                let (e0, e1) = (curve.t_start, curve.t_end);
                if !(e0.is_finite() && e1.is_finite()) {
                    continue;
                }
                let edge_eval: &dyn CurveEval = &curve.geometry;
                let mut prev = edge_eval.point_at(e0);
                for i in 1..=64 {
                    let q = edge_eval.point_at(e0 + (e1 - e0) * i as f64 / 64.0);
                    let (d, foot) = point_segment_distance(p_end, prev, q);
                    if best.as_ref().map(|(bd, _, _)| d < *bd).unwrap_or(true) {
                        best = Some((d, edge, foot));
                    }
                    prev = q;
                }
            }
            if let Some((d, edge, point)) = best
                && d <= band
            {
                if trace {
                    eprintln!(
                        "    [crossings] endpoint contact at arc t {t_missing:.1} \
                         d {d:.3e} — synthesized (exact-corner class)"
                    );
                }
                hits.push((t_missing, edge, point, d));
            }
        }
        hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }
    if hits.len() != 2 {
        return Err(format!(
            "expected 2 arc/boundary crossings, found {} ({:?})",
            hits.len(),
            hits.iter().map(|h| (h.0, h.3)).collect::<Vec<_>>()
        ));
    }
    Ok([
        (hits[0].0, hits[0].1, hits[0].2),
        (hits[1].0, hits[1].1, hits[1].2),
    ])
}

/// Min distance from `p` to segment `a..b`, with the foot point.
fn point_segment_distance(p: Point3, a: Point3, b: Point3) -> (f64, Point3) {
    let ab = b - a;
    let len2 = ab.dot(&ab);
    let s = if len2 > 0.0 {
        ((p - a).dot(&ab) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let foot = Point3::new(a.x + s * ab.x, a.y + s * ab.y, a.z + s * ab.z);
    (p.distance_to(&foot), foot)
}

/// Fraction of the arc's interior samples that lie INSIDE the
/// face's trim, judged in the face surface's UV chart (even-odd ray
/// cast against the projected outer loop). Distinguishes trims that
/// share one underlying surface — surface proximity cannot. Returns
/// `None` when the face has no surface or projection fails.
pub(super) fn arc_inside_trim_fraction(
    brep: &BRep,
    face: FaceHandle,
    cubic: &[Point3; 4],
) -> Option<f64> {
    let surface = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|sof| sof.geometry.clone())?;
    let u_period = match &surface {
        popkern_core::geometry::SurfaceGeometry::Cylinder { .. }
        | popkern_core::geometry::SurfaceGeometry::Torus { .. }
        | popkern_core::geometry::SurfaceGeometry::Sphere { .. } => Some(std::f64::consts::TAU),
        _ => None,
    };
    let arc = bezier_cubic_curve(cubic);
    let arc_eval: &dyn CurveEval = &popkern_core::geometry::CurveGeometry::Nurbs(arc);
    let arc_uv = project_polyline(&surface, u_period, |s| arc_eval.point_at(s), 16)?;
    let arc_u_mean = mean_u(&arc_uv);

    // Assemble the outer loop as one closed UV polygon, unwrapped by
    // continuity point-to-point and mean-shifted to the arc's branch.
    let mut polygon: Vec<(f64, f64)> = Vec::new();
    for he in face_loop_half_edges(brep, face).ok()? {
        let half_edge = brep.half_edges.get(he)?;
        let edge = brep.edges.get(half_edge.edge)?;
        let curve = edge.curve.as_ref()?;
        if !(curve.t_start.is_finite() && curve.t_end.is_finite()) {
            return None;
        }
        let forward = edge.half_edge == he;
        let (a, b) = if forward {
            (curve.t_start, curve.t_end)
        } else {
            (curve.t_end, curve.t_start)
        };
        let eval: &dyn CurveEval = &curve.geometry;
        for i in 0..16 {
            let s = i as f64 / 16.0;
            let p = eval.point_at(a + (b - a) * s);
            let (mut u, v) =
                popkern_core::geometry::uv_projection::project_point_to_uv(&surface, &p)?;
            if let (Some(period), Some(prev)) = (u_period, polygon.last()) {
                let k = ((prev.0 - u) / period).round();
                u += k * period;
            }
            polygon.push((u, v));
        }
    }
    if polygon.len() < 3 {
        return None;
    }
    if let Some(period) = u_period {
        let poly_mean = polygon.iter().map(|p| p.0).sum::<f64>() / polygon.len() as f64;
        let shift = ((arc_u_mean - poly_mean) / period).round() * period;
        for p in &mut polygon {
            p.0 += shift;
        }
    }
    // Interior samples only — the arc's ends live at the trim corner
    // by design and would count half-in/half-out as noise.
    let interior: Vec<&(f64, f64, f64)> =
        arc_uv.iter().filter(|p| p.2 > 0.2 && p.2 < 0.8).collect();
    if interior.is_empty() {
        return None;
    }
    let inside = interior
        .iter()
        .filter(|p| point_in_polygon((p.0, p.1), &polygon))
        .count();
    Some(inside as f64 / interior.len() as f64)
}

fn point_in_polygon(point: (f64, f64), polygon: &[(f64, f64)]) -> bool {
    let (x, y) = point;
    let mut inside = false;
    let n = polygon.len();
    for i in 0..n {
        let (x0, y0) = polygon[i];
        let (x1, y1) = polygon[(i + 1) % n];
        if (y0 > y) != (y1 > y) {
            let x_cross = x0 + (y - y0) / (y1 - y0) * (x1 - x0);
            if x < x_cross {
                inside = !inside;
            }
        }
    }
    inside
}

/// Project a 3-D curve into the surface's UV chart as a polyline
/// over `s ∈ [0, 1]`, unwrapped by continuity in `u` when the chart
/// is periodic. Returns `None` if any sample fails to project.
fn project_polyline(
    surface: &popkern_core::geometry::SurfaceGeometry,
    u_period: Option<f64>,
    point_at: impl Fn(f64) -> Point3,
    samples: usize,
) -> Option<Vec<(f64, f64, f64)>> {
    let mut out: Vec<(f64, f64, f64)> = Vec::with_capacity(samples + 1);
    for i in 0..=samples {
        let s = i as f64 / samples as f64;
        let p = point_at(s);
        let (mut u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(surface, &p)?;
        if let (Some(period), Some(prev)) = (u_period, out.last()) {
            // Continuity unwrap: pick the branch nearest the
            // previous sample.
            let k = ((prev.0 - u) / period).round();
            u += k * period;
        }
        out.push((u, v, s));
    }
    Some(out)
}

fn mean_u(polyline: &[(f64, f64, f64)]) -> f64 {
    if polyline.is_empty() {
        return 0.0;
    }
    polyline.iter().map(|p| p.0).sum::<f64>() / polyline.len() as f64
}

/// All transversal intersections between two UV polylines, returned
/// as `(s_a, s_b)` curve parameters interpolated within the crossing
/// segments.
fn polyline_intersections(a: &[(f64, f64, f64)], b: &[(f64, f64, f64)]) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for wa in a.windows(2) {
        let (ax0, ay0, as0) = wa[0];
        let (ax1, ay1, as1) = wa[1];
        for wb in b.windows(2) {
            let (bx0, by0, bs0) = wb[0];
            let (bx1, by1, bs1) = wb[1];
            let d = (ax1 - ax0) * (by1 - by0) - (ay1 - ay0) * (bx1 - bx0);
            if d.abs() < 1e-30 {
                continue;
            }
            let s = ((bx0 - ax0) * (by1 - by0) - (by0 - ay0) * (bx1 - bx0)) / d;
            let t = ((bx0 - ax0) * (ay1 - ay0) - (by0 - ay0) * (ax1 - ax0)) / d;
            if (0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&t) {
                out.push((as0 + (as1 - as0) * s, bs0 + (bs1 - bs0) * t));
            }
        }
    }
    out
}
