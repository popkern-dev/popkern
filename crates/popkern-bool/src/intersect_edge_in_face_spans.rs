//! The spans of an edge lying in a plane that fall inside a planar
//! face's region (see `intersect_edge_in_face`).
//!
//! Span ends are where the edge crosses the face's boundary: crossings
//! of the edge with the face's boundary edges in the plane, and the
//! face's vertices lying on the edge. Each is located on the 2D
//! polylines, then polished onto both curves by alternating closest-
//! point projection, so a span ends at the junction the graph splits
//! the edge at, not at a sample. Between consecutive ends the region
//! verdict is constant; three interior verdicts decide it.

use popkern_core::geometry::curve::{CurveEval, closest_point_in_window};
use popkern_core::geometry::{CurveGeometry, Point3, PortableMath, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Polyline samples along the edge.
const EDGE_POLY: usize = 128;
/// Polyline samples along each boundary edge of the face.
const BOUNDARY_POLY: usize = 64;

/// A plane chart: origin and unit axes.
pub(crate) struct PlaneChart {
    pub origin: Point3,
    pub u_axis: Vec3,
    pub v_axis: Vec3,
}

impl PlaneChart {
    fn uv(&self, p: &Point3) -> (f64, f64) {
        let d = *p - self.origin;
        (d.dot(&self.u_axis), d.dot(&self.v_axis))
    }
}

/// The parameter windows of `curve` over `[t0, t1]` inside face `fh`'s
/// region, in ascending order.
pub(crate) fn inside_spans(
    curve: &CurveGeometry,
    (t0, t1): (f64, f64),
    brep: &BRep,
    fh: FaceHandle,
    chart: &PlaneChart,
    tolerance: &Tolerance,
) -> Vec<(f64, f64)> {
    let lin = tolerance.linear();
    let inside = |t: f64| {
        let (u, v) = chart.uv(&curve.point_at(t));
        crate::classify_uv_loops::point_in_face_region_exact(brep, fh, u, v)
    };
    let boundary = boundary_curves(brep, fh);
    // On the face's own boundary: such a span runs along a boundary
    // edge (edge-on-edge coincidence, owned elsewhere) and divides
    // nothing.
    let on_boundary = |t: f64| {
        let p = curve.point_at(t);
        boundary.iter().any(|(bc, (b0, b1))| {
            bc.point_at(closest_point_in_window(bc, &p, *b0, *b1))
                .distance_to(&p)
                <= lin
        })
    };
    let mut breaks = vec![t0, t1];
    let edge_params: Vec<f64> = (0..=EDGE_POLY)
        .map(|k| t0 + (t1 - t0) * k as f64 / EDGE_POLY as f64)
        .collect();
    let edge_uv: Vec<(f64, f64)> = edge_params
        .iter()
        .map(|&t| chart.uv(&curve.point_at(t)))
        .collect();
    for (bc, (b0, b1)) in &boundary {
        let (bc, b0, b1) = (bc, *b0, *b1);
        // Face vertices on the edge.
        for tb in [b0, b1] {
            let q = bc.point_at(tb);
            let te = closest_point_in_window(curve, &q, t0, t1);
            if curve.point_at(te).distance_to(&q) <= lin {
                breaks.push(te);
            }
        }
        // Proper crossings of the two polylines.
        let bparams: Vec<f64> = (0..=BOUNDARY_POLY)
            .map(|k| b0 + (b1 - b0) * k as f64 / BOUNDARY_POLY as f64)
            .collect();
        let buv: Vec<(f64, f64)> = bparams.iter().map(|&t| chart.uv(&bc.point_at(t))).collect();
        for i in 0..EDGE_POLY {
            for j in 0..BOUNDARY_POLY {
                let Some((s, r)) = segment_cross(edge_uv[i], edge_uv[i + 1], buv[j], buv[j + 1])
                else {
                    continue;
                };
                let te = edge_params[i] + s * (edge_params[i + 1] - edge_params[i]);
                let tb = bparams[j] + r * (bparams[j + 1] - bparams[j]);
                let (te, tb) = polish(curve, (t0, t1), bc, (b0, b1), te, tb);
                if curve.point_at(te).distance_to(&bc.point_at(tb)) <= lin {
                    breaks.push(te);
                }
            }
        }
    }
    breaks.sort_by(f64::total_cmp);
    breaks.dedup_by(|x, y| curve.point_at(*x).distance_to(&curve.point_at(*y)) <= lin);
    let mut spans: Vec<(f64, f64)> = Vec::new();
    for w in breaks.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        if curve.point_at(lo).distance_to(&curve.point_at(hi)) <= lin {
            continue;
        }
        let probes = [0.25, 0.5, 0.75].map(|f| lo + (hi - lo) * f);
        let all_in = probes.iter().all(|&t| inside(t) == Some(true));
        if !all_in || probes.iter().all(|&t| on_boundary(t)) {
            continue;
        }
        match spans.last_mut() {
            Some(last) if last.1 == lo => last.1 = hi,
            _ => spans.push((lo, hi)),
        }
    }
    spans
}

/// The curves of face `fh`'s boundary edges, each over its edge window.
fn boundary_curves(brep: &BRep, fh: FaceHandle) -> Vec<(CurveGeometry, (f64, f64))> {
    let mut out = Vec::new();
    let Some(face) = brep.faces.get(fh) else {
        return out;
    };
    for &lh in std::iter::once(&face.outer_loop).chain(face.inner_loops.iter()) {
        let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
            continue;
        };
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(h) = brep.half_edges.get(he) else {
                break;
            };
            if let Some(c) = brep.edges.get(h.edge).and_then(|e| e.curve.as_ref()) {
                out.push((c.geometry.clone(), (c.t_start, c.t_end)));
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    out
}

/// The fractions along `p0`–`p1` and `q0`–`q1` where they cross, if
/// they cross (endpoints included).
fn segment_cross(
    p0: (f64, f64),
    p1: (f64, f64),
    q0: (f64, f64),
    q1: (f64, f64),
) -> Option<(f64, f64)> {
    let r = (p1.0 - p0.0, p1.1 - p0.1);
    let s = (q1.0 - q0.0, q1.1 - q0.1);
    let den = r.0 * s.1 - r.1 * s.0;
    if den.abs() <= f64::EPSILON * (r.0.phypot(r.1) * s.0.phypot(s.1)).max(f64::MIN_POSITIVE) {
        return None;
    }
    let d = (q0.0 - p0.0, q0.1 - p0.1);
    let a = (d.0 * s.1 - d.1 * s.0) / den;
    let b = (d.0 * r.1 - d.1 * r.0) / den;
    ((0.0..=1.0).contains(&a) && (0.0..=1.0).contains(&b)).then_some((a, b))
}

/// Alternating closest-point projection between two crossing curves,
/// from a polyline estimate, within each curve's window.
fn polish(
    e: &CurveGeometry,
    ew: (f64, f64),
    b: &CurveGeometry,
    bw: (f64, f64),
    mut te: f64,
    mut tb: f64,
) -> (f64, f64) {
    for _ in 0..60 {
        let nb = closest_point_in_window(b, &e.point_at(te), bw.0, bw.1);
        let ne = closest_point_in_window(e, &b.point_at(nb), ew.0, ew.1);
        let moved = (ne - te).abs() + (nb - tb).abs();
        te = ne;
        tb = nb;
        if moved <= f64::EPSILON * (1.0 + te.abs() + tb.abs()) {
            break;
        }
    }
    (te, tb)
}
