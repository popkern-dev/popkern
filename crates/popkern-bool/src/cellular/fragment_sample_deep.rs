//! A decisive fragment sample: one clear of the other operand by more
//! than the tolerance (W3 V cells, tangency design §4.8).
//!
//! Where two faces cross at their touch, the fragments meeting there
//! lie within the tolerance of the other operand near the touch and
//! along their records, and a sample placed there classifies by noise:
//! across the generated corpus's tangent crossings, thirteen scenes
//! sampled a fragment within one tolerance of the other operand (2e-8
//! to 8e-7), read it on the wrong side, and their loops broke where the
//! kept and dropped fragments met. A fragment is one class throughout,
//! so any point of it decides; the point of a chart grid over it that
//! clears the other operand farthest decides best.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;

use super::fragment_classify_one::ClassifyCtx;
use crate::merge::graph::Operand;
use crate::merge::walk::Fragment;

/// Grid stations per chart direction over the fragment's box.
const GRID: usize = 16;
/// A sample this many tolerances off the other operand is decisive.
const DECISIVE: f64 = 4.0;

/// `sample`, or the point of `fragment` clearing the other operand
/// farthest when `sample` lies within [`DECISIVE`] tolerances of it.
/// Fragments with holes keep their sample (their samplers place it off
/// the holes), as does any fragment whose chart polygon is unavailable
/// or wraps its chart's period.
pub(super) fn decisive_sample(
    ctx: &ClassifyCtx<'_>,
    fragment: &Fragment,
    sample: Point3,
) -> Point3 {
    let (own, other, other_aabbs) = match fragment.face.operand {
        Operand::A => (ctx.brep_a, ctx.brep_b, ctx.aabbs_b),
        _ /* B */ => (ctx.brep_b, ctx.brep_a, ctx.aabbs_a),
    };
    let band = DECISIVE * ctx.tolerance.linear();
    let reach = loop_extent(ctx, fragment).max(band);
    let clearance = |p: &Point3| -> f64 {
        let probe =
            popkern_core::geometry::intersection::surface_subdivision::Aabb { min: *p, max: *p };
        other_aabbs
            .iter()
            .filter(|(_, b)| b.overlaps(&probe, reach))
            .filter_map(|(face, _)| {
                let s = &other.faces.get(*face)?.surface.as_ref()?.geometry;
                let (u, v) = s.closest_point(p)?;
                // The face, not its surface's extension: a sample on the
                // plane a neighbouring face extends is clear of it.
                let on_face = !matches!(
                    other.face_uv_trim_containment(*face, s, u, v),
                    Some(popkern_core::validation::FaceUvContainment::Outside)
                );
                on_face.then(|| s.point_at(u, v).distance_to(p))
            })
            .fold(f64::INFINITY, f64::min)
    };
    // A coincident-overlap face lies on the other operand throughout:
    // its offset-pair classification reads it, and no point clears it
    // (a swept bar's runs joined along shared faces read a moved sample
    // on the wrong side).
    let on_overlap = match fragment.face.operand {
        Operand::A => ctx.overlap_a.contains(&fragment.face.face),
        _ /* B */ => ctx.overlap_b.contains(&fragment.face.face),
    };
    let here = clearance(&sample);
    if here >= band || on_overlap || !fragment.inner_loops.is_empty() {
        return sample;
    }
    let Some(surface) = own
        .faces
        .get(fragment.face.face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    else {
        return sample;
    };
    let Some(chart) = crate::merge::fragment_sample_interior::chart_polygon(
        ctx.graph,
        &fragment.outer_loop,
        &fragment.face,
        surface,
    ) else {
        return sample;
    };
    let poly = &chart.polygon;
    let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
    for &(u, v) in poly {
        lo = (lo.0.min(u), lo.1.min(v));
        hi = (hi.0.max(u), hi.1.max(v));
    }
    let (pu, pv) = crate::clip::classification_chart_periods(surface);
    let wraps = |span: f64, period: Option<f64>, k: f64| {
        period.is_some_and(|p| p.is_finite() && p > 0.0 && span > 0.5 * p * k)
    };
    if wraps(hi.0 - lo.0, pu, chart.ku) || wraps(hi.1 - lo.1, pv, chart.kv) {
        return sample;
    }
    let mut best = (here, sample);
    for i in 0..GRID {
        for j in 0..GRID {
            let u = lo.0 + (hi.0 - lo.0) * (i as f64 + 0.5) / GRID as f64;
            let v = lo.1 + (hi.1 - lo.1) * (j as f64 + 0.5) / GRID as f64;
            if !contains(poly, (u, v)) {
                continue;
            }
            let p = surface.point_at(u / chart.ku, v / chart.kv);
            let c = clearance(&p);
            if c > best.0 {
                best = (c, p);
            }
        }
    }
    // A point past the linear reach of the other operand reads its side
    // decisively (one 3.8 tolerances clear decided a sector that its
    // first sample, 2.8e-7 off, read wrong).
    if best.0 > ctx.tolerance.linear_reach() {
        best.1
    } else {
        sample
    }
}

/// The diagonal of the fragment's outer loop nodes' box.
fn loop_extent(ctx: &ClassifyCtx<'_>, fragment: &Fragment) -> f64 {
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for &h in &fragment.outer_loop {
        let q = ctx.graph.nodes[ctx.graph.half_edges[h].origin].position;
        for (k, c) in [q.x, q.y, q.z].into_iter().enumerate() {
            lo[k] = lo[k].min(c);
            hi[k] = hi[k].max(c);
        }
    }
    (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt()
}

/// Even-odd containment of `q` in `poly`.
fn contains(poly: &[(f64, f64)], q: (f64, f64)) -> bool {
    let mut inside = false;
    let n = poly.len();
    let mut j = n - 1;
    for i in 0..n {
        let ((ui, vi), (uj, vj)) = (poly[i], poly[j]);
        if (vi > q.1) != (vj > q.1) && q.0 < (uj - ui) * (q.1 - vi) / (vj - vi) + ui {
            inside = !inside;
        }
        j = i;
    }
    inside
}
