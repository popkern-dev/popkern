//! Seed classification for the region path (PK-141): arc-length
//! span fractions for descendant edges against the full trim curve.

use popkern_core::geometry::NurbsCurve;
use popkern_core::topology::{BRep, EdgeHandle};

/// ARC-LENGTH fraction along `trim` (0..1) of the midpoint of edge
/// `e`, with the distance to the trim. Parameter fractions are
/// unusable here — trim carriers can be padded/rational, so equal
/// parameter steps are not equal lengths (the param-fraction filter
/// classified an endpoint as frac 0.5). `None` when handles are
/// stale.
pub(in crate::construct) fn descendant_trim_fraction(
    brep: &BRep,
    e: EdgeHandle,
    trim: &NurbsCurve,
) -> Option<(f64, f64)> {
    use popkern_core::geometry::curve::CurveEval;
    let ed = brep.edges.get(e)?;
    let he = brep.half_edges.get(ed.half_edge)?;
    let tw = brep.half_edges.get(he.twin)?;
    let a = brep.vertices.get(he.origin)?.position;
    let b = brep.vertices.get(tw.origin)?.position;
    let m = popkern_core::geometry::Point3 {
        x: 0.5 * (a.x + b.x),
        y: 0.5 * (a.y + b.y),
        z: 0.5 * (a.z + b.z),
    };
    let (t0, t1) = trim.parameter_range();
    let n = 128;
    let mut lengths = Vec::with_capacity(n + 1);
    let mut total = 0.0;
    let mut prev = trim.point_at(t0);
    lengths.push(0.0);
    let mut best = (0usize, f64::MAX);
    for i in 0..=n {
        let t = t0 + (t1 - t0) * (i as f64) / (n as f64);
        let q = trim.point_at(t);
        if i > 0 {
            total += prev.distance_to(&q);
            lengths.push(total);
            prev = q;
        }
        let dm = m.distance_to(&q);
        if dm < best.1 {
            best = (i, dm);
        }
    }
    if total < 1e-12 {
        return None;
    }
    Some((lengths[best.0] / total, best.1))
}

/// Descendants of `fillet_desc` whose midpoints project to the
/// INTERIOR of the trim span (see [`descendant_trim_fraction`]) —
/// the span classification shared by region seeding and the direct
/// consumed identification.
pub(in crate::construct) fn mid_span_descendants(
    brep: &BRep,
    fillet_desc: &[EdgeHandle],
    trim: &NurbsCurve,
) -> Vec<EdgeHandle> {
    fillet_desc
        .iter()
        .copied()
        .filter(|&e| {
            descendant_trim_fraction(brep, e, trim)
                .map(|(frac, _)| (0.08..=0.92).contains(&frac))
                .unwrap_or(false)
        })
        .collect()
}
