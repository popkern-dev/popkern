//! Covered records: a record lying along another one of its face pair
//! is that one's duplicate authority and drops (see `drop_covered`).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval as _;
use popkern_core::tolerance::Tolerance;

use super::IntersectionRecord;

/// Drop, among the face pair's records `idxs` (with `info` their start,
/// end, middle and extent), each one covered by another.
pub(super) fn drop_covered(
    records: &[IntersectionRecord],
    idxs: &[usize],
    info: &[(Point3, Point3, Point3, f64)],
    gate: f64,
    on_curve: &dyn Fn(usize, &Point3) -> bool,
    drop: &mut [bool],
    tolerance: &Tolerance,
) {
    // COVERED records: junction-adjacent marching can deliver a
    // span AND overlapping micro-fragments of the same span from
    // different seeds (measured on the rotated torus × torus
    // union: pair (9v1,7v1) carried [0→0.025], [0→0.003] and
    // [0.005→0.025] of one junction arc — the walk chewed the
    // overlap into duplicate edge pairs and stranded halves).
    // A record whose start, mid and end all lie ON another
    // record's sampled polyline, with smaller-or-equal extent, is
    // duplicate authority and drops.
    for a in 0..idxs.len() {
        if drop[idxs[a]] {
            continue;
        }
        let poly_a: Vec<Point3> = {
            let c = &records[idxs[a]].curve.curve_3d;
            let (t0, t1) = c.parameter_range();
            (0..=32)
                .map(|k| c.point_at(t0 + (t1 - t0) * k as f64 / 32.0))
                .collect()
        };
        // Arc-length position of `p`'s closest on-polyline witness
        // within the gate, or None when off the polyline.
        let arc_pos = |p: &Point3| -> Option<f64> {
            let mut best: Option<(f64, f64)> = None; // (dist, arc)
            let mut acc = 0.0;
            for w in poly_a.windows(2) {
                let d = w[1] - w[0];
                let l = d.dot(&d).max(1e-300).sqrt();
                let t = ((*p - w[0]).dot(&d) / (l * l)).clamp(0.0, 1.0);
                let q = w[0] + d * t;
                let dist = q.distance_to(p);
                if dist < gate && best.is_none_or(|(bd, _)| dist < bd) {
                    best = Some((dist, acc + t * l));
                }
                acc += l;
            }
            best.map(|(_, arc)| arc)
        };
        for b in 0..idxs.len() {
            if a == b || drop[idxs[b]] {
                continue;
            }
            if info[b].3 > info[a].3 + gate {
                continue;
            }
            // A closed loop is never covered by an open record: the
            // record lacks the stretch the loop closes. A ball's
            // equator on a revolved cap came as a whole circle and
            // two arcs short of 0.02 rad at the coincident seams;
            // the circle's ends and middle lie on the longer arc,
            // it dropped, and the arcs left a hole (PK-1128).
            let closed = |k: usize| info[k].0.distance_to(&info[k].1) < gate && info[k].3 > gate;
            if closed(b) && !closed(a) {
                continue;
            }
            let (Some(p0), Some(p2), Some(p1)) = (
                arc_pos(&info[b].0),
                arc_pos(&info[b].2),
                arc_pos(&info[b].1),
            ) else {
                continue;
            };
            if !on_curve(idxs[a], &info[b].2) {
                continue;
            }
            // A duplicate OVERLAPS the polyline: its projections
            // spread along it comparably to its own extent. A
            // CHAINED continuation shorter than the gate projects
            // all three points onto one clamped endpoint (spread
            // ≈ 0) — dropping it severed a seam-sliver record and
            // stranded the assembly (measured: the spring-grind
            // section arc's 5.6e-4 piece across the wall seam,
            // PK-443).
            let spread = p0.max(p1).max(p2) - p0.min(p1).min(p2);
            if spread + tolerance.linear() >= 0.5 * info[b].3 {
                drop[idxs[b]] = true;
            }
        }
    }
}
