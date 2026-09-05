//! Closed-record parameterization rotation for the cellular boolean.
//!
//! A closed intersection curve must CLOSE at an edge-face hit point
//! when one lies on it: the loop walk enters/bounces off an
//! intersection circle only at its closure node, so a record whose
//! closure sits elsewhere (an offset walk start, a swept boss with a
//! remote seam) is rotated to that convention before graph insertion.
//! Split from `cellular/boolean.rs` under the structural line budget.

/// Rotate a CLOSED intersection record's parameterization so its
/// closure point coincides with an edge-face hit lying on the curve.
/// Open records and closed records with no hit on them pass through
/// unchanged.
pub(super) fn rotate_closed_record_to_hits(
    record: &crate::intersect::IntersectionRecord,
    hits: &[crate::intersect_edge_face::EdgeFaceHit],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> crate::intersect::IntersectionRecord {
    use popkern_core::geometry::curve::CurveEval as _;
    let c = &record.curve.curve_3d;
    let (t0, t1) = c.parameter_range();
    // Rotation is only meaningful for a genuinely CLOSED record. The
    // record's own flag is authoritative — the clip stage keeps it
    // honest (a multi-window clip emits open arcs flagged open). The
    // geometric gap check remains as a torn-closure guard, NOT as the
    // closure test itself: a pure distance gate at containment scale
    // (linear × 100) swallowed open corner arcs whose tip gap was
    // 0.17 and "rotated" them into gap-bridging garbage curves
    // (measured: cube∪sphere grew 12 arcs running through out-of-face
    // regions and the assembly went non-manifold).
    let closed = record.curve.is_closed
        && c.point_at(t0).distance_to(&c.point_at(t1)) < tolerance.linear() * 100.0;
    if !closed {
        return record.clone();
    }
    if popkern_core::trace::trace_enabled("boolean") {
        let gap = c.point_at(t0).distance_to(&c.point_at(t1));
        eprintln!(
            "  [rotate] closed record considered for hit rotation: gap={gap:.3e} t[{t0:.3},{t1:.3}]"
        );
    }
    for hit in hits {
        // Nearest parameter on the curve to the hit point.
        let mut best = (f64::INFINITY, t0);
        for k in 0..=256 {
            let t = t0 + (t1 - t0) * k as f64 / 256.0;
            let d = c.point_at(t).distance_to(&hit.point);
            if d < best.0 {
                best = (d, t);
            }
        }
        let (mut lo, mut hi) = (
            (best.1 - (t1 - t0) / 256.0).max(t0),
            (best.1 + (t1 - t0) / 256.0).min(t1),
        );
        for _ in 0..48 {
            let m1 = lo + (hi - lo) / 3.0;
            let m2 = hi - (hi - lo) / 3.0;
            if c.point_at(m1).distance_to(&hit.point) < c.point_at(m2).distance_to(&hit.point) {
                hi = m2;
            } else {
                lo = m1;
            }
        }
        let t_hit = 0.5 * (lo + hi);
        let d_hit = c.point_at(t_hit).distance_to(&hit.point);
        if popkern_core::trace::trace_enabled("boolean") {
            eprintln!(
                "  [rotate] hit ({:.3},{:.3},{:.3}) nearest t={t_hit:.4} d={d_hit:.3e}",
                hit.point.x, hit.point.y, hit.point.z
            );
        }
        if d_hit > tolerance.linear() * 100.0 {
            continue;
        }
        // Interior rotation only.
        if (t_hit - t0).abs() < (t1 - t0) * 1e-4 || (t1 - t_hit).abs() < (t1 - t0) * 1e-4 {
            continue;
        }
        let mut out = record.clone();
        let rotate = |curve: &popkern_core::geometry::NurbsCurve|
            -> Option<popkern_core::geometry::NurbsCurve> {
            let (a, b) = curve.split_at(t_hit).ok()?;
            // New curve = b then a, reparameterized to [0, 1] with
            // the join at the original closure (multiplicity-degree
            // interior knot — legal C0 join).
            let (b0, b1) = {
                let n = &b;
                (n.knots[0], *n.knots.last()?)
            };
            let (a0, a1) = {
                let n = &a;
                (n.knots[0], *n.knots.last()?)
            };
            let len_b = b1 - b0;
            let len_a = a1 - a0;
            let total = len_a + len_b;
            if total < 1e-12 {
                return None;
            }
            let mut knots: Vec<f64> = b
                .knots
                .iter()
                .map(|&k| (k - b0) / total)
                .collect();
            let join = len_b / total;
            // Skip a's leading clamp knots (degree+1 of them) — the
            // join point carries multiplicity = degree from b's tail
            // clamp minus the shared control point.
            let a_knots_tail: Vec<f64> = a
                .knots
                .iter()
                .skip(a.degree + 1)
                .map(|&k| join + (k - a0) / total)
                .collect();
            // Replace b's trailing clamp (degree+1 copies) with
            // multiplicity-degree at the join.
            for _ in 0..1 {
                knots.pop();
            }
            knots.extend(a_knots_tail);
            let mut control_points = b.control_points.clone();
            // JUNCTION CONTINUITY: the splice joins the original END
            // (b's tail) to the original START (a's head). In 3D
            // those coincide, but a UV curve on a periodic chart
            // unwraps continuously, so its two ends differ by one
            // period — splicing raw control points leaves a full-
            // period jump mid-curve and every downstream chart
            // embedding of the rotated record is garbage (measured:
            // the cross-drill pierce loop's drill-chart pcurve).
            // Translate the appended piece so the junction matches;
            // for 3D curves the delta is ~zero and this is a no-op.
            let delta = match (b.control_points.last(), a.control_points.first()) {
                (Some(bl), Some(af)) => *bl - *af,
                _ => popkern_core::geometry::Vec3::zero(),
            };
            control_points.extend(a.control_points.iter().skip(1).map(|&p| p + delta));
            let mut weights = b.weights.clone();
            weights.extend(a.weights.iter().skip(1).copied());
            Some(popkern_core::geometry::NurbsCurve {
                degree: b.degree,
                knots,
                control_points,
                weights,
            })
        };
        if let (Some(c3), Some(u1), Some(u2)) = (
            rotate(&record.curve.curve_3d),
            rotate(&record.curve.curve_uv1),
            rotate(&record.curve.curve_uv2),
        ) {
            out.curve.curve_3d = c3;
            out.curve.curve_uv1 = u1;
            out.curve.curve_uv2 = u2;
            // Weld the new closure exactly onto the hit point.
            if let Some(cp) = out.curve.curve_3d.control_points.first_mut() {
                *cp = hit.point;
            }
            if let Some(cp) = out.curve.curve_3d.control_points.last_mut() {
                *cp = hit.point;
            }
            return out;
        }
    }
    record.clone()
}
