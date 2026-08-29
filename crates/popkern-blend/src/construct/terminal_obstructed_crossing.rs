//! Terminal-arc × cap-boundary crossing search for the obstructed
//! strip-end repair. Split from `terminal_obstructed_end.rs` to
//! honour the structural line ratchet.

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::terminal_arc_trim::CornerCandidate;
use crate::types::TerminalArcTrimDeclineReason;

/// Single crossing of the open terminal arc with the cap face's
/// straight boundary edges (stub chain and off-anchor edges
/// excluded). The arc-side angle is found by sign-change bisection
/// on the crossed edge's line; the hit must land inside the edge
/// segment. Exactly one crossing is the supported class.
pub(super) fn find_arc_crossing(
    brep: &BRep,
    cap_face: FaceHandle,
    c: &CornerCandidate,
    off_v: VertexHandle,
    arc_at: &dyn Fn(f64) -> Point3,
    sweep: f64,
    pos_tol: f64,
) -> Result<(popkern_core::topology::EdgeHandle, Point3), TerminalArcTrimDeclineReason> {
    let stub_edges: [popkern_core::topology::EdgeHandle; 2] = [
        brep.half_edges
            .get(c.he_in)
            .map(|he| he.edge)
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?,
        brep.half_edges
            .get(c.he_out)
            .map(|he| he.edge)
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?,
    ];
    let f = brep
        .faces
        .get(cap_face)
        .ok_or(TerminalArcTrimDeclineReason::NoCommonPlanarCap)?;
    let mut hits: Vec<(popkern_core::topology::EdgeHandle, Point3)> = Vec::new();
    let loops = std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied());
    for lh in loops {
        let Some(l) = brep.loops.get(lh) else {
            continue;
        };
        let start = l.half_edge;
        let mut cur = start;
        for _ in 0..brep.half_edges.len() + 1 {
            let Some(he) = brep.half_edges.get(cur) else {
                break;
            };
            let next_origin = brep
                .half_edges
                .get(he.next)
                .map(|n| n.origin)
                .unwrap_or(he.origin);
            let advance = he.next;
            let edge = he.edge;
            let skip = stub_edges.contains(&edge)
                || he.origin == c.corner
                || next_origin == c.corner
                || he.origin == off_v
                || next_origin == off_v;
            if !skip
                && let (Some(p0), Some(p1)) = (
                    brep.vertices.get(he.origin).map(|v| v.position),
                    brep.vertices.get(next_origin).map(|v| v.position),
                )
                && let Some(hit) = arc_segment_hit(arc_at, sweep, p0, p1, pos_tol)
                && !hits.iter().any(|(e, _)| *e == edge)
            {
                hits.push((edge, hit));
            }
            cur = advance;
            if cur == start {
                break;
            }
        }
    }
    match hits.len() {
        1 => Ok(hits.remove(0)),
        _ => Err(TerminalArcTrimDeclineReason::ObstructedEnd),
    }
}

/// Intersect the open arc `θ ∈ (0, sweep)` with segment `p0→p1` by
/// bisection on the signed distance to the segment's line, then
/// verify the hit lies inside the segment.
fn arc_segment_hit(
    arc_at: &dyn Fn(f64) -> Point3,
    sweep: f64,
    p0: Point3,
    p1: Point3,
    pos_tol: f64,
) -> Option<Point3> {
    let seg = p1 - p0;
    let seg_len = seg.length();
    if seg_len < 1e-12 {
        return None;
    }
    // Signed lateral offset from the segment's line, measured in the
    // plane containing the arc (use the arc chord to build a side
    // direction that is stable for end-plane geometry).
    let a0 = arc_at(1e-9);
    let a1 = arc_at(sweep - 1e-9);
    let chord = a1 - a0;
    let side_dir = seg.cross(&chord.cross(&seg));
    let side = |p: Point3| -> f64 { (p - p0).dot(&side_dir) };
    const SAMPLES: usize = 128;
    let mut prev_theta = sweep * 1e-9;
    let mut prev_side = side(arc_at(prev_theta));
    for i in 1..=SAMPLES {
        let theta = sweep * (i as f64) / (SAMPLES as f64 + 1e-9);
        let s = side(arc_at(theta));
        if prev_side == 0.0 || s.signum() != prev_side.signum() {
            // Bisect.
            let (mut lo, mut hi) = (prev_theta, theta);
            let (mut flo, _) = (prev_side, s);
            for _ in 0..64 {
                let mid = 0.5 * (lo + hi);
                let fm = side(arc_at(mid));
                if fm == 0.0 {
                    lo = mid;
                    hi = mid;
                    break;
                }
                if fm.signum() == flo.signum() {
                    lo = mid;
                    flo = fm;
                } else {
                    hi = mid;
                }
            }
            let hit = arc_at(0.5 * (lo + hi));
            // Inside the segment (with a small end margin), and
            // actually ON the segment's line.
            let t = (hit - p0).dot(&seg) / (seg_len * seg_len);
            let foot = p0 + seg * t;
            if (0.001..=0.999).contains(&t) && foot.distance_to(&hit) <= pos_tol.max(1e-7) * 10.0 {
                return Some(hit);
            }
        }
        prev_theta = theta;
        prev_side = s;
    }
    None
}
