//! Collinear-chain continuation for the anchor-split stage (PK-171).
//! Split from `perp_anchor_positions.rs` (file-size limit).

use popkern_core::geometry::{NurbsCurve, Point3};

/// The anchor on a COLLINEAR CONTINUATION of the walked edge, when the
/// anchor distance overshoots it: a prior boolean may have split the
/// support rim mid-span (measured: an L-bracket union splits the base
/// top rim 2.0 from the corner; radius 4 clamped the anchor onto the
/// junction vertex, the strip side trim slanted to it, the terminal
/// stub settle declined on the off-cylinder anchor, and the stale
/// corner stub shipped with deviation exactly r(√2−1) — check-9 catch,
/// PK-171). Walks the loop's segments `loop_segs` from `loop_index`
/// while they stay collinear with the walked edge `p0 → p1`, trying
/// the trim crossing (`trim`, when trusted) on each piece and spending
/// the remaining radius. Returns the placement and its slot quality
/// (2 = trim crossing, 1 = walk), or `None` when the chain ends at a
/// genuine corner first.
pub(super) fn collinear_chain_anchor(
    loop_segs: &[Option<(Point3, Point3)>],
    loop_index: usize,
    forward: bool,
    (p0, p1): (Point3, Point3),
    radius: f64,
    trim: Option<&NurbsCurve>,
) -> Option<(Point3, u8)> {
    let edge_len = p0.distance_to(&p1);
    let dir0 = (p1 - p0) * (1.0 / edge_len);
    let mut remaining = radius - edge_len;
    let n = loop_segs.len();
    for step in 1..n {
        let j = if forward {
            (loop_index + step) % n
        } else {
            (loop_index + n - step) % n
        };
        let (q0, q1) = loop_segs[j]?;
        // Orient the piece along the walk.
        let (a, b) = if forward { (q0, q1) } else { (q1, q0) };
        let len = a.distance_to(&b);
        if len <= 1e-15 {
            continue;
        }
        let dir = (b - a) * (1.0 / len);
        let walk_dir = if forward { dir0 } else { dir0 * -1.0 };
        if dir.dot(&walk_dir) < 1.0 - 1e-9 {
            return None; // chain ends at a genuine corner
        }
        if let Some(curve) = trim
            && let Some((s, _)) =
                crate::construct::cyl_plane_intersect::intersect_curve_with_segment(curve, q0, q1)
            && (0.001..0.999).contains(&s)
        {
            return Some((
                Point3::new(
                    q0.x + s * (q1.x - q0.x),
                    q0.y + s * (q1.y - q0.y),
                    q0.z + s * (q1.z - q0.z),
                ),
                2,
            ));
        }
        if remaining <= len * 0.999 {
            let f = (remaining / len).max(0.001);
            return Some((
                Point3::new(
                    a.x + f * (b.x - a.x),
                    a.y + f * (b.y - a.y),
                    a.z + f * (b.z - a.z),
                ),
                1,
            ));
        }
        remaining -= len;
    }
    None
}
