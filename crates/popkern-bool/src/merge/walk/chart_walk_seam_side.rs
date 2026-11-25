//! A ray at a seam branch point lies on its own side of the seam.
//!
//! A node on a periodic chart's seam splits into two branch points a
//! period apart in u (the seam's low and high copies), and the face
//! lies on one side of each: above u at the low copy, below it at the
//! high one. Every ray the wedge matching sorts there leaves into that
//! half-plane. An endpoint-refined δ-sample of a curve running along
//! the seam within its projection's resolution can land a hair across
//! it, and the ray then points out of the chart: sorted through the
//! outside half-plane, it falls on the wrong side of the seam's own
//! rays and the parenthesis matching pairs its arrival with its own
//! twin (generated corpus SW-V scene 836: a tangent crossing on a
//! sphere's seam whose branch along the seam pointed 0.08 rad across
//! it; the walk turned back on it and the loop broke). Such a ray runs
//! along the seam, so it is read along the seam, a hair inside.

use super::super::chart_subdivision::{BranchPoint, FaceChartSubdivision};

/// How far inside the half-plane a ray read along the seam is turned:
/// enough to order it strictly inside the seam's own ray, far under
/// any angle the matching otherwise distinguishes.
const INSIDE_TURN: f64 = 1.0e-9;

/// The u direction into the face at `bp` (+1 at the seam's low copy,
/// −1 at its high copy), when `bp` is one of a seam node's two copies;
/// `None` elsewhere.
pub(super) fn seam_side_inward(sub: &FaceChartSubdivision, bp: &BranchPoint) -> Option<f64> {
    let period = sub.u_period?;
    let u_at = |b: &BranchPoint| -> Option<f64> {
        sub.departures
            .get(b)?
            .iter()
            .find_map(|c| sub.embeddings.get(c).and_then(|e| e.uv.first()))
            .map(|uv| uv.0)
    };
    let u = u_at(bp)?;
    sub.departures
        .keys()
        .filter(|b| b.node == bp.node && b.branch != bp.branch)
        .filter_map(u_at)
        .find_map(|other| {
            if other > u + 0.5 * period {
                Some(1.0)
            } else if other < u - 0.5 * period {
                Some(-1.0)
            } else {
                None
            }
        })
}

/// `d` when it leaves into the side `inward` names; read along the
/// seam, a hair inside, when it points across it.
pub(super) fn into_seam_side(d: (f64, f64), inward: f64) -> (f64, f64) {
    if d.0 * inward >= 0.0 || d.1 == 0.0 {
        return d;
    }
    let along = d.1.signum();
    let n = (INSIDE_TURN * INSIDE_TURN + 1.0).sqrt();
    (inward * INSIDE_TURN / n, along / n)
}

#[cfg(test)]
mod tests {
    use super::into_seam_side;

    #[test]
    fn a_ray_inside_its_side_is_kept() {
        assert_eq!(into_seam_side((0.3, 0.9), 1.0), (0.3, 0.9));
        assert_eq!(into_seam_side((-0.3, 0.9), -1.0), (-0.3, 0.9));
        assert_eq!(into_seam_side((0.0, -1.0), 1.0), (0.0, -1.0));
    }

    #[test]
    fn a_ray_across_the_seam_is_read_along_it_inside() {
        let (du, dv) = into_seam_side((-0.08, 0.997), 1.0);
        assert!(du > 0.0 && du < 1e-8, "{du}");
        assert!((dv - 1.0).abs() < 1e-12, "{dv}");
        let (du, dv) = into_seam_side((0.2, -0.98), -1.0);
        assert!(du < 0.0 && du > -1e-8, "{du}");
        assert!((dv + 1.0).abs() < 1e-12, "{dv}");
    }
}
