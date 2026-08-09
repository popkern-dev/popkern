//! Terminal-sliver-free conic arc segmentation (PK-257).
//!
//! Upstream angle authors carry atan2 noise (~1e-8 rad off exact
//! quarter angles even when the section vertices are exact). The
//! historical arc/ellipse conversions split a fixed quarter-knot
//! full circle/ellipse down to the requested range, so an endpoint
//! within that noise of the quarter grid stranded a nano piece with
//! a full-multiplicity break — sub-tolerance geometry that poisoned
//! every downstream knot union (PK-255: fillet-end seam arcs of
//! width π/2 ± 2e-8 collapsed the exact-boundary derivative rows).
//! Colliding ranges take [`uniform_segmented_conic_arc`]; everything
//! else keeps the historical byte-identical constructions their
//! consumers pair against.

use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Vec3};

/// Angular window (radians) for the quarter-grid collision tests.
/// Well above the ~1e-8 atan2 noise class carried by upstream angle
/// authors, well below any legitimate arc feature; a stranded piece
/// this size is sub-tolerance geometry for kernel-scale radii.
const QUARTER_COLLISION_WINDOW: f64 = 1.0e-6;

/// True when a quarter-grid knot anchored at `grid_origin` falls
/// strictly inside `(a0, a1)` within [`QUARTER_COLLISION_WINDOW`] of
/// either end — the historical split construction would strand a
/// terminal sliver there.
pub(super) fn quarter_grid_collides(grid_origin: f64, a0: f64, a1: f64) -> bool {
    let quarter = std::f64::consts::TAU / 4.0;
    (1..4).any(|index| {
        let knot = grid_origin + index as f64 * quarter;
        knot > a0
            && knot < a1
            && ((knot - a0).abs() <= QUARTER_COLLISION_WINDOW
                || (a1 - knot).abs() <= QUARTER_COLLISION_WINDOW)
    })
}

/// Uniformly segmented rational-quadratic conic arc: `n = ceil(span /
/// (π/2))` pieces of `span/n` each over exactly `[a0, a1]`, knots in
/// absolute angle, shoulder construction per piece (exact on circles
/// and — by affine invariance — ellipses). Sliver spans are
/// impossible by construction: every piece has width `span/n`.
pub(super) fn uniform_segmented_conic_arc(
    center: Point3,
    u_basis: Vec3,
    v_basis: Vec3,
    a0: f64,
    a1: f64,
) -> NurbsCurve {
    use std::f64::consts::TAU;
    let quarter = TAU / 4.0;
    let span = a1 - a0;
    let segments = ((span / quarter).ceil() as usize).max(1);
    let piece = span / segments as f64;
    let point_at_angle = |angle: f64, scale: f64| -> Point3 {
        center + u_basis * (scale * angle.cos()) + v_basis * (scale * angle.sin())
    };
    let shoulder_scale = {
        let half = 0.5 * piece;
        if half.cos().abs() > 1.0e-300 {
            1.0 / half.cos()
        } else {
            1.0
        }
    };
    let weight = (0.5 * piece).cos();
    let mut control_points = Vec::with_capacity(2 * segments + 1);
    let mut weights = Vec::with_capacity(2 * segments + 1);
    let mut knots = Vec::with_capacity(2 * segments + 4);
    knots.push(a0);
    for segment in 0..segments {
        let seg_start = a0 + piece * segment as f64;
        knots.push(seg_start);
        knots.push(seg_start);
        control_points.push(point_at_angle(seg_start, 1.0));
        weights.push(1.0);
        control_points.push(point_at_angle(seg_start + 0.5 * piece, shoulder_scale));
        weights.push(weight);
    }
    control_points.push(point_at_angle(a1, 1.0));
    weights.push(1.0);
    knots.push(a1);
    knots.push(a1);
    knots.push(a1);
    NurbsCurve {
        degree: 2,
        knots,
        control_points,
        weights,
    }
}
