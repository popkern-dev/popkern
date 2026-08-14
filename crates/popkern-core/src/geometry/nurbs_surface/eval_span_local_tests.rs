//! Bit-identity pin for the span-local tensor-product evaluation
//! (PK-281): the local-window path must reproduce the full-grid
//! de Boor evaluation EXACTLY (`==`, no tolerance), because the
//! recursion runs the same lerps on the same values in the same
//! order — only the wasted rows/columns are gone. Any future edit
//! that breaks bit-equality changes every downstream trajectory
//! (projection footpoints, pcurve fits) and must be a deliberate
//! decision, not a silent drift.

use super::NurbsSurface;
use crate::geometry::Point3;
use crate::geometry::nurbs_ops::{Homogeneous, de_boor_1d};
use crate::geometry::surface::SurfaceEval;

/// Reference: the pre-PK-281 full-grid evaluation (every row
/// evaluated via `de_boor_1d`, which does its own span search).
fn full_grid_point_at(s: &NurbsSurface, u: f64, v: f64) -> Point3 {
    let mut intermediate: Vec<Homogeneous> = Vec::with_capacity(s.num_u);
    for i_u in 0..s.num_u {
        let row: Vec<Homogeneous> = (0..s.num_v)
            .map(|i_v| Homogeneous::from_point_weight(s.cp(i_u, i_v), s.wt(i_u, i_v)))
            .collect();
        intermediate.push(de_boor_1d(s.degree_v, &s.knots_v, &row, v));
    }
    de_boor_1d(s.degree_u, &s.knots_u, &intermediate, u).to_point3()
}

/// A rational multi-span surface: degree (3, 2), 7×5 control points,
/// interior knots in both directions, non-uniform weights. Multi-span
/// is the case where span-local and full-grid could diverge if the
/// window or span bookkeeping were wrong.
fn multi_span_rational() -> NurbsSurface {
    let num_u = 7;
    let num_v = 5;
    let mut control_points = Vec::with_capacity(num_u * num_v);
    let mut weights = Vec::with_capacity(num_u * num_v);
    for i in 0..num_u {
        for j in 0..num_v {
            let (fi, fj) = (i as f64, j as f64);
            control_points.push(Point3::new(
                fi + 0.3 * (fj * 1.7).sin(),
                fj + 0.2 * (fi * 2.3).cos(),
                0.5 * ((fi - 3.0) * (fj - 2.0)) + 0.1 * fi,
            ));
            weights.push(1.0 + 0.4 * ((i * 5 + j * 3) % 7) as f64 / 7.0);
        }
    }
    NurbsSurface {
        degree_u: 3,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 0.4, 0.7, 1.0, 1.0, 1.0],
        control_points,
        num_u,
        num_v,
        weights,
    }
}

#[test]
fn span_local_eval_bit_identical_to_full_grid() {
    let s = multi_span_rational();
    let n = 23;
    for iu in 0..=n {
        for iv in 0..=n {
            let u = iu as f64 / n as f64;
            let v = iv as f64 / n as f64;
            let fast = s.point_at(u, v);
            let reference = full_grid_point_at(&s, u, v);
            assert!(
                fast.x == reference.x && fast.y == reference.y && fast.z == reference.z,
                "span-local eval diverged from full-grid at (u={u}, v={v}): \
                 fast=({}, {}, {}) reference=({}, {}, {})",
                fast.x,
                fast.y,
                fast.z,
                reference.x,
                reference.y,
                reference.z
            );
        }
    }
}

#[test]
fn span_local_eval_bit_identical_at_knots_and_ends() {
    let s = multi_span_rational();
    // Exact knot values and domain endpoints exercise the span-clamp
    // branches of find_span.
    for &u in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        for &v in &[0.0, 0.4, 0.7, 1.0] {
            let fast = s.point_at(u, v);
            let reference = full_grid_point_at(&s, u, v);
            assert!(
                fast.x == reference.x && fast.y == reference.y && fast.z == reference.z,
                "span-local eval diverged at knot (u={u}, v={v})"
            );
        }
    }
}
