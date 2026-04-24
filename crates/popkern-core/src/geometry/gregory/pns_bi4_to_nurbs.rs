//! Convert a PnS bi-4 Bézier patch to a `NurbsSurface`.
//!
#![allow(clippy::needless_range_loop)]

//! A bi-4 Bézier patch on `[0, 1]²` is exactly a degree-(4, 4)
//! NURBS with 5 × 5 control points, a clamped uniform knot
//! vector, and unit weights. This is a one-liner-ish conversion
//! — no denominator-clearing (Gregory patches need that because
//! of the rational blend; PnS's bi-4 is plain Bézier).

use crate::geometry::Point3;
use crate::geometry::nurbs_surface::NurbsSurface;

/// Convert a PnS bi-4 cap sector patch `[[Point3; 5]; 5]` (row
/// index = `u_c`, column index = `v_c`) into a clamped
/// degree-(4, 4) NURBS surface on `[0, 1]²`.
///
/// The NURBS surface's `(u, v)` parameter matches the bi-4 patch's
/// canonical `(u_c, v_c)`.
pub fn pns_bi4_to_nurbs(patch: &[[Point3; 5]; 5]) -> NurbsSurface {
    // Clamped knot vector for degree 4 with 5 CPs: 5 zeros + 5 ones.
    let clamped_knots: Vec<f64> = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0].to_vec();

    // Flatten row-major: idx(u_c, v_c) = 5·u_c + v_c.
    let mut control_points: Vec<Point3> = Vec::with_capacity(25);
    for u_c in 0..5 {
        for v_c in 0..5 {
            control_points.push(patch[u_c][v_c]);
        }
    }
    let weights: Vec<f64> = vec![1.0; 25];

    NurbsSurface {
        degree_u: 4,
        degree_v: 4,
        knots_u: clamped_knots.clone(),
        knots_v: clamped_knots,
        control_points,
        num_u: 5,
        num_v: 5,
        weights,
    }
}

#[cfg(test)]
#[path = "pns_bi4_to_nurbs_tests.rs"]
mod tests;
