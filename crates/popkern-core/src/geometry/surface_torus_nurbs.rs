//! Exact rational torus→NURBS conversion for partial angular windows
//! (PK-312).
//!
//! The full-period conversion lives in `surface.rs`
//! (`exact_full_torus_nurbs`, PK-310) and stays byte-stable; this
//! module extends exactness to sub-period windows, which previously
//! fell to the `sample_to_nurbs` fallback whose interpolation sag
//! shrank patch hulls up to 0.13 below the true surface and blinded
//! the SSI seeder's AABB cull (the measured PK-310 wrong-serve class).
//!
//! Construction: each axis's angular span decomposes into the minimum
//! number of equal arcs of at most 90°, each an exact rational
//! quadratic (on-curve endpoints weight 1; tangent-intersection
//! middle control point at radius `1/cos(half_angle)` with weight
//! `cos(half_angle)`). The two axis patterns tensor into a
//! surface-of-revolution patch exactly as in the full-period case.
//! The surface point set is exact (rational arcs lie on the circle)
//! and the control hull circumscribes the patch, so seeding boxes are
//! conservative. The NURBS parameter tracks the angle exactly at arc
//! boundaries and arc midpoints and deviates by a bounded skew in
//! between — the same contract the shipped full-period arm has, and
//! sufficient for every consumer (hull culling plus heuristic Newton
//! seeds refined against the original surfaces).

use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::{Point3, Vec3};

/// One axis's rational-arc pattern: control-point plane coordinates
/// `(cos·scale, sin·scale)`, weights, and the clamped quadratic knot
/// vector over the requested range.
struct AxisArcPattern {
    /// `(cos(θ)·scale, sin(θ)·scale)` per control point.
    plane: Vec<(f64, f64)>,
    weights: Vec<f64>,
    knots: Vec<f64>,
}

/// Decompose `[range.0, range.1]` (span ≤ 2π + 1e-9) into equal
/// rational-quadratic arcs of at most 90° each.
fn axis_arc_pattern(range: (f64, f64)) -> AxisArcPattern {
    let (a, b) = range;
    let span = b - a;
    let arcs = (span / std::f64::consts::FRAC_PI_2 - 1.0e-12)
        .ceil()
        .max(1.0) as usize;
    let step = span / arcs as f64;
    let w_mid = (step * 0.5).cos();
    let mid_scale = 1.0 / w_mid;
    let mut plane = Vec::with_capacity(2 * arcs + 1);
    let mut weights = Vec::with_capacity(2 * arcs + 1);
    for m in 0..arcs {
        let start = a + step * m as f64;
        plane.push((start.cos(), start.sin()));
        weights.push(1.0);
        let mid = start + step * 0.5;
        plane.push((mid.cos() * mid_scale, mid.sin() * mid_scale));
        weights.push(w_mid);
    }
    plane.push((b.cos(), b.sin()));
    weights.push(1.0);

    let mut knots = Vec::with_capacity(2 * arcs + 4);
    knots.extend([a, a, a]);
    for m in 1..arcs {
        let t = a + step * m as f64;
        knots.extend([t, t]);
    }
    knots.extend([b, b, b]);
    AxisArcPattern {
        plane,
        weights,
        knots,
    }
}

/// Exact rational NURBS of a torus patch over partial angular windows
/// (either span below the full period; spans up to 2π are accepted).
/// See the module docs for the construction and the exactness
/// contract.
pub(super) fn exact_torus_arc_nurbs(
    center: Point3,
    axis: Vec3,
    major_radius: f64,
    minor_radius: f64,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> NurbsSurface {
    let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(axis);
    let pu = axis_arc_pattern(u_range);
    let pv = axis_arc_pattern(v_range);
    let mut control_points = Vec::with_capacity(pu.plane.len() * pv.plane.len());
    let mut weights = Vec::with_capacity(pu.plane.len() * pv.plane.len());
    // Storage is u-major: index = i_u * num_v + i_v (see
    // `NurbsSurface::control_point`). u is the major angle, v the
    // profile angle.
    for (cu, wu) in pu.plane.iter().zip(pu.weights.iter()) {
        for (cv, wv) in pv.plane.iter().zip(pv.weights.iter()) {
            // Profile control point in (radial, axial) coordinates.
            let rho = major_radius + minor_radius * cv.0;
            let z = minor_radius * cv.1;
            let p = center + (perp1 * cu.0 + perp2 * cu.1) * rho + axis * z;
            control_points.push(p);
            weights.push(wu * wv);
        }
    }
    NurbsSurface {
        degree_u: 2,
        degree_v: 2,
        knots_u: pu.knots,
        knots_v: pv.knots,
        num_u: pu.plane.len(),
        num_v: pv.plane.len(),
        control_points,
        weights,
    }
}
