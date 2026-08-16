//! Sampled-analytic solver for general-position plane×torus
//! intersections — the quartic spiric sections (PK-279).
//!
//! A generic spiric section has genus one and admits no exact rational
//! parameterization, so unlike the circle classes it cannot be an
//! exact NURBS. What CAN be exact is every sample: for a fixed axis
//! angle `u` the plane cuts the tube circle where
//! `K(u)·cos(v − φ(u)) = A(u)` with
//!
//! - `m(u) = n1·cos u + n2·sin u` (plane normal's radial component),
//! - `A(u) = d0 − R·m(u)`, `K(u) = r·√(m(u)² + na²)`,
//! - `φ(u) = atan2(na, m(u))`,
//!
//! giving closed-form tube angles `v±(u) = φ ± acos(A/K)` wherever the
//! feasibility function `F(u) = A² − K²` is non-positive. Branch
//! enumeration is therefore COMPLETE by construction: the maximal
//! `F ≤ 0` arcs of the axis circle are found by dense sampling plus
//! Newton–bisection root refinement (the `line_torus` pattern), each
//! arc closes into exactly one loop (the `v+` and `v−` branches meet
//! at the arc ends where `|A| = K`), and a root-free negative `F`
//! yields the two-loop section. The marching path, by contrast, finds
//! only the branches its seeding happens to reach.
//!
//! Loops are emitted through [`fit_march_samples`] — the same fitting
//! the marching path uses — with samples densified until the fitted
//! curve lies on BOTH surfaces within half the linear tolerance.
//! Failure to converge declines the arm (`None`) and the caller keeps
//! the marching fallback; loops smaller than 10× linear tolerance are
//! dropped as grazing-contact blobs per the tangency-silent contract.

use crate::geometry::Point3;
use crate::geometry::intersection::marching::MarchSample;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::intersection::surface_surface::fit_march_samples;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;
use std::f64::consts::TAU;

use super::plane_torus_spiric_frame::SpiricFrame;

/// Feasibility scan resolution around the axis circle. 2048 matches
/// the `line_torus` quartic bracketing density; `F` is a degree-2
/// trigonometric polynomial, so consecutive same-sign samples at this
/// density cannot hide a sign-change pair.
pub(super) const SCAN_SAMPLES: usize = 2048;

/// Per-branch sample counts tried in order until the fitted loop
/// passes the serving verification gate. The ladder starts at 513:
/// measured on the R=2/r=0.5 fixtures, 129 samples fit to ~4e-6
/// (always failing the half-tolerance gate) while 513 reaches ~1e-8,
/// so a 129 rung only wasted a fit per loop. 2049 is the last resort
/// for small-feature pcurve accuracy (fit control points cap at 200).
const BRANCH_SAMPLES: [usize; 3] = [513, 1025, 2049];

/// Solve the general-position plane×torus intersection. `swapped` =
/// the torus is operand 1. Returns `None` only when fitting or
/// numeric guards decline — the caller falls back to marching.
pub(super) fn try_spiric_sections(
    plane: &SurfaceGeometry,
    torus: &SurfaceGeometry,
    torus_windows: ((f64, f64), (f64, f64)),
    swapped: bool,
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let SurfaceGeometry::Plane {
        origin: po,
        normal,
        u_axis,
        v_axis,
    } = plane
    else {
        return None;
    };
    let SurfaceGeometry::Torus {
        center,
        axis,
        major_radius,
        minor_radius,
    } = torus
    else {
        return None;
    };
    let n_hat = normal.normalized()?;
    let axis_n = axis.normalized()?;
    let frame = SpiricFrame::new(
        *po,
        n_hat,
        *u_axis,
        *v_axis,
        *center,
        axis_n,
        *major_radius,
        *minor_radius,
    )?;

    // Locate the F sign structure around the axis circle.
    let roots = frame.feasibility_roots();
    let arcs = if roots.is_empty() {
        if frame.f_at(0.0) > 0.0 {
            return Some(Vec::new()); // plane misses the torus
        }
        Vec::new() // negative everywhere: two full-period loops
    } else {
        frame.active_arcs(&roots)
    };

    let mut curves = Vec::new();
    let mut emit = |mut samples: Vec<MarchSample>| -> Option<bool> {
        // Land the torus chart chain in the CALLER's branch: the
        // marching path wraps samples into the caller windows and
        // downstream UV clipping matches against face trim domains in
        // that branch (see `align_pcurve_to_windows`). One constant
        // whole-period shift per axis preserves chain continuity.
        align_chain_to_windows(&mut samples, torus_windows, swapped);
        // Grazing-contact blob: drop, per the tangency-silent contract.
        if loop_extent(&samples) <= 10.0 * tolerance.linear() {
            return Some(false);
        }
        let curve = fit_march_samples(&samples, true).ok()?;
        if !frame.fitted_curve_serves(&curve, swapped, 0.5 * tolerance.linear().max(1.0e-9)) {
            return None;
        }
        curves.push(SurfaceSurfaceIntersection::Curve(curve));
        Some(true)
    };

    if roots.is_empty() {
        // Two disjoint loops, one per tube-angle branch.
        for branch in [1.0f64, -1.0] {
            let mut accepted = false;
            for &n in &BRANCH_SAMPLES {
                let samples = frame.full_period_branch(branch, n, swapped);
                match emit(samples) {
                    Some(true) | Some(false) => {
                        accepted = true;
                        break;
                    }
                    None => continue, // densify and retry
                }
            }
            if !accepted {
                return None;
            }
        }
    } else {
        for &(ua, ub) in &arcs {
            let mut accepted = false;
            for &n in &BRANCH_SAMPLES {
                let samples = frame.arc_loop(ua, ub, n, swapped);
                match emit(samples) {
                    Some(true) | Some(false) => {
                        accepted = true;
                        break;
                    }
                    None => continue,
                }
            }
            if !accepted {
                return None;
            }
        }
    }
    Some(curves)
}

/// Shift the torus-side (u, v) of a sample chain by whole periods so
/// the chain's centre lands as close as possible to the caller's
/// window centres.
fn align_chain_to_windows(
    samples: &mut [MarchSample],
    torus_windows: ((f64, f64), (f64, f64)),
    swapped: bool,
) {
    for axis in [0usize, 1] {
        let window = if axis == 0 {
            torus_windows.0
        } else {
            torus_windows.1
        };
        let value = |s: &MarchSample| match (swapped, axis) {
            (true, 0) => s.u1,
            (true, 1) => s.v1,
            (false, 0) => s.u2,
            (false, 1) => s.v2,
            _ => unreachable!(),
        };
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for s in samples.iter() {
            lo = lo.min(value(s));
            hi = hi.max(value(s));
        }
        if !lo.is_finite() || !hi.is_finite() {
            continue;
        }
        let k = ((0.5 * (window.0 + window.1) - 0.5 * (lo + hi)) / TAU).round();
        if k != 0.0 {
            for s in samples.iter_mut() {
                match (swapped, axis) {
                    (true, 0) => s.u1 += k * TAU,
                    (true, 1) => s.v1 += k * TAU,
                    (false, 0) => s.u2 += k * TAU,
                    (false, 1) => s.v2 += k * TAU,
                    _ => unreachable!(),
                }
            }
        }
    }
}

/// Bounding-box diagonal of a sample chain.
fn loop_extent(samples: &[MarchSample]) -> f64 {
    let mut lo = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for s in samples {
        lo.x = lo.x.min(s.point.x);
        lo.y = lo.y.min(s.point.y);
        lo.z = lo.z.min(s.point.z);
        hi.x = hi.x.max(s.point.x);
        hi.y = hi.y.max(s.point.y);
        hi.z = hi.z.max(s.point.z);
    }
    hi.distance_to(&lo)
}

#[cfg(test)]
#[path = "plane_torus_spiric_tests.rs"]
mod tests;
