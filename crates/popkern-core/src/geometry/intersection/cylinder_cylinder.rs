//! Analytic cylinder×cylinder intersection.
//!
//! The marching solver fragments the everyday cross-drill class
//! (perpendicular bore through a cylinder wall) into overlapping
//! partial traces and can miss whole branches (measured: the
//! cross-drilled boss produced 20 raw records consolidating to 7
//! overlapping pieces of ONE pierce loop, the second loop absent).
//! This path constructs the curve exactly: every intersection point
//! lies on a ruling of the smaller cylinder, and each ruling meets
//! the other cylinder in the roots of a plain quadratic. Sweeping the
//! ruling angle θ over one period yields the complete branch
//! structure by construction — two disjoint loops when the quadratic
//! stays positive (full pierce), tangent-joined loops per
//! positive-discriminant interval otherwise.

use popkern_math::PortableMath;
use std::f64::consts::TAU;

use crate::error::CoreError;
use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::marching::MarchSample;
use crate::geometry::intersection::surface_surface::{
    SurfaceSurfaceIntersection, fit_march_samples,
};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Samples per full pierce loop. The samples are EXACT surface
/// points; the cubic fit through them carries interpolation error
/// O(h⁴), far below linear tolerance at this density.
const LOOP_SAMPLES: usize = 512;

/// Try the analytic cylinder×cylinder path. Returns `None` only for
/// non-cylinder pairs and for near-coincident parallel pairs (the
/// coincident detector owns those); parallel offset axes take the
/// exact line arm below.
pub(crate) fn try_cylinder_cylinder(
    surface1: &SurfaceGeometry,
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let (
        SurfaceGeometry::Cylinder {
            origin: o1,
            axis: a1,
            radius: r1,
        },
        SurfaceGeometry::Cylinder {
            origin: o2,
            axis: a2,
            radius: r2,
        },
    ) = (surface1, surface2)
    else {
        return None;
    };
    let a1n = a1.normalized()?;
    let a2n = a2.normalized()?;
    if a1n.cross(&a2n).length() < tolerance.angular() {
        // POPKERN_ANALYTIC_PARCYL_SSI=0 pins the marching path for
        // A/B measurement and containment (PK-239 precedent; a
        // per-call read like the torus switch).
        let arm_enabled = std::env::var_os("POPKERN_ANALYTIC_PARCYL_SSI")
            .map(|v| v != "0")
            .unwrap_or(true);
        if !arm_enabled {
            return None;
        }
        return try_parallel_axes(
            *o1, *a1, a1n, *r1, v1_range, *o2, *a2, a2n, *r2, v2_range, tolerance,
        );
    }
    // Parameterize on the SMALLER cylinder (its rulings cross the
    // larger transversally over the whole pierce).
    let (so, sa, sr, to, ta, tr) = if r2 < r1 {
        (*o2, a2n, *r2, *o1, a1n, *r1)
    } else {
        (*o1, a1n, *r1, *o2, a2n, *r2)
    };
    let (sp1, sp2) = orthonormal_basis(sa);

    // Ruling direction is constant: precompute its component
    // perpendicular to the target axis.
    let d_dot_a = sa.dot(&ta);
    let d_perp = sa - ta * d_dot_a;
    let a_coef = d_perp.dot(&d_perp);
    if a_coef < tolerance.angular() * tolerance.angular() {
        return None; // effectively parallel
    }

    // Quadratic along the ruling at angle θ:
    //   |dp_perp(θ) + t·d_perp|² = tr²
    let quad = |theta: f64| -> (f64, f64, f64) {
        let origin = so + sp1 * (sr * theta.pcos()) + sp2 * (sr * theta.psin());
        let dp = origin - to;
        let dp_perp = dp - ta * dp.dot(&ta);
        let b = 2.0 * dp_perp.dot(&d_perp);
        let c = dp_perp.dot(&dp_perp) - tr * tr;
        (b * b - 4.0 * a_coef * c, b, c)
    };
    let point_at = |theta: f64, t: f64| -> Point3 {
        so + sp1 * (sr * theta.pcos()) + sp2 * (sr * theta.psin()) + sa * t
    };

    // Scan the discriminant over one period.
    let n = LOOP_SAMPLES;
    let discs: Vec<f64> = (0..n).map(|i| quad(i as f64 / n as f64 * TAU).0).collect();
    let all_positive = discs.iter().all(|&d| d > 0.0);
    let all_negative = discs.iter().all(|&d| d <= 0.0);
    if all_negative {
        return Some(Vec::new()); // no intersection (or grazing below resolution)
    }

    let root = |theta: f64, low: bool| -> f64 {
        let (disc, b, _) = quad(theta);
        let s = disc.max(0.0).sqrt();
        if low {
            (-b - s) / (2.0 * a_coef)
        } else {
            (-b + s) / (2.0 * a_coef)
        }
    };

    let mut loops: Vec<Vec<(f64, f64)>> = Vec::new(); // (θ, t) sequences
    if all_positive {
        // Full pierce: two disjoint closed loops, one per root branch.
        for low in [true, false] {
            let pts: Vec<(f64, f64)> = (0..=n)
                .map(|i| {
                    let theta = i as f64 / n as f64 * TAU;
                    (theta, root(theta, low))
                })
                .collect();
            loops.push(pts);
        }
    } else {
        // Partial pierce: one closed loop per positive-discriminant
        // interval — out along the low root, back along the high
        // root, joined at the tangency angles (disc = 0), refined by
        // bisection.
        let refine = |mut neg: f64, mut pos: f64| -> f64 {
            for _ in 0..60 {
                let mid = 0.5 * (neg + pos);
                if quad(mid).0 > 0.0 {
                    pos = mid;
                } else {
                    neg = mid;
                }
            }
            0.5 * (neg + pos)
        };
        // Find sign-change boundaries walking the scan circularly.
        let mut intervals: Vec<(f64, f64)> = Vec::new();
        for i in 0..n {
            let j = (i + 1) % n;
            if discs[i] <= 0.0 && discs[j] > 0.0 {
                // entering positive: find matching exit
                let theta_in = refine(i as f64 / n as f64 * TAU, j as f64 / n as f64 * TAU);
                let mut k = j;
                while discs[(k + 1) % n] > 0.0 {
                    k = (k + 1) % n;
                    if k == i {
                        break;
                    }
                }
                let exit_neg = ((k + 1) % n) as f64 / n as f64 * TAU;
                let exit_pos = k as f64 / n as f64 * TAU;
                let theta_out = refine(
                    if exit_neg > exit_pos {
                        exit_neg
                    } else {
                        exit_neg + TAU
                    },
                    if exit_pos >= theta_in {
                        exit_pos
                    } else {
                        exit_pos + TAU
                    },
                );
                let theta_out = if theta_out > theta_in {
                    theta_out
                } else {
                    theta_out + TAU
                };
                intervals.push((theta_in, theta_out));
            }
        }
        for (theta_in, theta_out) in intervals {
            let span = theta_out - theta_in;
            let m = ((span / TAU * n as f64).ceil() as usize).max(16);
            let mut pts: Vec<(f64, f64)> = Vec::with_capacity(2 * m + 2);
            for i in 0..=m {
                let theta = theta_in + span * i as f64 / m as f64;
                pts.push((theta, root(theta, true)));
            }
            for i in (0..=m).rev() {
                let theta = theta_in + span * i as f64 / m as f64;
                pts.push((theta, root(theta, false)));
            }
            // Close back to the start point (same tangency point).
            loops.push(pts);
        }
    }

    // Chart parameters for a point on each cylinder, continuously
    // unwrapped along the sample sequence.
    let chart = |surface: &SurfaceGeometry, p: Point3, prev_u: Option<f64>| -> (f64, f64) {
        let SurfaceGeometry::Cylinder { origin, axis, .. } = surface else {
            return (0.0, 0.0);
        };
        let an = axis.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
        let (p1, p2) = orthonormal_basis(an);
        let rel = p - *origin;
        // point_at is origin + axis·v (raw axis): invert with the
        // squared length so non-unit stored axes stay exact.
        let v = rel.dot(axis) / axis.length_squared();
        let mut u = rel.dot(&p2).patan2(rel.dot(&p1));
        if let Some(pu) = prev_u {
            while u - pu > std::f64::consts::PI {
                u -= TAU;
            }
            while pu - u > std::f64::consts::PI {
                u += TAU;
            }
        }
        (u, v)
    };

    let pad1 = (v1_range.1 - v1_range.0).abs() * 0.5 + tolerance.linear();
    let pad2 = (v2_range.1 - v2_range.0).abs() * 0.5 + tolerance.linear();
    let mut results = Vec::new();
    for pts in loops {
        if pts.len() < 4 {
            continue;
        }
        let mut samples: Vec<MarchSample> = Vec::with_capacity(pts.len());
        let (mut u1_prev, mut u2_prev) = (None, None);
        for &(theta, t) in &pts {
            let p = point_at(theta, t);
            let (u1, v1) = chart(surface1, p, u1_prev);
            let (u2, v2) = chart(surface2, p, u2_prev);
            u1_prev = Some(u1);
            u2_prev = Some(u2);
            samples.push(MarchSample {
                u1,
                v1,
                u2,
                v2,
                point: p,
            });
        }
        // Skip loops entirely outside both requested v-windows (the
        // caller's clip stage owns exact trimming; this only avoids
        // fitting curves with no chance of surviving it).
        let v1_hit = samples
            .iter()
            .any(|s| s.v1 > v1_range.0 - pad1 && s.v1 < v1_range.1 + pad1);
        let v2_hit = samples
            .iter()
            .any(|s| s.v2 > v2_range.0 - pad2 && s.v2 < v2_range.1 + pad2);
        if !v1_hit || !v2_hit {
            continue;
        }
        match fit_march_samples(&samples, true) {
            Ok(curve) => results.push(SurfaceSurfaceIntersection::Curve(curve)),
            Err(CoreError::PreconditionFailed(_)) => continue,
            Err(_) => continue,
        }
    }
    Some(results)
}

/// Parallel-axis cylinder×cylinder: the cross-section is a plain 2D
/// circle×circle problem, so the surfaces meet in zero, one
/// (tangency) or two lines parallel to the shared axis direction.
/// The marching path solved this class numerically — the everyday
/// boss-beside-fillet booleans of a shelled box spent ~400 ms per
/// union fitting NURBS through straight lines (measured, PK-427) —
/// and the fitted spellings then miss every exact-line downstream
/// serve. Tangency is silent (the PK-252 contract: a tangent contact
/// is not a crossing curve) and near-coincident pairs return `None`
/// so the coincident detector keeps owning them.
#[allow(clippy::too_many_arguments)]
fn try_parallel_axes(
    o1: Point3,
    a1: Vec3,
    a1n: Vec3,
    r1: f64,
    v1_range: (f64, f64),
    o2: Point3,
    a2: Vec3,
    a2n: Vec3,
    r2: f64,
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let lin = tolerance.linear();
    // Axis separation measured in the cross-section plane of axis 1.
    let delta = o2 - o1;
    let w = delta - a1n * delta.dot(&a1n);
    let d = w.length();
    if d < lin && (r1 - r2).abs() < lin {
        return None; // coincident surfaces — the coincident detector owns them
    }
    if d > r1 + r2 - lin {
        // Separated, or externally tangent (silent).
        return Some(Vec::new());
    }
    if d < (r1 - r2).abs() + lin {
        // Nested, or internally tangent (silent).
        return Some(Vec::new());
    }

    // 2D circle×circle roots in the (e, perp) frame of axis 1.
    let e = w * (1.0 / d);
    let perp = a1n.cross(&e);
    let along = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let h_sq = r1 * r1 - along * along;
    if h_sq <= 0.0 {
        return Some(Vec::new());
    }
    let h = h_sq.sqrt();
    let mid = o1 + e * along;

    // Span both callers' axial windows with headroom, in arclength
    // along the shared direction (the clip stage owns exact
    // trimming; matching `case_parallel`'s pad convention).
    let len1 = a1.length().max(1e-12);
    let len2 = a2.length().max(1e-12);
    let start2 = delta.dot(&a1n); // arclength offset of o2's plane
    let (mut t_lo, mut t_hi) = (v1_range.0 * len1, v1_range.1 * len1);
    for t2 in [
        start2 + v2_range.0 * len2 * a2n.dot(&a1n),
        start2 + v2_range.1 * len2 * a2n.dot(&a1n),
    ] {
        t_lo = t_lo.min(t2);
        t_hi = t_hi.max(t2);
    }
    let pad = (t_hi - t_lo).abs().max(1.0);
    let (t_lo, t_hi) = (t_lo - pad, t_hi + pad);

    let (p1a, p1b) = orthonormal_basis(a1n);
    let (p2a, p2b) = orthonormal_basis(a2n);
    let mut results = Vec::with_capacity(2);
    for sign in [1.0, -1.0] {
        let q = mid + perp * (sign * h);
        let radial1 = q - o1;
        let u1 = norm_angle(f64::patan2(radial1.dot(&p1b), radial1.dot(&p1a)));
        let radial2 = (q - o2) - a2n * (q - o2).dot(&a2n);
        let u2 = norm_angle(f64::patan2(radial2.dot(&p2b), radial2.dot(&p2a)));
        let a = q + a1n * t_lo;
        let b = q + a1n * t_hi;
        // Chart v satisfies point = origin + … + axis·v, so
        // v = (p − o)·axis / |axis|².
        let v_of = |p: Point3, o: Point3, axis: Vec3, len: f64| (p - o).dot(&axis) / (len * len);
        let curve_3d = nurbs_line(a, b);
        let curve_uv1 = nurbs_line(
            Point3::new(u1, v_of(a, o1, a1, len1), 0.0),
            Point3::new(u1, v_of(b, o1, a1, len1), 0.0),
        );
        let curve_uv2 = nurbs_line(
            Point3::new(u2, v_of(a, o2, a2, len2), 0.0),
            Point3::new(u2, v_of(b, o2, a2, len2), 0.0),
        );
        results.push(SurfaceSurfaceIntersection::Curve(CurveOnSurface {
            curve_3d,
            curve_uv1,
            curve_uv2,
            is_closed: false,
        }));
    }
    Some(results)
}

/// Normalize an angle into the cylinder chart's [0, 2π) branch.
fn norm_angle(a: f64) -> f64 {
    let mut r = a % TAU;
    if r < 0.0 {
        r += TAU;
    }
    r
}

/// Degree-1 NURBS line whose parameter is arclength.
fn nurbs_line(p0: Point3, p1: Point3) -> NurbsCurve {
    let len = p0.distance_to(&p1).max(1e-12);
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, len, len],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    }
}

#[cfg(test)]
#[path = "cylinder_cylinder_tests.rs"]
mod tests;
