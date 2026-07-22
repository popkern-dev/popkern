//! Analytic cylinder×cylinder intersection (non-parallel axes).
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

use std::f64::consts::TAU;

use crate::error::CoreError;
use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::intersection::marching::MarchSample;
use crate::geometry::intersection::surface_surface::{
    SurfaceSurfaceIntersection, fit_march_samples,
};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Samples per full pierce loop. The samples are EXACT surface
/// points; the cubic fit through them carries interpolation error
/// O(h⁴), far below linear tolerance at this density.
const LOOP_SAMPLES: usize = 512;

/// Try the analytic cylinder×cylinder path. Returns `None` for
/// non-cylinder pairs and for (near-)parallel axes — parallel
/// cylinders meet in lines or coincide, which the line/coincident
/// detectors own.
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
        return None;
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
        let origin = so + sp1 * (sr * theta.cos()) + sp2 * (sr * theta.sin());
        let dp = origin - to;
        let dp_perp = dp - ta * dp.dot(&ta);
        let b = 2.0 * dp_perp.dot(&d_perp);
        let c = dp_perp.dot(&dp_perp) - tr * tr;
        (b * b - 4.0 * a_coef * c, b, c)
    };
    let point_at = |theta: f64, t: f64| -> Point3 {
        so + sp1 * (sr * theta.cos()) + sp2 * (sr * theta.sin()) + sa * t
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
        let mut u = rel.dot(&p2).atan2(rel.dot(&p1));
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::CurveEval as _;

    #[test]
    fn cross_drill_pierce_loops_exact() {
        let body = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 16.0,
        };
        let drill = SurfaceGeometry::Cylinder {
            origin: Point3::new(-40.0, 0.0, 20.0),
            axis: Vec3::new(1.0, 0.0, 0.0),
            radius: 6.0,
        };
        let tol = Tolerance::default();
        let curves = try_cylinder_cylinder(&body, (0.0, 40.0), &drill, (0.0, 80.0), &tol)
            .expect("cylinder pair must take the analytic path");
        assert_eq!(curves.len(), 2, "full pierce must yield two loops");
        for c in &curves {
            let SurfaceSurfaceIntersection::Curve(cos) = c else {
                panic!("expected curve");
            };
            let (t0, t1) = cos.curve_3d.parameter_range();
            let mut worst = 0.0_f64;
            for i in 0..=200 {
                let t = t0 + (t1 - t0) * f64::from(i) / 200.0;
                let p = cos.curve_3d.point_at(t);
                let d_body = ((p.x * p.x + p.y * p.y).sqrt() - 16.0).abs();
                let d_drill = ((p.y * p.y + (p.z - 20.0) * (p.z - 20.0)).sqrt() - 6.0).abs();
                worst = worst.max(d_body.max(d_drill));
            }
            assert!(
                worst < 1e-6,
                "pierce loop deviates {worst:.3e} from the surfaces"
            );
        }
    }
}
