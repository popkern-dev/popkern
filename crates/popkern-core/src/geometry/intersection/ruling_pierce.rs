//! The ruling scan shared by the cylinder×quadric pierces: along the
//! cylinder's ruling at angle θ the other surface (a cylinder, a
//! sphere) is met at the roots of a quadratic `a·t² + b·t + c`, its
//! discriminant `D(θ)` deciding the branch structure. Sweeping θ over
//! one period yields every branch by construction: two disjoint loops
//! when `D` stays positive (a full pierce), one loop per positive
//! interval joined at the tangency angles otherwise, and open branches
//! between singular crossings (see `cylinder_cylinder_singular`).

use popkern_math::PortableMath;
use std::f64::consts::TAU;

use super::cylinder_cylinder_singular::{Quad, crossing_branches};
use crate::geometry::Point3;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::surface::SurfaceGeometry;

/// The `(θ, t)` sample sequences of every branch at `density` samples
/// per turn, and whether they are open branches (between crossings)
/// rather than closed loops. `discs` is `D` sampled over one turn.
pub(super) fn ruling_loops(
    quad: Quad<'_>,
    discs: &[f64],
    a_coef: f64,
    crossings: &[f64],
    density: usize,
) -> (Vec<Vec<(f64, f64)>>, bool) {
    let n = discs.len();
    let all_positive = discs.iter().all(|&d| d > 0.0);
    let root = |theta: f64, low: bool| -> f64 {
        let (disc, b, _) = quad(theta);
        let s = disc.max(0.0).sqrt();
        if low {
            (-b - s) / (2.0 * a_coef)
        } else {
            (-b + s) / (2.0 * a_coef)
        }
    };
    let mut loops: Vec<Vec<(f64, f64)>> = Vec::new();
    let crossing_split = if crossings.is_empty() {
        None
    } else {
        crossing_branches(&quad, discs, crossings, a_coef, density)
    };
    let open_branches = crossing_split.is_some();
    if let Some(branches) = crossing_split {
        loops = branches;
    } else if all_positive {
        // Full pierce: two disjoint closed loops, one per root branch.
        for low in [true, false] {
            let pts: Vec<(f64, f64)> = (0..=density)
                .map(|i| {
                    let theta = i as f64 / density as f64 * TAU;
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
                // The exit bracket [last positive, first non-positive]
                // unwrapped TOGETHER past θ_in: an interval running
                // through θ = 0 bracketed [26°, 385°] when each end
                // was unwrapped on its own, and the loop ran the full
                // turn through negative discriminants (PK-782: 9.53
                // off the second cylinder, densified to 32770 points).
                let mut exit_pos = k as f64 / n as f64 * TAU;
                if exit_pos < theta_in {
                    exit_pos += TAU;
                }
                let theta_out = refine(exit_pos + TAU / n as f64, exit_pos);
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
            let m = ((span / TAU * density as f64).ceil() as usize).max(16);
            // Cosine spacing: at the tangency ends the roots move as
            // √(θ − θ_end), which no uniform θ sampling fits — a
            // 50° loop of two r5 cylinders densified to the cap,
            // 32770 points (PK-782). In φ, with θ the cosine map,
            // √disc ∝ sin φ and the loop is smooth.
            let at = |i: usize| {
                let phi = std::f64::consts::PI * i as f64 / m as f64;
                theta_in + span * 0.5 * (1.0 - phi.pcos())
            };
            let mut pts: Vec<(f64, f64)> = Vec::with_capacity(2 * m + 2);
            for i in 0..=m {
                let theta = at(i);
                pts.push((theta, root(theta, true)));
            }
            for i in (0..=m).rev() {
                let theta = at(i);
                pts.push((theta, root(theta, false)));
            }
            // Close back to the start point (same tangency point).
            loops.push(pts);
        }
    }
    (loops, open_branches)
}

/// The largest distance from `curve` (sampled at four points per fit
/// sample) to either surface, measured on the implicit surfaces
/// (cylinders and spheres; any other surface reads 0), and of either
/// pcurve's image from it.
pub(super) fn fit_deviation(
    curve: &CurveOnSurface,
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    samples: usize,
) -> f64 {
    use crate::geometry::curve::CurveEval;
    use crate::geometry::surface::SurfaceEval;
    let off = |surface: &SurfaceGeometry, p: Point3| -> f64 {
        match surface {
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } => {
                let Some(a) = axis.normalized() else {
                    return 0.0;
                };
                let rel = p - *origin;
                ((rel - a * rel.dot(&a)).length() - radius).abs()
            }
            SurfaceGeometry::Sphere { center, radius, .. } => {
                (p.distance_to(center) - radius).abs()
            }
            _ => 0.0,
        }
    };
    // The pcurves' images too (the fit gives the 3D curve and both
    // pcurves one parameter): near a sphere's pole a loop's chart u
    // turns fast while θ steps evenly, and a fit densified on the 3D
    // curve alone left its sphere pcurve far off, the boolean keeping
    // the wrong side of it.
    let image_off =
        |pcurve: &crate::geometry::NurbsCurve, surface: &SurfaceGeometry, t: f64, p: Point3| {
            let (s0, s1) = pcurve.parameter_range();
            let (c0, c1) = curve.curve_3d.parameter_range();
            let s = s0 + (s1 - s0) * (t - c0) / (c1 - c0);
            let q = pcurve.point_at(s);
            surface.point_at(q.x, q.y).distance_to(&p)
        };
    let (t0, t1) = curve.curve_3d.parameter_range();
    let k = 4 * samples.max(1);
    (0..=k)
        .map(|i| {
            let t = t0 + (t1 - t0) * i as f64 / k as f64;
            let p = curve.curve_3d.point_at(t);
            off(surface1, p)
                .max(off(surface2, p))
                .max(image_off(&curve.curve_uv1, surface1, t, p))
                .max(image_off(&curve.curve_uv2, surface2, t, p))
        })
        .fold(0.0, f64::max)
}
