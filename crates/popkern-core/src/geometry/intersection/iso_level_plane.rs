//! Exact intersections of a plane with a NURBS surface whose control
//! rows each stand at one height above the plane (PK-1209).
//!
//! When every control point `P_ij` of one row `j` (one index of the
//! other direction held) lies at the same signed height `h_j` above a
//! plane and the weights separate (`w_ij = a_i b_j`), the surface's
//! height above the plane at `(u, v)` is
//! `Σ_j N_j(v) b_j h_j / Σ_j N_j(v) b_j`, a function of `v` alone: the
//! plane meets the surface along the iso-curves `v = v*` at the roots
//! of that one rational function. A NURBS ellipsoid (a sphere's net
//! scaled along each axis) cut square to an axis is the everyday case:
//! the marcher traced its section, an iso-curve closed round the seam,
//! as records whose ends did not close there, and the complex declined
//! on loop gaps at most cut heights (placement luck: z = 0, ±0.5 and
//! −0.7 declined, −0.3 served).
//!
//! The test reads the control net, so it is exact by the tensor-product
//! form, not sampled: the heights are equal within the linear tolerance
//! and the weight ratios within rounding. Each curve is the surface's
//! own iso-curve, on its chart the line `v = v*`; on the plane its exact
//! affine image. A root at a collapsed row (a pole, the curve a point)
//! and a touch (no sign change) emit nothing, as every arm does:
//! contact discovery types them. A stretch of rows lying in the plane
//! is a coincidence, left to the coincidence detector. Anything else
//! keeps the marching path.

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::circle_of_revolution::pcurve_on;
use crate::geometry::intersection::revolved_nurbs::{domain, transpose};
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

/// A parameter window `((u0, u1), (v0, v1))`.
type Window = ((f64, f64), (f64, f64));

/// The plane × level-row NURBS answer, in argument order, or `None`
/// when the pair is not of that class.
pub(crate) fn try_iso_level_plane(
    surface1: &SurfaceGeometry,
    w1: Window,
    surface2: &SurfaceGeometry,
    w2: Window,
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let (nurbs, wn, plane, wp, swapped) = match (surface1, surface2) {
        (SurfaceGeometry::Nurbs(n), p @ SurfaceGeometry::Plane { .. }) => (n, w1, p, w2, false),
        (p @ SurfaceGeometry::Plane { .. }, SurfaceGeometry::Nurbs(n)) => (n, w2, p, w1, true),
        _ => return None,
    };
    let SurfaceGeometry::Plane { origin, normal, .. } = plane else {
        return None;
    };
    let normal = normal.normalized()?;
    let height = |p: Point3| (p - *origin).dot(&normal);
    // Height a function of u: the curves are iso-u. A function of v is
    // read through the transpose, so the curves are always the view's
    // `iso_curve_at_u`.
    let transposed = transpose(nurbs);
    let (view, level_window, curves_in_v) = if level_rows(nurbs, &height, tolerance) {
        (nurbs, wn.0, false)
    } else if level_rows(&transposed, &height, tolerance) {
        (&transposed, wn.1, true)
    } else {
        return None;
    };
    let roots = level_roots(view, level_window, &height, tolerance)?;
    let (v0, v1) = domain(&view.knots_v, view.degree_v);
    let mut out = Vec::with_capacity(roots.len());
    for u in roots {
        let curve_3d = view.iso_curve_at_u(u);
        // A collapsed row (a pole) is a point, discovery's to type.
        let span = (0..=8)
            .map(|k| curve_3d.point_at(v0 + (v1 - v0) * f64::from(k) / 8.0))
            .collect::<Vec<_>>();
        if span
            .iter()
            .all(|p| p.distance_to(&span[0]) <= tolerance.linear())
        {
            continue;
        }
        let ends = [(u, v0), (u, v1)].map(|(p, a)| {
            if curves_in_v {
                Point3::new(a, p, 0.0)
            } else {
                Point3::new(p, a, 0.0)
            }
        });
        let on_nurbs = NurbsCurve {
            degree: 1,
            knots: vec![v0, v0, v1, v1],
            control_points: ends.to_vec(),
            weights: vec![1.0, 1.0],
        };
        let on_plane = pcurve_on(plane, &curve_3d, wp.0, wp.1, tolerance)?;
        let is_closed =
            curve_3d.point_at(v0).distance_to(&curve_3d.point_at(v1)) <= tolerance.linear();
        let (curve_uv1, curve_uv2) = if swapped {
            (on_plane, on_nurbs)
        } else {
            (on_nurbs, on_plane)
        };
        out.push(SurfaceSurfaceIntersection::Curve(CurveOnSurface {
            curve_3d,
            curve_uv1,
            curve_uv2,
            is_closed,
        }));
    }
    Some(out)
}

/// Whether each row `i` of `s` (its u index held) stands at one height
/// above the plane, the weights separable and the heights differing
/// between rows: then the surface's height above the plane is a
/// function of `u` alone.
fn level_rows(s: &NurbsSurface, height: &impl Fn(Point3) -> f64, tolerance: &Tolerance) -> bool {
    let lin = tolerance.linear();
    if s.num_u < 2 || s.num_v < 2 {
        return false;
    }
    let at = |i: usize, j: usize| s.control_points[i * s.num_v + j];
    let wt = |i: usize, j: usize| s.weights[i * s.num_v + j];
    // Height varies across i somewhere (else the surface is the plane's
    // parallel, no curve or a coincidence: not this arm's).
    let h0: Vec<f64> = (0..s.num_u).map(|i| height(at(i, 0))).collect();
    for (i, &row_height) in h0.iter().enumerate() {
        for j in 0..s.num_v {
            if (height(at(i, j)) - row_height).abs() > lin {
                return false;
            }
            let ratio = wt(i, j) / wt(i, 0);
            let first = wt(0, j) / wt(0, 0);
            if !ratio.is_finite() || (ratio - first).abs() > 1e-9 * first.abs().max(1.0) {
                return false;
            }
        }
    }
    let (lo, hi) = h0
        .iter()
        .fold((f64::MAX, f64::MIN), |(a, b), &h| (a.min(h), b.max(h)));
    hi - lo > lin
}

/// The `u` parameters, within `u_window`, where the height above the
/// plane (a function of `u` alone) crosses zero; `None` when a stretch
/// of it lies in the plane (a coincidence).
fn level_roots(
    s: &NurbsSurface,
    u_window: (f64, f64),
    height: &impl Fn(Point3) -> f64,
    tolerance: &Tolerance,
) -> Option<Vec<f64>> {
    let lin = tolerance.linear();
    let (k0, k1) = domain(&s.knots_u, s.degree_u);
    let (lo, hi) = (
        u_window.0.min(u_window.1).max(k0),
        u_window.0.max(u_window.1).min(k1),
    );
    if hi <= lo {
        return Some(Vec::new());
    }
    let v_mid = {
        let (v0, v1) = domain(&s.knots_v, s.degree_v);
        0.5 * (v0 + v1)
    };
    let g = |u: f64| height(s.point_at(u, v_mid));
    let n = (32 * (s.num_u - s.degree_u).max(1)).clamp(64, 4096);
    let us: Vec<f64> = (0..=n)
        .map(|k| lo + (hi - lo) * k as f64 / n as f64)
        .collect();
    let gs: Vec<f64> = us.iter().map(|&u| g(u)).collect();
    if gs.windows(3).any(|w| w.iter().all(|x| x.abs() <= lin)) {
        return None;
    }
    let mut roots = Vec::new();
    for k in 0..n {
        let (ga, gb) = (gs[k], gs[k + 1]);
        if ga == 0.0 {
            roots.push(us[k]);
            continue;
        }
        if ga * gb > 0.0 {
            continue;
        }
        let (mut a, mut b, mut g_a) = (us[k], us[k + 1], ga);
        for _ in 0..100 {
            let m = 0.5 * (a + b);
            let gm = g(m);
            if gm * g_a <= 0.0 {
                b = m;
            } else {
                a = m;
                g_a = gm;
            }
        }
        roots.push(0.5 * (a + b));
    }
    if gs[n] == 0.0 {
        roots.push(us[n]);
    }
    roots.dedup_by(|a, b| (*a - *b).abs() <= 1e-12 * (hi - lo).max(1.0));
    // A root on the domain's boundary row is the plane lying along the
    // surface's edge (a sheet's boundary flush on a block's side): an
    // edge coincidence, not a section, kept for the marcher and the
    // coincidence detector as before.
    let edge_reach = 1e-9 * (k1 - k0).max(1.0);
    if roots
        .iter()
        .any(|&u| (u - k0).abs() <= edge_reach || (u - k1).abs() <= edge_reach)
    {
        return None;
    }
    Some(roots)
}

#[cfg(test)]
#[path = "iso_level_plane_tests.rs"]
mod tests;
