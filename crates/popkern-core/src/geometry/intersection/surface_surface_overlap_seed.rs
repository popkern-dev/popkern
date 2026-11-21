//! Overlap rescue seeding: two surfaces overlapping by a little past the
//! tolerance at a near touch meet in a small loop round the deepest point
//! (semi-axes `√(2|g| / λ)` for the gap `g` and the relative curvature's
//! eigenvalues), often smaller than every grid the starting-point search
//! lays: two cones overlapping by ten tolerances at an external touch met
//! in no seed and their common was served empty (W3 near-miss cells).
//!
//! The contact search finds the deepest point as a critical point of the
//! signed gap (a near miss of negative gap). From it, the loop is crossed
//! along rays in the tangent plane: the gap, negative at the point, is
//! followed outward by doubling and then bisected to its zero, and the
//! seed solved there in the plane holding the ray and the normal. Scoped
//! as a rescue: the caller runs it only when a pair emitted nothing.

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::intersection::contact_classify::ContactScale;
use crate::geometry::intersection::contact_discovery::find_surface_contacts;
use crate::geometry::intersection::marching::MarchSample;
use crate::geometry::intersection::marching::marching_newton::march_project_with_iters;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::uv_projection::project_point_to_uv;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

/// Rays per overlap (the march traces the loop from any one seed; the
/// others are consumed by it).
const RAYS: usize = 4;

/// A parameter window.
type Window = ((f64, f64), (f64, f64));

/// A point's parameters.
type Uv = (f64, f64);

/// The rescue rungs of a pair that emitted nothing, from `*rung` on:
/// implicit-distance sign changes on a knot-span grid
/// (`surface_surface_implicit_seed`), then the loops round overlaps.
/// `true` when a rung pushed seeds (march them and come back).
pub(super) fn push_rescue_seeds(
    rung: &mut usize,
    seeds: &mut Vec<(MarchSample, bool)>,
    (a, b): ((&SurfaceGeometry, Window), (&SurfaceGeometry, Window)),
    tolerance: &Tolerance,
) -> bool {
    while seeds.is_empty() && *rung < 2 {
        *rung += 1;
        match *rung {
            1 => super::surface_surface_implicit_seed::push_implicit_grid_seeds(
                seeds, a.0, a.1.0, a.1.1, b.0, b.1.0, b.1.1, tolerance,
            ),
            _ if !contact_search_affordable(a.0, b.0) => {}
            _ => push_overlap_seeds(seeds, a, b, tolerance),
        }
    }
    !seeds.is_empty()
}

/// Whether the contact search is cheap enough to run on every pair of
/// `s1` and `s2` that emits nothing. It is a closest-point solve per grid
/// station each way, iterative on a free-form surface: two Gregory
/// sub-patches meeting in nothing took 8.9 s against a 1 s budget, and a
/// G2-rounded box's own fillet construction (plane × Gregory corner, 72
/// pairs, no seed) went from 316 s to 1319 s in debug. So not between two
/// NURBS surfaces, nor with a Gregory or polyhedral-spline patch; a cheaper
/// search for those pairs is PK-1291.
fn contact_search_affordable(s1: &SurfaceGeometry, s2: &SurfaceGeometry) -> bool {
    let patch = |s: &SurfaceGeometry| {
        matches!(
            s,
            SurfaceGeometry::Gregory(_) | SurfaceGeometry::PolyhedralSpline(_)
        )
    };
    let both_nurbs = matches!(
        (s1, s2),
        (SurfaceGeometry::Nurbs(_), SurfaceGeometry::Nurbs(_))
    );
    !(both_nurbs || patch(s1) || patch(s2))
}

/// Push a seed on the loop round each overlap past the tolerance the
/// contact search finds between `s1` and `s2` over their windows.
fn push_overlap_seeds(
    seeds: &mut Vec<(MarchSample, bool)>,
    (s1, w1): (&SurfaceGeometry, Window),
    (s2, w2): (&SurfaceGeometry, Window),
    tolerance: &Tolerance,
) {
    let extent = window_extent(s1, w1).max(tolerance.linear());
    let scale = ContactScale {
        linear: tolerance.linear(),
        extent,
    };
    let search = find_surface_contacts(s1, w1, s2, w2, scale);
    for miss in search.near_misses {
        if tolerance.within_linear(-miss.gap) || miss.gap > 0.0 {
            continue;
        }
        for seed in loop_seeds(s1, s2, miss.point, extent, tolerance) {
            seeds.push((seed, false));
        }
    }
}

/// The diagonal of `s` sampled at its window's corners and middle.
fn window_extent(s: &SurfaceGeometry, ((u0, u1), (v0, v1)): Window) -> f64 {
    let pts: Vec<Point3> = [0.0, 0.5, 1.0]
        .iter()
        .flat_map(|a| [0.0, 0.5, 1.0].map(move |b| (*a, b)))
        .map(|(a, b)| s.point_at(u0 + (u1 - u0) * a, v0 + (v1 - v0) * b))
        .collect();
    let span = |f: fn(&Point3) -> f64| {
        let (lo, hi) = pts
            .iter()
            .map(f)
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), x| {
                (l.min(x), h.max(x))
            });
        hi - lo
    };
    Vec3::new(span(|p| p.x), span(|p| p.y), span(|p| p.z)).length()
}

/// The signed gap of `s1`'s foot of `p` from `s2`, with both feet.
fn gap_at(s1: &SurfaceGeometry, s2: &SurfaceGeometry, p: Point3) -> Option<(f64, Uv, Uv)> {
    let uv1 = project_point_to_uv(s1, &p)?;
    let q = s1.point_at(uv1.0, uv1.1);
    let uv2 = s2.closest_point(&q)?;
    let n = s2.normal_at(uv2.0, uv2.1).normalized()?;
    Some(((q - s2.point_at(uv2.0, uv2.1)).dot(&n), uv1, uv2))
}

/// Seeds where rays from the overlap's deepest point `at` cross the loop.
fn loop_seeds(
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
    at: Point3,
    extent: f64,
    tolerance: &Tolerance,
) -> Vec<MarchSample> {
    let Some((g0, uv1, _)) = gap_at(s1, s2, at) else {
        return Vec::new();
    };
    let Some(n) = s1.normal_at(uv1.0, uv1.1).normalized() else {
        return Vec::new();
    };
    let (e1, e2) = orthonormal_basis(n);
    let mut out = Vec::new();
    for k in 0..RAYS {
        let a = std::f64::consts::TAU * k as f64 / RAYS as f64;
        let d = e1 * a.pcos() + e2 * a.psin();
        if let Some(seed) = ray_crossing(s1, s2, at, d, g0.signum(), extent, tolerance)
            .and_then(|s| solve_seed(s1, s2, at, s, d, n, tolerance))
        {
            out.push(seed);
        }
    }
    out
}

/// How far along `d` from `at` the gap changes sign from `sign`: doubled
/// from ten tolerances until it does (within `extent`), then bisected.
fn ray_crossing(
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
    at: Point3,
    d: Vec3,
    sign: f64,
    extent: f64,
    tolerance: &Tolerance,
) -> Option<f64> {
    let past = |s: f64| gap_at(s1, s2, at + d * s).map(|(g, _, _)| g * sign <= 0.0);
    let (mut inside, mut outside) = (0.0, 10.0 * tolerance.linear());
    while !past(outside)? {
        inside = outside;
        outside *= 2.0;
        if outside > extent {
            return None;
        }
    }
    for _ in 0..60 {
        let mid = 0.5 * (inside + outside);
        if past(mid)? {
            outside = mid;
        } else {
            inside = mid;
        }
    }
    Some(0.5 * (inside + outside))
}

/// The seed on both surfaces at distance `s` along `d` from `at`, solved
/// in the plane through it holding the ray and the normal `n`.
fn solve_seed(
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
    at: Point3,
    s: f64,
    d: Vec3,
    n: Vec3,
    tolerance: &Tolerance,
) -> Option<MarchSample> {
    let target = at + d * s;
    let (_, uv1, uv2) = gap_at(s1, s2, target)?;
    let plane_normal = n.cross(&d).normalized()?;
    let (seed, _) = march_project_with_iters(
        s1,
        uv1.0,
        uv1.1,
        s2,
        uv2.0,
        uv2.1,
        &target,
        &plane_normal,
        tolerance,
    )?;
    Some(seed)
}

#[cfg(test)]
#[path = "surface_surface_overlap_seed_tests.rs"]
mod tests;
