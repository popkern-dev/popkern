//! Newton projection solvers used by the differential marching algorithm.

use super::marching_linalg::solve_4x4;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

use super::MarchSample;

/// Project onto both surfaces using Newton 4×4.
///
/// Solves the system:
/// - `S1(u1,v1) - S2(u2,v2) = 0`  (3 scalar equations)
/// - `(S1(u1,v1) - target) · march_dir = 0`  (advancing front constraint)
///
/// Returns `None` if the solver fails to converge or the Jacobian is singular.
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn march_project(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
    target: &Point3,
    march_dir: &Vec3,
    tolerance: &Tolerance,
) -> Option<MarchSample> {
    newton_project(
        surface1, u1, v1, surface2, u2, v2, target, march_dir, tolerance,
    )
    .map(|(sample, _)| sample)
}

/// Shared Newton 4×4 core.
///
/// Iterates run UNCLAMPED by design: convergence to interior solutions
/// routinely swings through out-of-domain parameters, and projecting
/// iterates onto the domain box stalls the equality-constrained system
/// (measured on the served pyramid4 strip×cap pair: clamping lost the
/// seed BOTH operand orders previously found). Callers that need
/// in-domain seeds must post-validate the converged sample instead.
#[allow(clippy::too_many_arguments)]
fn newton_project(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
    target: &Point3,
    march_dir: &Vec3,
    tolerance: &Tolerance,
) -> Option<(MarchSample, usize)> {
    let max_iter = 20;
    let mut u1 = u1;
    let mut v1 = v1;
    let mut u2 = u2;
    let mut v2 = v2;

    for iter in 0..max_iter {
        let p1 = surface1.point_at(u1, v1);
        let p2 = surface2.point_at(u2, v2);

        // Rows 0-2 residual: S1 - S2
        let f = p1 - p2;
        // Row 3 residual: (S1 - target) · march_dir
        let f4 = (p1 - *target).dot(march_dir);

        if crate::adaptive_precision::adaptive_dist_sq_lt(
            (p1.x, p1.y, p1.z),
            (p2.x, p2.y, p2.z),
            tolerance.linear() * tolerance.linear(),
        ) && crate::adaptive_precision::adaptive_abs_plane_distance_lt(
            (p1.x, p1.y, p1.z),
            (target.x, target.y, target.z),
            (march_dir.x, march_dir.y, march_dir.z),
            tolerance.linear(),
        ) {
            return Some((
                MarchSample {
                    u1,
                    v1,
                    u2,
                    v2,
                    point: p1,
                },
                iter + 1,
            ));
        }

        let s1u = surface1.partial_u(u1, v1);
        let s1v = surface1.partial_v(u1, v1);
        let s2u = surface2.partial_u(u2, v2);
        let s2v = surface2.partial_v(u2, v2);

        // Build 4×4 Jacobian J such that J * [du1, dv1, du2, dv2]^T = [f.x, f.y, f.z, f4]:
        //
        // Row 0: [S1_u.x,         S1_v.x,         -S2_u.x,  -S2_v.x ]
        // Row 1: [S1_u.y,         S1_v.y,         -S2_u.y,  -S2_v.y ]
        // Row 2: [S1_u.z,         S1_v.z,         -S2_u.z,  -S2_v.z ]
        // Row 3: [T·S1_u,         T·S1_v,          0,        0       ]
        //
        // Derivation of row 3:
        //   F4 = (S1(u1,v1) - target) · T
        //   ∂F4/∂u1 = S1_u · T,  ∂F4/∂v1 = S1_v · T,  ∂F4/∂u2 = 0,  ∂F4/∂v2 = 0
        let j = [
            [s1u.x, s1v.x, -s2u.x, -s2v.x],
            [s1u.y, s1v.y, -s2u.y, -s2v.y],
            [s1u.z, s1v.z, -s2u.z, -s2v.z],
            [march_dir.dot(&s1u), march_dir.dot(&s1v), 0.0, 0.0],
        ];
        let rhs = [f.x, f.y, f.z, f4];

        let delta = solve_4x4(j, rhs)?;

        // WORLD-METRIC trust region — not domain clamping (that
        // variant is measured-falsified above): a healthy Newton
        // step's world-space motion is residual-scale, while a
        // near-singular Jacobian row can fling one variable across
        // the whole chart in a single iteration (measured: a
        // thread-root band × chamfer-cone solve stepped the band's
        // path parameter 16 WINDOWS in one iteration, and whether it
        // happened flipped on a 5e-17 tilt of the cone axis — a
        // knife-edge, not a convergence property). Scaling the step
        // so neither surface moves more than a few residuals leaves
        // out-of-domain excursions and quadratic tail convergence
        // untouched and turns the blowup into a recoverable move.
        let w1 = (s1u * delta[0] + s1v * delta[1]).length();
        let w2 = (s2u * delta[2] + s2v * delta[3]).length();
        let motion = w1.max(w2);
        let residual = (f.x * f.x + f.y * f.y + f.z * f.z).sqrt() + f4.abs();
        let cap = 8.0 * residual + tolerance.linear();
        let scale = if motion > cap { cap / motion } else { 1.0 };

        u1 -= delta[0] * scale;
        v1 -= delta[1] * scale;
        u2 -= delta[2] * scale;
        v2 -= delta[3] * scale;
    }

    None
}

/// Like [`march_project`], but also returns the number of Newton iterations used.
///
/// Returns `None` if the solver fails to converge.
#[allow(clippy::too_many_arguments)]
pub(crate) fn march_project_with_iters(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
    target: &Point3,
    march_dir: &Vec3,
    tolerance: &Tolerance,
) -> Option<(MarchSample, usize)> {
    newton_project(
        surface1, u1, v1, surface2, u2, v2, target, march_dir, tolerance,
    )
}

/// Which chart coordinate a boundary-constrained refinement pins, and
/// to what value (a face iso-boundary in that surface's chart).
#[derive(Debug, Clone, Copy)]
pub struct BoundaryPin {
    /// Pin a coordinate of surface 1's chart (else surface 2's).
    pub on_surface_one: bool,
    /// Pin the u coordinate (else v).
    pub axis_u: bool,
    /// The pinned chart value.
    pub value: f64,
}

/// Newton 4×4 refinement of a surface×surface intersection point onto
/// a chart iso-boundary: rows 0–2 are `S1 − S2`, row 3 pins one chart
/// coordinate to `pin.value`. The pinned coordinate is seeded AT the
/// value, and the pin row is linear, so every iterate satisfies the
/// pin exactly — the system is effectively a well-posed 3-unknown
/// solve for the boundary crossing of the TRUE intersection curve
/// (fitted-curve clip crossings carry the fit's lateral slop instead;
/// measured up to ~1e-3 at grazing junction clusters, PK-313).
///
/// Returns `None` when the Jacobian goes singular (e.g. the curve is
/// tangent to the pinned iso-line) or convergence fails — callers keep
/// their unrefined endpoint in that case.
#[allow(clippy::too_many_arguments)]
pub(crate) fn newton_boundary_crossing(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
    pin: BoundaryPin,
    tolerance: &Tolerance,
) -> Option<MarchSample> {
    let max_iter = 20;
    let mut u1 = u1;
    let mut v1 = v1;
    let mut u2 = u2;
    let mut v2 = v2;
    match (pin.on_surface_one, pin.axis_u) {
        (true, true) => u1 = pin.value,
        (true, false) => v1 = pin.value,
        (false, true) => u2 = pin.value,
        (false, false) => v2 = pin.value,
    }
    for _ in 0..max_iter {
        let p1 = surface1.point_at(u1, v1);
        let p2 = surface2.point_at(u2, v2);
        let f = p1 - p2;
        if crate::adaptive_precision::adaptive_dist_sq_lt(
            (p1.x, p1.y, p1.z),
            (p2.x, p2.y, p2.z),
            tolerance.linear() * tolerance.linear(),
        ) {
            return Some(MarchSample {
                u1,
                v1,
                u2,
                v2,
                point: p1,
            });
        }
        let s1u = surface1.partial_u(u1, v1);
        let s1v = surface1.partial_v(u1, v1);
        let s2u = surface2.partial_u(u2, v2);
        let s2v = surface2.partial_v(u2, v2);
        let mut pin_row = [0.0; 4];
        pin_row[match (pin.on_surface_one, pin.axis_u) {
            (true, true) => 0,
            (true, false) => 1,
            (false, true) => 2,
            (false, false) => 3,
        }] = 1.0;
        let pin_residual = match (pin.on_surface_one, pin.axis_u) {
            (true, true) => u1 - pin.value,
            (true, false) => v1 - pin.value,
            (false, true) => u2 - pin.value,
            (false, false) => v2 - pin.value,
        };
        let j = [
            [s1u.x, s1v.x, -s2u.x, -s2v.x],
            [s1u.y, s1v.y, -s2u.y, -s2v.y],
            [s1u.z, s1v.z, -s2u.z, -s2v.z],
            pin_row,
        ];
        let rhs = [f.x, f.y, f.z, pin_residual];
        let delta = solve_4x4(j, rhs)?;
        u1 -= delta[0];
        v1 -= delta[1];
        u2 -= delta[2];
        v2 -= delta[3];
    }
    None
}

/// Refine a surface×surface intersection point exactly onto a chart
/// iso-boundary (see [`BoundaryPin`]): Newton on `S1 = S2` with one
/// chart coordinate pinned to the boundary value. Used by boolean
/// junction-endpoint refinement — per-pair clip crossings evaluate
/// FITTED record curves and carry the fit's lateral slop (~1e-4–1e-3
/// measured at grazing junction clusters), while this solve lands on
/// the TRUE intersection∩boundary point. Returns the refined chart
/// coordinates of both surfaces and the 3D point, or `None` when the
/// solve is singular (curve tangent to the boundary) or does not
/// converge — callers keep the unrefined endpoint then.
#[allow(clippy::too_many_arguments)]
pub fn refine_boundary_crossing(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
    pin: BoundaryPin,
    tolerance: &Tolerance,
) -> Option<(f64, f64, f64, f64, Point3)> {
    newton_boundary_crossing(surface1, u1, v1, surface2, u2, v2, pin, tolerance)
        .map(|s| (s.u1, s.v1, s.u2, s.v2, s.point))
}

/// A point converged onto a surface-pair intersection: the 3D point
/// plus its uv registration on each surface.
pub type SurfacePairProjection = (Point3, (f64, f64), (f64, f64));

/// Project a 3D point onto the intersection curve of two surfaces.
///
/// Public entry for junction-anchored consumers (the boolean's gap
/// closure, PK-465): given a point NEAR the intersection of two
/// surfaces, converge it onto the curve with the same constrained
/// Newton the SSI march uses, seeded by per-surface projections. The
/// uv guesses go through the knot-span multi-start on free-form
/// charts (a cold projection can walk an adjacent turn's basin).
/// Returns the converged point and both uv registrations, or `None`
/// when either projection or the Newton fails.
pub fn project_point_to_surface_pair(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    near: &Point3,
    tolerance: &Tolerance,
) -> Option<SurfacePairProjection> {
    let project = |surface: &SurfaceGeometry| -> Option<(f64, f64)> {
        let gate = tolerance.linear().max(1e-12);
        let cold = crate::geometry::uv_projection::project_point_to_uv(surface, near);
        if let Some((u, v)) = cold
            // A wrong-turn basin on a coiled chart reads 0.05-0.9 of
            // surface away while still "succeeding" — anything beyond
            // 1e4x linear falls to the knot-seeded multi-start, and
            // the better of the two answers wins.
            && surface.point_at(u, v).distance_to(near) < gate * 1.0e4
        {
            return Some((u, v));
        }
        if let SurfaceGeometry::Nurbs(nurbs) = surface {
            let seeds = nurbs.projection_seeds()?;
            let multi = crate::geometry::surface::closest_point_from_seeds(
                nurbs,
                near,
                &seeds,
                (false, false),
                6,
            );
            return match (cold, multi) {
                (Some(c), Some(m)) => {
                    if surface.point_at(m.0, m.1).distance_to(near)
                        <= surface.point_at(c.0, c.1).distance_to(near)
                    {
                        Some(m)
                    } else {
                        Some(c)
                    }
                }
                (c, m) => m.or(c),
            };
        }
        cold
    };
    let (u1, v1) = project(surface1)?;
    let (u2, v2) = project(surface2)?;
    let n1 = surface1.normal_at(u1, v1);
    let n2 = surface2.normal_at(u2, v2);
    let dir = {
        let cross = n1.cross(&n2);
        if cross.length() > 1e-12 {
            cross
        } else {
            surface1.partial_u(u1, v1)
        }
    };
    if let Some(sample) = march_project(surface1, u1, v1, surface2, u2, v2, near, &dir, tolerance) {
        return Some((sample.point, (sample.u1, sample.v1), (sample.u2, sample.v2)));
    }
    // (continued below)
    // GRAZING fallback: the constrained Newton conditions on n1×n2
    // and exhausts where the surfaces graze — exactly the pairs the
    // junction-anchored callers exist for. Alternating projection
    // (project to A, then to B, repeat) converges linearly onto the
    // intersection neighborhood regardless of crossing angle; each
    // re-projection is seeded at the previous uv so multi-sheet
    // charts keep their basin.
    let mut q = *near;
    let (mut cu1, mut cv1) = (u1, v1);
    let (mut cu2, mut cv2) = (u2, v2);
    for _ in 0..60 {
        let (nu1, nv1) =
            crate::geometry::uv_projection::project_point_to_uv_near(surface1, &q, cu1, cv1)?;
        let p1 = surface1.point_at(nu1, nv1);
        let (nu2, nv2) =
            crate::geometry::uv_projection::project_point_to_uv_near(surface2, &p1, cu2, cv2)?;
        let p2 = surface2.point_at(nu2, nv2);
        let step = p2.distance_to(&q);
        (cu1, cv1, cu2, cv2) = (nu1, nv1, nu2, nv2);
        q = p2;
        if step < tolerance.linear() * 0.01 {
            break;
        }
    }
    let d1 = surface1.point_at(cu1, cv1).distance_to(&q);
    let (fu2, fv2) =
        crate::geometry::uv_projection::project_point_to_uv_near(surface2, &q, cu2, cv2)?;
    let d2 = surface2.point_at(fu2, fv2).distance_to(&q);
    if d1 > tolerance.linear() * 10.0 || d2 > tolerance.linear() * 10.0 {
        return None;
    }
    Some((q, (cu1, cv1), (fu2, fv2)))
}
