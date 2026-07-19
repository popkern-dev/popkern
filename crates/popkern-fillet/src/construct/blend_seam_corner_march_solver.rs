//! PK-152 stage 2 marching solver: the strip×strip intersection
//! curve, marched exactly per sample from the far anchor to the
//! support-plane exit (see `blend_seam_corner_march` for the repair
//! orchestration and authoring).

use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::types::BlendSeamCornerDeclineReason as Reason;

use super::super::terminal_arc_trim::StripCylinder;

/// Support plane extracted from a face.
pub(super) struct SupportPlane {
    pub(super) origin: Point3,
    pub(super) normal: Vec3,
    pub(super) u_axis: Vec3,
    pub(super) v_axis: Vec3,
}

pub(super) fn face_plane(brep: &BRep, face: FaceHandle) -> Option<SupportPlane> {
    match brep.faces.get(face)?.surface.as_ref().map(|s| &s.geometry) {
        Some(SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        }) => Some(SupportPlane {
            origin: *origin,
            normal: *normal,
            u_axis: *u_axis,
            v_axis: *v_axis,
        }),
        _ => None,
    }
}

/// March the strip×strip intersection from `start` (on both
/// cylinders) until it crosses the support plane; returns exact
/// per-sample points from `start` to the bisected exit point.
pub(super) fn march_seam(
    cyl_on: &StripCylinder,
    cyl_off: &StripCylinder,
    start: Point3,
    toward: Point3,
    plane: &SupportPlane,
    tolerance: &Tolerance,
) -> Result<Vec<Point3>, Reason> {
    let (p1, p2) = orthonormal_basis(cyl_on.axis);
    let r = cyl_on.radius;
    let a_on = cyl_on.axis;
    let a_off = cyl_off.axis;
    let c = a_on.dot(&a_off);
    let cross_len = a_on.cross(&a_off).length();
    if cross_len < 1e-6 {
        return Err(Reason::ParallelAxes);
    }
    let point_on = |u: f64, t: f64| -> Point3 {
        cyl_on.origin + p1 * (r * u.cos()) + p2 * (r * u.sin()) + a_on * t
    };
    // The strip_off membership constraint is quadratic in the axial
    // coordinate t at fixed chart angle u.
    let solve_roots = |u: f64| -> Option<(f64, f64)> {
        let base = cyl_on.origin + p1 * (r * u.cos()) + p2 * (r * u.sin());
        let d0 = base - cyl_off.origin;
        let a_ = 1.0 - c * c;
        let b_ = d0.dot(&a_on) - c * d0.dot(&a_off);
        let c0 = d0.dot(&d0) - d0.dot(&a_off).powi(2) - cyl_off.radius * cyl_off.radius;
        let disc = b_ * b_ - a_ * c0;
        if disc < 0.0 {
            return None;
        }
        let sq = disc.sqrt();
        Some(((-b_ + sq) / a_, (-b_ - sq) / a_))
    };
    let solve_t = |u: f64, t_prev: f64| -> Option<f64> {
        let (t1, t2) = solve_roots(u)?;
        Some(if (t1 - t_prev).abs() <= (t2 - t_prev).abs() {
            t1
        } else {
            t2
        })
    };
    // The start anchor is the SINGULAR point of the intersection
    // quartic (both strips tangent to the shared support there): the
    // first step's two roots are symmetric around the start station
    // and continuity cannot disambiguate them. Resolve the first
    // step toward the corner region instead.
    let solve_t_first = |u: f64| -> Option<f64> {
        let (t1, t2) = solve_roots(u)?;
        let base = |t: f64| cyl_on.origin + p1 * (r * u.cos()) + p2 * (r * u.sin()) + a_on * t;
        Some(
            if (base(t1) - toward).length() <= (base(t2) - toward).length() {
                t1
            } else {
                t2
            },
        )
    };
    let d_start = start - cyl_on.origin;
    let radial_start = d_start - a_on * d_start.dot(&a_on);
    let u0 = radial_start.dot(&p2).atan2(radial_start.dot(&p1));
    let t0 = d_start.dot(&a_on);
    let s_of = |p: Point3| -> f64 { (p - plane.origin).dot(&plane.normal) };
    let s_start = s_of(start);
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if s_start.abs() <= pos_tol {
        return Err(Reason::MarchUnresolved);
    }
    let side = s_start.signum();
    // Pick the march direction that closes on the support plane.
    let step = 0.02_f64;
    let mut dir = 0.0;
    for cand in [1.0_f64, -1.0] {
        if let Some(t) = solve_t(u0 + cand * step, t0) {
            let p = point_on(u0 + cand * step, t);
            if s_of(p) * side < s_start * side {
                dir = cand;
                break;
            }
        }
    }
    if dir == 0.0 {
        return Err(Reason::MarchUnresolved);
    }
    // Exit angle. A fillet SUPPORT plane is tangent to the strip
    // along the rail generator, so the seam meets it tangentially —
    // there is no sign change to bisect. The tangency angle is
    // analytic: the radial direction from the axis toward the plane.
    // A transversal plane (axis closer than r) keeps the
    // sign-crossing exit instead.
    let n_ax = plane.normal - a_on * plane.normal.dot(&a_on);
    let s_axis = (cyl_on.origin - plane.origin).dot(&plane.normal);
    let n_perp = n_ax.normalized();
    let ue = match n_perp {
        Some(n_perp) if (s_axis.abs() - r).abs() <= pos_tol.max(r * 1e-6) => {
            // Tangent support: exit at the rail generator's angle,
            // unwrapped in the march direction.
            let d_t = -n_perp * s_axis.signum();
            let u_rail = d_t.dot(&p2).atan2(d_t.dot(&p1));
            let mut du = (u_rail - u0) % std::f64::consts::TAU;
            if du < 0.0 {
                du += std::f64::consts::TAU;
            }
            if dir < 0.0 {
                du -= std::f64::consts::TAU;
            }
            if du.abs() < 1e-9 {
                return Err(Reason::MarchUnresolved);
            }
            u0 + du
        }
        _ => {
            // Transversal plane: step until the seam crosses it,
            // then bisect.
            let mut u_prev = u0;
            let mut t_prev = t0;
            let mut u_exit = None;
            for k in 1..=400 {
                let u = u0 + dir * step * k as f64;
                let Some(t) = solve_t(u, t_prev) else {
                    return Err(Reason::MarchUnresolved);
                };
                if s_of(point_on(u, t)) * side <= 0.0 {
                    let (mut lo, mut hi) = (u_prev, u);
                    let mut t_seed = t_prev;
                    for _ in 0..80 {
                        let mid = 0.5 * (lo + hi);
                        let Some(tm) = solve_t(mid, t_seed) else {
                            return Err(Reason::MarchUnresolved);
                        };
                        t_seed = tm;
                        if s_of(point_on(mid, tm)) * side <= 0.0 {
                            hi = mid;
                        } else {
                            lo = mid;
                        }
                    }
                    u_exit = Some(0.5 * (lo + hi));
                    break;
                }
                u_prev = u;
                t_prev = t;
            }
            u_exit.ok_or(Reason::MarchUnresolved)?
        }
    };
    // Uniform exact resample from start to exit; the seam must stay
    // on the body side of the plane the whole way.
    let n = 65;
    let mut pts = Vec::with_capacity(n);
    pts.push(start);
    let mut t_seed = t0;
    for k in 1..n {
        let u = u0 + (ue - u0) * k as f64 / (n - 1) as f64;
        let solved = if k == 1 {
            solve_t_first(u)
        } else {
            solve_t(u, t_seed)
        };
        let Some(t) = solved else {
            return Err(Reason::MarchUnresolved);
        };
        t_seed = t;
        let p = point_on(u, t);
        if s_of(p) * side < -pos_tol {
            return Err(Reason::MarchUnresolved);
        }
        pts.push(p);
    }
    // The exit point must land on the plane.
    let w = *pts.last().ok_or(Reason::MarchUnresolved)?;
    if s_of(w).abs() > pos_tol {
        return Err(Reason::MarchUnresolved);
    }
    Ok(pts)
}
