//! UV parameter-space projection utilities.
//!
//! Projects 3D points onto surface UV spaces and computes bounded parameter
//! ranges for faces with infinite parameter domains.

use std::f64::consts::TAU;

use crate::geometry::Point3;
use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Project a 3D point onto a surface's UV parameter space.
///
/// Returns `Some((u, v))` for analytic surface types (closed-form inverse)
/// and for NURBS surfaces (iterative Gauss-Newton projection).
/// Returns `None` only when no starting point converges close enough.
pub fn project_point_to_uv(surface: &SurfaceGeometry, point: &Point3) -> Option<(f64, f64)> {
    match surface {
        SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        } => {
            let d = *point - *origin;
            let uu = u_axis.dot(u_axis);
            let uv = u_axis.dot(v_axis);
            let vv = v_axis.dot(v_axis);
            let du = d.dot(u_axis);
            let dv = d.dot(v_axis);
            let det = uu * vv - uv * uv;
            if det.abs() <= 1e-30 {
                return None;
            }
            let u = (du * vv - dv * uv) / det;
            let v = (dv * uu - du * uv) / det;
            Some((u, v))
        }
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius: _,
        } => {
            let (perp1, perp2) = orthonormal_basis(*axis);
            let d = *point - *origin;
            let mut u = f64::atan2(d.dot(&perp2), d.dot(&perp1));
            if u < 0.0 {
                u += TAU;
            }
            let v = d.dot(axis);
            Some((u, v))
        }
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle: _,
        } => {
            let (perp1, perp2) = orthonormal_basis(*axis);
            let d = *point - *apex;
            let mut u = f64::atan2(d.dot(&perp2), d.dot(&perp1));
            if u < 0.0 {
                u += TAU;
            }
            let v = d.dot(axis);
            Some((u, v))
        }
        SurfaceGeometry::Sphere { center, radius } => {
            let d = *point - *center;
            let r = d.length();
            // Avoid division by zero at the centre.
            if r < 1e-15 {
                return Some((0.0, 0.0));
            }
            let dz = d.z / r;
            // Clamp to avoid NaN from acos outside [-1, 1].
            let v = dz.clamp(-1.0, 1.0).acos();
            let mut u = f64::atan2(d.y, d.x);
            if u < 0.0 {
                u += TAU;
            }
            let _ = radius; // used only for the forward map
            Some((u, v))
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius: _,
        } => {
            let (perp1, perp2) = orthonormal_basis(*axis);
            let d = *point - *center;
            // Project onto the equatorial plane to find the major angle u.
            let radial_1 = d.dot(&perp1);
            let radial_2 = d.dot(&perp2);
            let mut u = f64::atan2(radial_2, radial_1);
            if u < 0.0 {
                u += TAU;
            }
            // Find the minor angle v: direction from the major-circle point
            // to the actual point, decomposed into (radial, axial).
            let h = d.dot(axis);
            let radial_dist = (radial_1 * radial_1 + radial_2 * radial_2).sqrt();
            let mut v = f64::atan2(h, radial_dist - major_radius);
            if v < 0.0 {
                v += TAU;
            }
            Some((u, v))
        }
        // NURBS: the clamped 5-start multistart, UNCHANGED for the
        // cold construction paths — trim/clip/pcurve authoring sites
        // depend on its basin behavior (swapping in the knot-seeded
        // multi-sheet projection, even as best-of-both, moved a
        // boss-plate fillet-strip projection onto another fold and
        // the served body carried a real 5.2e-3 curve-surface lie).
        // Multi-sheet-correct DISTANCE measurement lives in
        // [`project_point_to_uv_multisheet`] for validation
        // consumers (PK-449).
        SurfaceGeometry::Nurbs(nurbs) => {
            uv_projection_clamped::nurbs_gauss_newton_project(surface, nurbs, point).or_else(|| {
                // Rescue on HARD failure only (all five starts ended
                // >1 unit away — every basin walked onto another
                // sheet): the knot-span-seeded multi-start resolves
                // coiled multi-sheet charts (PK-446/PK-449). Scoped
                // to the None path so every succeeding projection
                // keeps its tuned basin byte-identically (measured,
                // PK-463: a groove-flank trim endpoint ON the chart
                // border projected >1 off through all five starts,
                // the pcurve recompute silently skipped the
                // half-edge, and the served body meshed torn).
                let seeds = nurbs.projection_seeds()?;
                crate::geometry::surface::closest_point_from_seeds(
                    nurbs,
                    point,
                    &seeds,
                    (false, false),
                    6,
                )
            })
        }
        // Gregory patches are non-analytic-inverse; fall through
        // to the generic Gauss-Newton projection. The default
        // [`SurfaceEval::closest_point`] dispatches via
        // `closest_point_gauss_newton` for any surface variant
        // without a closed-form inverse.
        SurfaceGeometry::Gregory(_) | SurfaceGeometry::PolyhedralSpline(_) => {
            surface.closest_point(point)
        }
    }
}

#[path = "uv_projection_clamped.rs"]
mod uv_projection_clamped;

/// Multi-sheet-correct closest-preimage projection, for VALIDATION
/// consumers that measure true distance-to-surface (check 9, gap
/// diagnostics). On a multi-span NURBS chart this takes the best of
/// the knot-span-seeded projection (resolves multi-sheet basins —
/// the 5-start grid walked every start into the wrong turn of a
/// 6.5-turn coiled wall and reported a 0.4168 phantom deviation for
/// an on-seam point, PK-446/PK-449) and the clamped 5-start grid
/// (broad exploration where coarse seed ranking cuts the right span
/// from the refined few). Construction paths keep
/// [`project_point_to_uv`]: they depend on its basin behavior, and a
/// closer preimage on ANOTHER FOLD is the wrong answer there.
pub fn project_point_to_uv_multisheet(
    surface: &SurfaceGeometry,
    point: &Point3,
) -> Option<(f64, f64)> {
    if let SurfaceGeometry::Nurbs(nurbs) = surface
        && let Some(seeds) = nurbs.projection_seeds()
    {
        let seeded = crate::geometry::surface::closest_point_from_seeds(
            nurbs,
            point,
            &seeds,
            (false, false),
            6,
        );
        let grid = uv_projection_clamped::nurbs_gauss_newton_project(surface, nurbs, point);
        return match (seeded, grid) {
            (Some(a), Some(b)) => {
                let da = point.distance_to(&surface.point_at(a.0, a.1));
                let db = point.distance_to(&surface.point_at(b.0, b.1));
                Some(if da <= db { a } else { b })
            }
            (a, b) => a.or(b),
        };
    }
    project_point_to_uv(surface, point)
}

/// Project a 3D point to UV using the previous UV as starting guess.
///
/// For analytic surfaces, delegates to [`project_point_to_uv`] (exact,
/// no iteration needed). For NURBS, uses single-start Gauss-Newton from
/// `(u_hint, v_hint)` instead of the multi-start search. This prevents
/// branch jumping when projecting consecutive points along a curve.
pub fn project_point_to_uv_near(
    surface: &SurfaceGeometry,
    point: &Point3,
    u_hint: f64,
    v_hint: f64,
) -> Option<(f64, f64)> {
    match surface {
        // Analytic surfaces have unique inverse — hint is irrelevant.
        SurfaceGeometry::Nurbs(nurbs) => {
            uv_projection_clamped::nurbs_gauss_newton_single(surface, nurbs, point, u_hint, v_hint)
        }
        // Gregory patches have no analytic inverse either, and the
        // global `closest_point` multistart can escape the [0, 1]²
        // domain near a degenerate pole corner (measured on the
        // hexagon6_variable setback fan: v walking to 1.45 with
        // 2.5e-2 off-surface images while an in-domain preimage at
        // 1e-5 exists). Hint-seeded, domain-clamped Gauss-Newton
        // honours the continuation contract.
        // The PnS patch shares the Gregory situation exactly: no
        // analytic inverse, unit-square domain, and a degenerate pole
        // corner the global multistart can escape from.
        SurfaceGeometry::Gregory(_) | SurfaceGeometry::PolyhedralSpline(_) => {
            uv_projection_clamped::gauss_newton_single_clamped(
                surface,
                point,
                u_hint,
                v_hint,
                (0.0, 1.0),
                (0.0, 1.0),
            )
        }
        _ => project_point_to_uv(surface, point),
    }
}

/// Compute bounded parameter ranges for a face's surface.
///
/// If the surface's native parameter range is already finite in both u and v,
/// returns it unchanged. Otherwise, projects the given vertex positions into UV
/// space, computes the axis-aligned bounding box, and adds 50 % padding
/// (minimum 1.0 in each direction) for intersection-solver headroom.
pub fn bounded_parameter_range(
    vertices: &[Point3],
    surface: &SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();

    // If both ranges are finite AND the surface is non-periodic,
    // the native range already describes the face extent. For
    // periodic surfaces (Torus, Cylinder, Sphere) the native
    // range is the full period — the face may only occupy a
    // strip, so we fall through to vertex-based projection.
    let is_periodic = matches!(
        surface,
        SurfaceGeometry::Torus { .. }
            | SurfaceGeometry::Cylinder { .. }
            | SurfaceGeometry::Sphere { .. }
    );
    if !is_periodic
        && u_min.is_finite()
        && u_max.is_finite()
        && v_min.is_finite()
        && v_max.is_finite()
    {
        return ((u_min, u_max), (v_min, v_max));
    }

    if vertices.is_empty() {
        return ((u_min, u_max), (v_min, v_max));
    }

    // Project each vertex to UV and accumulate the bounding box.
    let mut uv_lo = (f64::INFINITY, f64::INFINITY);
    let mut uv_hi = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut any_projected = false;

    for pos in vertices {
        if let Some((u, v)) = project_point_to_uv(surface, pos) {
            uv_lo.0 = uv_lo.0.min(u);
            uv_lo.1 = uv_lo.1.min(v);
            uv_hi.0 = uv_hi.0.max(u);
            uv_hi.1 = uv_hi.1.max(v);
            any_projected = true;
        }
    }

    if !any_projected {
        return ((u_min, u_max), (v_min, v_max));
    }

    // 50 % padding (minimum 1.0) so downstream intersection solvers have
    // enough room; precise 3-D trimming happens later.
    let u_pad = ((uv_hi.0 - uv_lo.0) * 0.5).max(1.0);
    let v_pad = ((uv_hi.1 - uv_lo.1) * 0.5).max(1.0);

    // Only override dimensions that were infinite.
    let u_range = if u_min.is_finite() && u_max.is_finite() {
        (u_min, u_max)
    } else {
        (uv_lo.0 - u_pad, uv_hi.0 + u_pad)
    };
    let v_range = if v_min.is_finite() && v_max.is_finite() {
        (v_min, v_max)
    } else {
        (uv_lo.1 - v_pad, uv_hi.1 + v_pad)
    };

    (u_range, v_range)
}

/// Compute tight bounded parameter ranges for a face's surface.
///
/// Same projection algorithm as [`bounded_parameter_range`] but with
/// 5 % padding (minimum 0.05) instead of 50 %. Suitable for
/// tessellation and AABB sampling where a small margin is needed.
pub fn tight_bounded_parameter_range(
    vertices: &[Point3],
    surface: &SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    vertex_bounded_parameter_range(vertices, surface, 0.05, 0.05)
}

/// Compute exact bounded parameter ranges from face vertices with no
/// padding.
///
/// Returns the exact UV bounding box of the given vertex projections.
/// Suitable for self-intersection validation where any padding would
/// cause false positives at shared edges.
pub fn exact_bounded_parameter_range(
    vertices: &[Point3],
    surface: &SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    vertex_bounded_parameter_range(vertices, surface, 0.0, 0.0)
}

/// Internal: compute bounded parameter ranges with configurable padding.
///
/// For analytic surfaces (Plane, Cylinder, Cone, Sphere, Torus), projects
/// face vertices into UV space and uses the tighter of (native range,
/// projected range + padding). This ensures faces using a subset of a
/// torus or sphere are tessellated only over the portion they occupy.
///
/// For NURBS surfaces, returns the native knot-vector range as-is. NURBS
/// surfaces from sweep/fillet pipelines have knot ranges that already
/// match the face extent, and vertex projection on NURBS with degenerate
/// topology (seam vertices) gives misleading UV ranges.
fn vertex_bounded_parameter_range(
    vertices: &[Point3],
    surface: &SurfaceGeometry,
    padding_fraction: f64,
    padding_minimum: f64,
) -> ((f64, f64), (f64, f64)) {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();

    // NURBS: trust the knot-vector range. The construction pipeline
    // (sweep, fillet, boolean) builds NURBS surfaces whose knot ranges
    // match the face extent. Vertex projection on NURBS with degenerate
    // seam topology produces misleadingly tight ranges.
    if matches!(surface, SurfaceGeometry::Nurbs(_)) {
        return ((u_min, u_max), (v_min, v_max));
    }

    if vertices.is_empty() {
        return ((u_min, u_max), (v_min, v_max));
    }

    // Project each vertex to UV and accumulate the bounding box.
    let mut uv_lo = (f64::INFINITY, f64::INFINITY);
    let mut uv_hi = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut any_projected = false;

    for pos in vertices {
        if let Some((u, v)) = project_point_to_uv(surface, pos) {
            uv_lo.0 = uv_lo.0.min(u);
            uv_lo.1 = uv_lo.1.min(v);
            uv_hi.0 = uv_hi.0.max(u);
            uv_hi.1 = uv_hi.1.max(v);
            any_projected = true;
        }
    }

    if !any_projected {
        return ((u_min, u_max), (v_min, v_max));
    }

    let u_pad = ((uv_hi.0 - uv_lo.0) * padding_fraction).max(padding_minimum);
    let v_pad = ((uv_hi.1 - uv_lo.1) * padding_fraction).max(padding_minimum);

    let projected_u = (uv_lo.0 - u_pad, uv_hi.0 + u_pad);
    let projected_v = (uv_lo.1 - v_pad, uv_hi.1 + v_pad);

    // For each dimension: use the tighter of (native range, projected range).
    // If native range is infinite, use projected. If native is finite but
    // wider than projected, use projected (clamped to native bounds).
    let u_range = if u_min.is_finite() && u_max.is_finite() {
        (projected_u.0.max(u_min), projected_u.1.min(u_max))
    } else {
        projected_u
    };
    let v_range = if v_min.is_finite() && v_max.is_finite() {
        (projected_v.0.max(v_min), projected_v.1.min(v_max))
    } else {
        projected_v
    };

    (u_range, v_range)
}

#[cfg(test)]
#[path = "uv_projection_tests.rs"]
mod tests;

pub use super::uv_projection_face_box::{face_bounded_parameter_range, face_exact_parameter_range};

/// Analytic, branch-aware point inversion for surfaces where the
/// generic Newton projection wanders across period branches.
/// `hint` (a previous sample's uv) selects the branch; without it,
/// the canonical branch is returned. Falls back to the generic
/// projection for surface kinds without an analytic inverse.
pub fn project_point_to_uv_analytic(
    surface: &SurfaceGeometry,
    p: &Point3,
    hint: Option<(f64, f64)>,
) -> Option<(f64, f64)> {
    let analytic = project_point_to_uv_analytic_unverified(surface, p, hint);
    // Round-trip contract: the analytic inversion assumes the
    // surface's chart uses `orthonormal_basis(axis)` with angle-u /
    // axial-v. Producers are not uniform (fillet strip cylinders
    // measured 5e-3 off through this inversion on the pentagon
    // setback gate), so verify by evaluation and fall back to the
    // generic projection on mismatch.
    if let Some((u, v)) = analytic {
        use crate::geometry::surface::SurfaceEval;
        if surface.point_at(u, v).distance_to(p) < 1.0e-9 {
            return Some((u, v));
        }
    }
    match hint {
        Some((hu, hv)) => project_point_to_uv_near(surface, p, hu, hv),
        None => project_point_to_uv(surface, p),
    }
}

fn project_point_to_uv_analytic_unverified(
    surface: &SurfaceGeometry,
    p: &Point3,
    hint: Option<(f64, f64)>,
) -> Option<(f64, f64)> {
    match surface {
        SurfaceGeometry::Cylinder { origin, axis, .. } => {
            let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(*axis);
            let d = *p - *origin;
            let mut u = d.dot(&perp2).atan2(d.dot(&perp1));
            if let Some((hu, _)) = hint {
                let tau = std::f64::consts::TAU;
                while u - hu > 0.5 * tau {
                    u -= tau;
                }
                while u - hu < -0.5 * tau {
                    u += tau;
                }
            }
            Some((u, d.dot(axis)))
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            ..
        } => {
            let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(*axis);
            let d = *p - *center;
            let mut u = d.dot(&perp2).atan2(d.dot(&perp1));
            let radial = d - *axis * d.dot(axis);
            let r_len = radial.length();
            let mut v = (d.dot(axis)).atan2(r_len - major_radius);
            if let Some((hu, hv)) = hint {
                let tau = std::f64::consts::TAU;
                while u - hu > 0.5 * tau {
                    u -= tau;
                }
                while u - hu < -0.5 * tau {
                    u += tau;
                }
                while v - hv > 0.5 * tau {
                    v -= tau;
                }
                while v - hv < -0.5 * tau {
                    v += tau;
                }
            }
            Some((u, v))
        }
        _ => match hint {
            Some((hu, hv)) => project_point_to_uv_near(surface, p, hu, hv),
            None => project_point_to_uv(surface, p),
        },
    }
}

#[cfg(test)]
#[path = "uv_projection_extended_tests.rs"]
mod extended_tests;
