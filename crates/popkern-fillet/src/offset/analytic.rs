//! Analytic (exact) surface offset for primitive surface types.
//!
//! Direction contract: `inward = true` means "offset along the stored
//! geometric normal's NEGATIVE direction" (-n), uniformly across all
//! arms. Radius-form surfaces (cylinder/sphere/torus) translate that
//! into a radius delta through the measured radial sign of their stored
//! normal — a parametrization whose normal points toward the axis must
//! GROW its radius under an along--n offset.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::error::OffsetError;

/// Fraction of the exact existence bound reported as `max_safe` in
/// radius-through-zero and spindle declines, matching the collapse
/// gate's convention that `max_safe` must itself be serveable.
const MAX_SAFE_FRACTION: f64 = 0.99;

/// Sign of the stored normal along the radius-growing direction: `+1.0`
/// when the surface normal points away from its axis/center (radius
/// grows along `+n`), `-1.0` when it points toward it. `None` for
/// surface classes without a radius form.
pub fn radial_normal_sign(surface: &SurfaceGeometry) -> Option<f64> {
    let radial = match surface {
        SurfaceGeometry::Cylinder { origin, axis, .. } => {
            let p = surface.point_at(0.0, 0.0);
            let axis_n = axis.normalized()?;
            let rel = p - *origin;
            rel - axis_n * rel.dot(&axis_n)
        }
        SurfaceGeometry::Sphere { center, .. } => {
            let p = surface.point_at(0.0, 0.5);
            p - *center
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            ..
        } => {
            let p = surface.point_at(0.0, 0.0);
            let axis_n = axis.normalized()?;
            let rel = p - *center;
            let ring_dir = (rel - axis_n * rel.dot(&axis_n)).normalized()?;
            let ring_center = *center + ring_dir * *major_radius;
            p - ring_center
        }
        _ => return None,
    };
    let n = match surface {
        SurfaceGeometry::Sphere { .. } => surface.normal_at(0.0, 0.5),
        _ => surface.normal_at(0.0, 0.0),
    };
    let dot = n.dot(&radial);
    if dot == 0.0 {
        return None;
    }
    Some(if dot > 0.0 { 1.0 } else { -1.0 })
}

/// Outcome of attempting a closed-form analytic offset.
#[derive(Debug, Clone)]
pub enum AnalyticOffset {
    /// The surface class has a closed-form offset and it exists at the
    /// requested distance.
    Exact(SurfaceGeometry),
    /// No closed-form arm exists for this surface class (NURBS,
    /// Gregory, polyhedral spline) — approximation is the honest route.
    Unsupported,
}

/// Compute an exact offset of an analytic surface, distinguishing
/// "no closed-form arm for this class" from "the offset does not exist".
///
/// For a plane, the origin is translated by `±distance * normal`. For a
/// cone, the apex is translated along the (normalized) axis by
/// `∓ distance / sin(half_angle)` — the same half-angle cone, shifted so
/// every point sits `distance` along the surface normal (offset along
/// `+n̂` moves the apex by `−d/sinα · axis`). The offset cone always
/// exists as a surface; a through-apex request (the offset consumes the
/// original apex region) is a *trim* concern handled downstream by the
/// recompute/collapse machinery, never a reason to approximate. For
/// cylinder, sphere, and torus the radius is adjusted by the offset
/// distance folded with the stored normal's measured radial sign, so
/// `inward` uniformly means "offset along −n" (see the module doc).
///
/// # Errors
///
/// `Err` is a genuine geometric decline: the offset surface does not
/// exist at the requested distance (a radius driven through zero, or a
/// torus minor radius grown to a self-intersecting spindle). Callers
/// must NOT fall back to approximation on `Err` — an approximation of a
/// nonexistent offset is inverted geometry. Surface classes without a
/// closed form return [`AnalyticOffset::Unsupported`] instead.
pub fn offset_analytic_classified(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
) -> Result<AnalyticOffset, OffsetError> {
    let sign = if inward { -1.0 } else { 1.0 };

    match surface {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        } => {
            let offset_vec: Vec3 = *normal * (sign * distance);
            Ok(AnalyticOffset::Exact(SurfaceGeometry::Plane {
                origin: *origin + offset_vec,
                normal: *normal,
                u_axis: *u_axis,
                v_axis: *v_axis,
            }))
        }

        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            let sin_a = half_angle.sin();
            let cos_a = half_angle.cos();
            let axis_n = axis
                .normalized()
                .filter(|_| sin_a > 0.0 && cos_a > 0.0)
                .ok_or(OffsetError::DegenerateGeometry { class: "cone" })?;
            // For v > 0 the stored normal (∂u × ∂v) points along
            // cosα·e_r − sinα·â; offsetting every surface point by
            // `sign·distance` along that unit normal reproduces the
            // same cone with the apex translated −(sign·distance)/sinα
            // along the unit axis.
            let new_apex = *apex - axis_n * (sign * distance / sin_a);
            Ok(AnalyticOffset::Exact(SurfaceGeometry::Cone {
                apex: new_apex,
                axis: *axis,
                half_angle: *half_angle,
            }))
        }

        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let radial = radial_normal_sign(surface)
                .ok_or(OffsetError::DegenerateGeometry { class: "cylinder" })?;
            let new_radius = radius + sign * radial * distance;
            if new_radius <= 0.0 {
                return Err(OffsetError::RadiusThroughZero {
                    class: "cylinder",
                    new_radius,
                    max_safe: radius * MAX_SAFE_FRACTION,
                });
            }
            Ok(AnalyticOffset::Exact(SurfaceGeometry::Cylinder {
                origin: *origin,
                axis: *axis,
                radius: new_radius,
            }))
        }

        SurfaceGeometry::Sphere { center, radius } => {
            let radial = radial_normal_sign(surface)
                .ok_or(OffsetError::DegenerateGeometry { class: "sphere" })?;
            let new_radius = radius + sign * radial * distance;
            if new_radius <= 0.0 {
                return Err(OffsetError::RadiusThroughZero {
                    class: "sphere",
                    new_radius,
                    max_safe: radius * MAX_SAFE_FRACTION,
                });
            }
            Ok(AnalyticOffset::Exact(SurfaceGeometry::Sphere {
                center: *center,
                radius: new_radius,
            }))
        }

        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let radial = radial_normal_sign(surface)
                .ok_or(OffsetError::DegenerateGeometry { class: "torus" })?;
            let new_minor = minor_radius + sign * radial * distance;
            if new_minor <= 0.0 {
                return Err(OffsetError::RadiusThroughZero {
                    class: "torus",
                    new_radius: new_minor,
                    max_safe: minor_radius * MAX_SAFE_FRACTION,
                });
            }
            if new_minor >= *major_radius {
                // Spindle torus: the tube swallows the axis and the
                // surface self-intersects. Declined until global
                // trimming (PK-327) can serve the trimmed result.
                return Err(OffsetError::SpindleTorus {
                    new_minor,
                    major: *major_radius,
                    max_safe: (major_radius - minor_radius) * MAX_SAFE_FRACTION,
                });
            }
            Ok(AnalyticOffset::Exact(SurfaceGeometry::Torus {
                center: *center,
                axis: *axis,
                major_radius: *major_radius,
                minor_radius: new_minor,
            }))
        }

        SurfaceGeometry::Nurbs(_)
        | SurfaceGeometry::Gregory(_)
        | SurfaceGeometry::PolyhedralSpline(_) => Ok(AnalyticOffset::Unsupported),
    }
}

/// Compute an exact offset of an analytic surface.
///
/// Compatibility wrapper over [`offset_analytic_classified`]: surface
/// classes without a closed-form arm are reported as `Err` alongside
/// genuine geometric declines. Dispatch code that must distinguish the
/// two (to avoid approximating a nonexistent offset) uses the
/// classified form.
pub fn offset_analytic(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
) -> Result<SurfaceGeometry, OffsetError> {
    match offset_analytic_classified(surface, distance, inward)? {
        AnalyticOffset::Exact(s) => Ok(s),
        AnalyticOffset::Unsupported => {
            let class = match surface {
                SurfaceGeometry::Nurbs(_) => "NURBS",
                SurfaceGeometry::Gregory(_) => "Gregory",
                SurfaceGeometry::PolyhedralSpline(_) => "polyhedral spline",
                _ => "analytic",
            };
            Err(OffsetError::Unsupported { class })
        }
    }
}

/// The parameter window on the *offset* surface geometrically
/// corresponding to `(u_range, v_range)` on the original.
///
/// Every exact arm except the cone preserves parametrization
/// (`S'(u,v) = S(u,v) + sign·d·n̂(u,v)`), so the window is unchanged.
/// The exact cone offset translates the apex, which shifts the axial
/// parameter: a point at `(u, v)` maps to `(u, v + sign·d·cos²α/sinα)`
/// on the offset cone. Callers that intersect offset surfaces inside
/// windows derived from the original faces (the fillet center-curve
/// construction) must route their windows through this map or steep
/// cones lose the intersection region.
pub fn corresponding_offset_window(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> ((f64, f64), (f64, f64)) {
    if let SurfaceGeometry::Cone { half_angle, .. } = surface {
        let sin_a = half_angle.sin();
        let cos_a = half_angle.cos();
        if sin_a > 0.0 {
            let sign = if inward { -1.0 } else { 1.0 };
            let shift = sign * distance * cos_a * cos_a / sin_a;
            return (u_range, (v_range.0 + shift, v_range.1 + shift));
        }
    }
    (u_range, v_range)
}

#[cfg(test)]
#[path = "analytic_tests.rs"]
mod tests;
