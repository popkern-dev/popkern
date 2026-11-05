//! Revolve end-condition resolution, extracted from `resolve.rs`.

use popkern_core::geometry::{CurveGeometry, Point3, Vec3};
use popkern_core::naming::EntityName;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_math::PortableMath;

use crate::error::ReplayError;
use crate::feature::{OperationMode, RevolveEnd};
use crate::parameter::ParameterTable;

use super::{resolve_face, resolve_param};

/// Context needed for [`RevolveEnd::UpToSurface`] resolution.
pub(crate) struct RevolveContext<'a> {
    /// Axis origin for the revolve.
    pub axis_origin: Point3,
    /// Axis direction (unit vector) for the revolve.
    pub axis_dir: Vec3,
    /// Centroid of the profile face being revolved.
    pub profile_centroid: Point3,
    /// Points of the profile face each end angle is checked from.
    pub profile_points: Vec<Point3>,
    /// BRep containing the profile.
    pub brep: &'a BRep,
    /// The body owning the target face, when it is another body (its
    /// boss merge host).
    pub host: Option<&'a BRep>,
    /// The revolve's mode: only a boss follows a face the section
    /// meets at several angles.
    pub mode: OperationMode,
    /// Geometric tolerance.
    pub tolerance: &'a Tolerance,
}

/// Resolve a [`RevolveEnd`] to a concrete angle in radians.
///
/// - [`RevolveEnd::Blind`]: looks up the angle parameter value.
/// - [`RevolveEnd::FullCircle`]: returns `2π` (τ).
/// - [`RevolveEnd::MidPlane`]: looks up the half_angle parameter value;
///   callers are responsible for the bidirectional revolve logic.
/// - [`RevolveEnd::UpToSurface`]: traces the profile centroid's circular
///   arc around the axis and finds where it intersects the target surface.
pub(crate) fn resolve_revolve_angle(
    end: &RevolveEnd,
    param_table: &ParameterTable,
    ctx: Option<&RevolveContext<'_>>,
) -> Result<f64, ReplayError> {
    match end {
        RevolveEnd::Blind { angle } => resolve_param(param_table, *angle),
        RevolveEnd::FullCircle => Ok(std::f64::consts::TAU),
        RevolveEnd::MidPlane { half_angle } => resolve_param(param_table, *half_angle),
        RevolveEnd::UpToSurface { face } => {
            let ctx = ctx.ok_or_else(|| ReplayError::FeatureFailed {
                label: "revolve".into(),
                error: "RevolveEnd::UpToSurface requires revolve context".into(),
            })?;
            resolve_revolve_up_to_surface(face, ctx)
        }
    }
}

/// The angle, about the axis and in the revolve's sense, from `p` to
/// where its circle first meets `surface` (`None` when it never does).
fn first_hit_angle(
    p: Point3,
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    range: ((f64, f64), (f64, f64)),
    ctx: &RevolveContext<'_>,
) -> Result<Option<f64>, ReplayError> {
    let axis = ctx.axis_dir;
    let along = (p - ctx.axis_origin).dot(&axis);
    let center = ctx.axis_origin + axis * along;
    let radial = p - center;
    let radius = radial.length();
    if radius < ctx.tolerance.linear() {
        return Ok(None);
    }
    let circle = CurveGeometry::Circle {
        center,
        normal: axis,
        radius,
    };
    let result = popkern_core::geometry::intersection::intersect_surface_curve(
        surface,
        range.0,
        range.1,
        &circle,
        (0.0, std::f64::consts::TAU),
        ctx.tolerance,
    )?;
    // Measured from `p` itself: the circle's parameter starts on an
    // arbitrary perpendicular of the axis (`circle_basis`), not at the
    // profile, so it is no angle of the revolve.
    let angle_to = |q: Point3| {
        let to = q - center;
        axis.dot(&radial.cross(&to))
            .patan2(radial.dot(&to))
            .rem_euclid(std::f64::consts::TAU)
    };
    let floor = ctx.tolerance.linear() / radius;
    Ok(result
        .intersections
        .iter()
        .filter_map(|i| match i {
            popkern_core::geometry::intersection::SurfaceCurveIntersection::Point {
                point, ..
            } => Some(angle_to(*point)),
            _ => None,
        })
        .filter(|a| *a > floor)
        .fold(None, |best: Option<f64>, a| {
            Some(best.map_or(a, |b| b.min(a)))
        }))
}

/// Resolve UpToSurface for a revolve: the angle at which the section
/// meets the target face. Every profile point turns on its own circle;
/// when they meet the face at one angle (a stop through the axis) that
/// is the end. Otherwise only a boss into the face's own body follows
/// it, turned past the farthest meeting with everything beyond the face
/// checked to lie in that body (as for an extrude, `extrude_conform`).
fn resolve_revolve_up_to_surface(
    face_name: &EntityName,
    ctx: &RevolveContext<'_>,
) -> Result<f64, ReplayError> {
    let target = ctx.host.unwrap_or(ctx.brep);
    let face_handle = resolve_face(face_name, target)?;
    let sof =
        target.faces[face_handle]
            .surface
            .as_ref()
            .ok_or_else(|| ReplayError::FeatureFailed {
                label: "revolve".into(),
                error: "target face has no surface geometry".into(),
            })?;
    let surface = &sof.geometry;
    let positions = target.face_outer_loop_positions(face_handle);
    let range = popkern_core::geometry::uv_projection::bounded_parameter_range(&positions, surface);
    let declined = crate::extrude_conform::declined;
    let at_centroid =
        first_hit_angle(ctx.profile_centroid, surface, range, ctx)?.ok_or_else(|| {
            ReplayError::FeatureFailed {
                label: "revolve".into(),
                error: "revolving profile does not intersect target surface".into(),
            }
        })?;
    // Each point's angle, and the radius its circle turns on.
    let mut hits = Vec::with_capacity(ctx.profile_points.len());
    for p in &ctx.profile_points {
        let r = (*p - ctx.axis_origin).cross(&ctx.axis_dir).length();
        if r < ctx.tolerance.linear() {
            continue;
        }
        match first_hit_angle(*p, surface, range, ctx)? {
            Some(a) => hits.push((*p, a, r)),
            None => return Err(declined("face_out_of_reach", at_centroid, at_centroid)),
        }
    }
    let near = hits.iter().map(|h| h.1).fold(at_centroid, f64::min);
    let far = hits.iter().map(|h| h.1).fold(at_centroid, f64::max);
    let reach = hits.iter().map(|h| h.2).fold(0.0, f64::max);
    // One angle to the tolerance at the outermost point.
    if (far - near) * reach <= ctx.tolerance.linear() {
        return Ok(at_centroid);
    }
    let Some(body) = ctx.host else {
        return Err(declined("curved_end_needs_boss", near, far));
    };
    if !matches!(ctx.mode, OperationMode::Boss) {
        return Err(declined("curved_end_needs_boss", near, far));
    }
    let end_at = (far + 0.05 * (far - near) + 10.0 * ctx.tolerance.linear() / reach.max(1e-300))
        .min(std::f64::consts::TAU);
    let boxes = popkern_bool::face_bounding_boxes(body);
    let inside = |q: Point3| {
        !matches!(
            popkern_bool::point_in_solid(&q, body, &boxes, ctx.tolerance),
            Ok(popkern_bool::PointClassification::Outside) | Err(_)
        )
    };
    let turned = |p: Point3, a: f64| {
        popkern_core::geometry::revolution::rotate_point(&p, &ctx.axis_origin, &ctx.axis_dir, a)
    };
    for (p, a, r) in &hits {
        let past = ctx.tolerance.linear() / r;
        if [(a + end_at) / 2.0, end_at]
            .iter()
            .any(|t| *t > a + past && !inside(turned(*p, *t)))
        {
            return Err(declined("thin_behind_face", near, far));
        }
    }
    Ok(end_at)
}
