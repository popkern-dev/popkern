//! Revolve end-condition resolution, extracted from `resolve.rs`.

use popkern_core::geometry::{CurveGeometry, Point3, Vec3};
use popkern_core::naming::EntityName;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::ReplayError;
use crate::feature::RevolveEnd;
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
    /// BRep containing the target surface.
    pub brep: &'a BRep,
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

/// Resolve UpToSurface for revolve by tracing the profile centroid's
/// circular arc and intersecting with the target surface.
fn resolve_revolve_up_to_surface(
    face_name: &EntityName,
    ctx: &RevolveContext<'_>,
) -> Result<f64, ReplayError> {
    let face_handle = resolve_face(face_name, ctx.brep)?;
    let sof = ctx.brep.faces[face_handle]
        .surface
        .as_ref()
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "revolve".into(),
            error: "target face has no surface geometry".into(),
        })?;
    let surface = &sof.geometry;

    // Build the circular arc of the profile centroid revolving around the axis.
    // Project centroid onto the axis to find the circle center and radius.
    let to_centroid = ctx.profile_centroid - ctx.axis_origin;
    let along_axis = to_centroid.dot(&ctx.axis_dir);
    let circle_center = ctx.axis_origin + ctx.axis_dir * along_axis;
    let radial = ctx.profile_centroid - circle_center;
    let radius = radial.length();

    if radius < ctx.tolerance.linear {
        return Err(ReplayError::FeatureFailed {
            label: "revolve".into(),
            error: "profile centroid is on the revolve axis".into(),
        });
    }

    // Build a Circle curve geometry for the centroid's path
    let circle_curve = CurveGeometry::Circle {
        center: circle_center,
        normal: ctx.axis_dir,
        radius,
    };

    // Get bounded UV range for the target surface
    let positions = ctx.brep.face_outer_loop_positions(face_handle);
    let (u_range, v_range) =
        popkern_core::geometry::uv_projection::bounded_parameter_range(&positions, surface);

    // Intersect the circle with the target surface (full revolution range)
    let result = popkern_core::geometry::intersection::intersect_surface_curve(
        surface,
        u_range,
        v_range,
        &circle_curve,
        (0.0, std::f64::consts::TAU),
        ctx.tolerance,
    )?;

    // Find the nearest forward intersection angle
    let mut min_t = f64::MAX;
    for isect in &result.intersections {
        if let popkern_core::geometry::intersection::SurfaceCurveIntersection::Point { t, .. } =
            isect
            && *t > ctx.tolerance.angular
            && *t < min_t
        {
            min_t = *t;
        }
    }

    if min_t == f64::MAX {
        return Err(ReplayError::FeatureFailed {
            label: "revolve".into(),
            error: "revolving profile does not intersect target surface".into(),
        });
    }

    Ok(min_t)
}
