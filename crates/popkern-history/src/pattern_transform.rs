//! Compute per-instance spatial transforms for pattern features.
//!
//! Each pattern type (linear, circular, curve-driven, mirror) produces a
//! [`Vec`] of [`Transform3`] values — one per instance.  Instance 0 is always
//! the identity (the seed position).
//!
//! The public-facing entry point is [`compute_pattern_transforms`], called by
//! [`crate::pattern_execute::execute_pattern`].

use crate::error::ReplayError;
use crate::feature::{AxisDesc, Feature};
use crate::parameter::ParameterTable;
use crate::pattern::MirrorPlaneDesc;
use crate::resolve;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;
use popkern_sweep::arc_length::ArcLengthTable;

/// Look up a parameter value by ID.
fn resolve_param(
    id: crate::parameter::ParameterId,
    table: &ParameterTable,
) -> Result<f64, ReplayError> {
    table
        .get(id)
        .map(|p| p.value)
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: String::new(),
            error: format!("parameter {id:?} not found"),
        })
}

/// Compute per-instance spatial transforms for a pattern feature.
///
/// Returns a tuple of ([`Vec`] of [`Transform3`], [`Vec`] of warning strings).
/// The transforms list has one entry per instance; instance 0 is always
/// [`Transform3::identity`] (the seed itself).  Warnings are non-fatal
/// messages (e.g. angle clamping) that callers should surface to the user.
///
/// # Errors
///
/// Returns [`ReplayError`] if entity resolution, axis computation, or
/// parameter lookup fails, or if `count == 0`.
pub(crate) fn compute_pattern_transforms(
    feature: &Feature,
    brep: &BRep,
    param_table: &ParameterTable,
) -> Result<(Vec<Transform3>, Vec<String>), ReplayError> {
    match feature {
        Feature::LinearPattern {
            primary, secondary, ..
        } => {
            let transforms =
                compute_linear_transforms(primary, secondary.as_ref(), brep, param_table)?;
            Ok((transforms, Vec::new()))
        }
        Feature::CircularPattern {
            axis,
            count,
            total_angle,
            ..
        } => compute_circular_transforms(axis, *count, *total_angle, brep, param_table),
        Feature::CurvePattern {
            curve,
            count,
            equal_spacing,
            align_to_curve,
            ..
        } => {
            let transforms = compute_curve_transforms(
                curve,
                *count,
                *equal_spacing,
                *align_to_curve,
                brep,
                param_table,
            )?;
            Ok((transforms, Vec::new()))
        }
        Feature::Mirror { plane, .. } => {
            let transforms = compute_mirror_transforms(plane, brep)?;
            Ok((transforms, Vec::new()))
        }
        _ => Err(ReplayError::FeatureFailed {
            label: String::new(),
            error: "compute_pattern_transforms called on non-pattern feature".into(),
        }),
    }
}

/// Compute transforms for a linear (rectangular grid) pattern.
fn compute_linear_transforms(
    primary: &crate::pattern::PatternDirection,
    secondary: Option<&crate::pattern::PatternDirection>,
    brep: &BRep,
    param_table: &ParameterTable,
) -> Result<Vec<Transform3>, ReplayError> {
    let dir1 = resolve::resolve_direction(&primary.direction, brep)?;
    let count1 = resolve_param(primary.count, param_table)?.round() as usize;
    let dist1 = resolve_param(primary.total_distance, param_table)?;

    if count1 == 0 {
        return Err(ReplayError::FeatureFailed {
            label: String::new(),
            error: "linear pattern primary count must be > 0".into(),
        });
    }

    let pitch1 = if count1 > 1 {
        dist1 / (count1 - 1) as f64
    } else {
        0.0
    };

    let (dir2, count2, pitch2) = if let Some(sec) = secondary {
        let d = resolve::resolve_direction(&sec.direction, brep)?;
        let c = resolve_param(sec.count, param_table)?.round() as usize;
        if c == 0 {
            return Err(ReplayError::FeatureFailed {
                label: String::new(),
                error: "linear pattern secondary count must be > 0".into(),
            });
        }
        let dist = resolve_param(sec.total_distance, param_table)?;
        let p = if c > 1 { dist / (c - 1) as f64 } else { 0.0 };
        (d, c, p)
    } else {
        (Vec3::zero(), 1, 0.0)
    };

    let total = count1 * count2;
    let mut transforms = Vec::with_capacity(total);
    // Row-major ordering: instance index = i * count2 + j.
    // Instance 0 is (i=0, j=0), instance 1 is (i=0, j=1), etc.
    for i in 0..count1 {
        for j in 0..count2 {
            let offset = dir1 * (i as f64 * pitch1) + dir2 * (j as f64 * pitch2);
            transforms.push(Transform3::translation(offset.x, offset.y, offset.z));
        }
    }
    Ok(transforms)
}

/// Compute transforms for a circular pattern.
///
/// Returns `(transforms, warnings)`.  If the requested total angle exceeds
/// 360 degrees, it is clamped and a warning is emitted.
fn compute_circular_transforms(
    axis: &AxisDesc,
    count_id: crate::parameter::ParameterId,
    total_angle_id: crate::parameter::ParameterId,
    brep: &BRep,
    param_table: &ParameterTable,
) -> Result<(Vec<Transform3>, Vec<String>), ReplayError> {
    let (origin, dir) = resolve::resolve_axis(axis, brep)?;
    let count = resolve_param(count_id, param_table)?.round() as usize;
    if count == 0 {
        return Err(ReplayError::FeatureFailed {
            label: String::new(),
            error: "circular pattern count must be > 0".into(),
        });
    }

    let raw_angle = resolve_param(total_angle_id, param_table)?;
    let mut warnings = Vec::new();
    let total_angle = if raw_angle.abs() > std::f64::consts::TAU {
        warnings.push(format!(
            "circular pattern total angle {:.1}\u{00b0} exceeds 360\u{00b0}; clamped to 360\u{00b0}",
            raw_angle.to_degrees(),
        ));
        raw_angle.signum() * std::f64::consts::TAU
    } else {
        raw_angle
    };

    let angle_step = if count > 1 {
        total_angle / (count - 1) as f64
    } else {
        0.0
    };

    let mut transforms = Vec::with_capacity(count);
    for i in 0..count {
        let angle = i as f64 * angle_step;
        let t = Transform3::rotation_around_axis(origin, dir, angle).map_err(|e| {
            ReplayError::FeatureFailed {
                label: String::new(),
                error: e.to_string(),
            }
        })?;
        transforms.push(t);
    }
    Ok((transforms, warnings))
}

/// Compute transforms for a curve-driven pattern.
///
/// Resolves the named edge, converts its geometry to a [`NurbsCurve`], and
/// distributes `count` instances along the curve.  When `equal_spacing` is
/// `true` the instances are placed at equal arc-length intervals; otherwise
/// they are placed at equal parameter intervals.
///
/// When `align_to_curve` is `true`, each instance is additionally rotated so
/// that its local frame tracks the curve tangent.
fn compute_curve_transforms(
    curve: &popkern_core::naming::EntityName,
    count_id: crate::parameter::ParameterId,
    equal_spacing: bool,
    align_to_curve: bool,
    brep: &BRep,
    param_table: &ParameterTable,
) -> Result<Vec<Transform3>, ReplayError> {
    let count = resolve_param(count_id, param_table)?.round() as usize;
    if count == 0 {
        return Err(ReplayError::FeatureFailed {
            label: String::new(),
            error: "curve pattern count must be > 0".into(),
        });
    }

    let edge_handle = resolve::resolve_edge(curve, brep)?;

    let curve_on_edge =
        brep.edges[edge_handle]
            .curve
            .as_ref()
            .ok_or_else(|| ReplayError::FeatureFailed {
                label: String::new(),
                error: "curve pattern edge has no curve geometry".into(),
            })?;

    // Convert to NurbsCurve for arc-length parameterization.
    let nurbs = curve_on_edge
        .geometry
        .to_nurbs_in_range((curve_on_edge.t_start, curve_on_edge.t_end));

    let (t0, t1) = nurbs.parameter_range();
    let start_pos = nurbs.point_at(t0);
    let start_tangent = nurbs.tangent_at(t0);

    // Build arc-length table for equal-spacing mode.
    let table = if equal_spacing {
        Some(ArcLengthTable::build(&nurbs, 100))
    } else {
        None
    };

    let mut transforms = Vec::with_capacity(count);
    for i in 0..count {
        let frac = if count > 1 {
            i as f64 / (count - 1) as f64
        } else {
            0.0
        };

        let t_param = if let Some(ref tbl) = table {
            let s = frac * tbl.total_length;
            tbl.s_to_t(s)
        } else {
            t0 + frac * (t1 - t0)
        };

        let pos = nurbs.point_at(t_param);
        let offset = pos - start_pos;
        let translation = Transform3::translation(offset.x, offset.y, offset.z);

        let transform = if align_to_curve {
            let current_tangent = nurbs.tangent_at(t_param);
            match alignment_rotation(start_tangent, current_tangent) {
                Some((axis, angle)) => {
                    let rotation =
                        Transform3::rotation_around_axis(pos, axis, angle).map_err(|e| {
                            ReplayError::FeatureFailed {
                                label: String::new(),
                                error: e.to_string(),
                            }
                        })?;
                    rotation.compose(&translation)
                }
                None => translation,
            }
        } else {
            translation
        };

        transforms.push(transform);
    }
    Ok(transforms)
}

/// Compute the (axis, angle) rotation that takes `from` to `to`.
///
/// Returns `None` when the vectors are (anti-)parallel or either is zero-length,
/// in which case no rotation is needed or determinable.
fn alignment_rotation(from: Vec3, to: Vec3) -> Option<(Vec3, f64)> {
    let cross = from.cross(&to);
    let sin_len = cross.length();
    if sin_len < 1e-10 {
        return None;
    }
    let dot = from.dot(&to).clamp(-1.0, 1.0);
    let angle = dot.acos();
    let axis = cross.normalized()?;
    Some((axis, angle))
}

/// Compute transforms for a mirror pattern.
///
/// Always returns exactly two transforms: [identity, reflection].
fn compute_mirror_transforms(
    plane: &MirrorPlaneDesc,
    brep: &BRep,
) -> Result<Vec<Transform3>, ReplayError> {
    let (origin, normal) = resolve_mirror_plane(plane, brep)?;
    let reflection =
        Transform3::reflection(origin, normal).map_err(|e| ReplayError::FeatureFailed {
            label: String::new(),
            error: e.to_string(),
        })?;
    Ok(vec![Transform3::identity(), reflection])
}

/// Resolve a [`MirrorPlaneDesc`] to an (origin, normal) pair.
fn resolve_mirror_plane(
    plane: &MirrorPlaneDesc,
    brep: &BRep,
) -> Result<(Point3, Vec3), ReplayError> {
    match plane {
        MirrorPlaneDesc::XY => Ok((Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))),
        MirrorPlaneDesc::XZ => Ok((Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0))),
        MirrorPlaneDesc::YZ => Ok((Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0))),
        MirrorPlaneDesc::OnFace(name) => {
            let face_handle = resolve::resolve_face(name, brep)?;
            let sof = brep.faces[face_handle].surface.as_ref().ok_or_else(|| {
                ReplayError::FeatureFailed {
                    label: String::new(),
                    error: "mirror face has no surface geometry".into(),
                }
            })?;

            // Compute face centroid from outer loop vertices.
            let positions = brep.face_outer_loop_positions(face_handle);
            if positions.is_empty() {
                return Err(ReplayError::FeatureFailed {
                    label: String::new(),
                    error: "mirror face has no vertices".into(),
                });
            }
            let n = positions.len() as f64;
            let cx = positions.iter().map(|p| p.x).sum::<f64>() / n;
            let cy = positions.iter().map(|p| p.y).sum::<f64>() / n;
            let cz = positions.iter().map(|p| p.z).sum::<f64>() / n;
            let centroid = Point3::new(cx, cy, cz);

            // Determine normal from the surface geometry.
            use popkern_core::geometry::surface::SurfaceEval;
            let normal = match &sof.geometry {
                popkern_core::geometry::SurfaceGeometry::Plane { normal, .. } => *normal,
                other => {
                    let (u_range, v_range) =
                        popkern_core::geometry::uv_projection::tight_bounded_parameter_range(
                            &positions, other,
                        );
                    let u_mid = (u_range.0 + u_range.1) / 2.0;
                    let v_mid = (v_range.0 + v_range.1) / 2.0;
                    other.normal_at(u_mid, v_mid)
                }
            };
            let normal = if sof.orientation { normal } else { -normal };

            Ok((centroid, normal))
        }
        MirrorPlaneDesc::Custom { origin, normal } => Ok((
            Point3::new(origin[0], origin[1], origin[2]),
            Vec3::new(normal[0], normal[1], normal[2]),
        )),
    }
}

#[cfg(test)]
#[path = "pattern_transform_tests.rs"]
mod tests;
