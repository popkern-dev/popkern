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
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;

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
/// Returns a [`Vec`] of [`Transform3`] values, one per instance.  Instance 0
/// is always [`Transform3::identity`] (the seed itself).
///
/// # Errors
///
/// Returns [`ReplayError`] if entity resolution, axis computation, or
/// parameter lookup fails, or if `count == 0`.
pub(crate) fn compute_pattern_transforms(
    feature: &Feature,
    brep: &BRep,
    param_table: &ParameterTable,
) -> Result<Vec<Transform3>, ReplayError> {
    match feature {
        Feature::LinearPattern {
            primary, secondary, ..
        } => compute_linear_transforms(primary, secondary.as_ref(), brep, param_table),
        Feature::CircularPattern {
            axis,
            count,
            total_angle,
            ..
        } => compute_circular_transforms(axis, *count, *total_angle, brep, param_table),
        Feature::CurvePattern { curve, count, .. } => {
            compute_curve_transforms(curve, *count, brep, param_table)
        }
        Feature::Mirror { plane, .. } => compute_mirror_transforms(plane, brep),
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
fn compute_circular_transforms(
    axis: &AxisDesc,
    count_id: crate::parameter::ParameterId,
    total_angle_id: crate::parameter::ParameterId,
    brep: &BRep,
    param_table: &ParameterTable,
) -> Result<Vec<Transform3>, ReplayError> {
    let (origin, dir) = resolve::resolve_axis(axis, brep)?;
    let count = resolve_param(count_id, param_table)?.round() as usize;
    if count == 0 {
        return Err(ReplayError::FeatureFailed {
            label: String::new(),
            error: "circular pattern count must be > 0".into(),
        });
    }

    let raw_angle = resolve_param(total_angle_id, param_table)?;
    // TODO: spec says "clamp to 360° with warning" — a structured warning
    // channel doesn't exist yet.  For now we clamp silently.
    let total_angle = raw_angle.clamp(-std::f64::consts::TAU, std::f64::consts::TAU);

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
    Ok(transforms)
}

/// Compute transforms for a curve-driven pattern.
///
/// Resolves the named edge to extract its start and end vertices, then
/// distributes `count` instances evenly along the chord (start → end
/// direction) as pure translations.
///
/// TODO: replace with full arc-length distribution using
/// `popkern_sweep::arc_length::ArcLengthTable` once that crate is available
/// as a dependency.  Frame alignment (aligning each instance tangent to the
/// curve) also needs the arc-length reparameterisation.
fn compute_curve_transforms(
    curve: &popkern_core::naming::EntityName,
    count_id: crate::parameter::ParameterId,
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

    // Guard: only linear (straight) edges are supported in v1.
    let curve_on_edge =
        brep.edges[edge_handle]
            .curve
            .as_ref()
            .ok_or_else(|| ReplayError::FeatureFailed {
                label: String::new(),
                error: "curve pattern edge has no curve geometry".into(),
            })?;
    if !matches!(
        curve_on_edge.geometry,
        popkern_core::geometry::curve::CurveGeometry::Line { .. }
    ) {
        return Err(ReplayError::FeatureFailed {
            label: String::new(),
            error: "CurvePattern along non-linear edges is not yet implemented; use a linear edge"
                .into(),
        });
    }

    let (start_pos, end_pos) = edge_chord(brep, edge_handle)?;

    let chord = end_pos - start_pos;
    let chord_len = chord.length();

    let mut transforms = Vec::with_capacity(count);
    for i in 0..count {
        let t = if count > 1 && chord_len > f64::EPSILON {
            let frac = i as f64 / (count - 1) as f64;
            let offset = chord * frac;
            Transform3::translation(offset.x, offset.y, offset.z)
        } else {
            Transform3::identity()
        };
        transforms.push(t);
    }
    Ok(transforms)
}

/// Return the start and end world-space positions of an edge's chord.
fn edge_chord(
    brep: &BRep,
    edge_handle: popkern_core::topology::EdgeHandle,
) -> Result<(Point3, Point3), ReplayError> {
    let edge = &brep.edges[edge_handle];
    let he = &brep.half_edges[edge.half_edge];
    let twin = &brep.half_edges[he.twin];
    let start = brep.vertices[he.origin].position;
    let end = brep.vertices[twin.origin].position;
    Ok((start, end))
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
