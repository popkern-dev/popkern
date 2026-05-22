use popkern_core::geometry::curve_inversion::invert_curve_at_point;
use popkern_core::geometry::edge_triple::{EdgeAnchors, EdgeTripleOptions, compute_edge_geometry};
use popkern_core::geometry::fitting::interpolate_curve;
use popkern_core::geometry::intersection::{SurfaceSurfaceIntersection, intersect_surfaces};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::super::error::ApexBuildError;
use super::partition_cylinder_geometry::{CylinderTrimSurface, point_lies_on_cylinder_face};
use super::partition_support_intersections::line_segment_lies_on_cylinder_and_plane;

pub(super) fn exact_cylinder_trim_curve_between(
    brep: &BRep,
    cylinder_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
) -> Result<CurveOnEdge, ApexBuildError> {
    let start_point = brep
        .vertices
        .get(start)
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                 cylinder trim {label}: stale start vertex {start:?}"
            ),
        })?;
    let end_point = brep
        .vertices
        .get(end)
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                 cylinder trim {label}: stale end vertex {end:?}"
            ),
        })?;
    let cylinder = CylinderTrimSurface::from_face(brep, cylinder_face, fillet_index)?;
    if cylinder
        .validate_generator_segment(brep, fillet_index, start, end)
        .is_ok()
    {
        return CurveOnEdge::line_segment(start_point, end_point, brep.tolerance.linear())
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                     generator trim {label}: degenerate endpoints {start:?}->{end:?}"
                ),
            });
    }
    exact_support_cylinder_intersection_curve_between(
        brep,
        cylinder_face,
        start,
        end,
        start_point,
        end_point,
        label,
        fillet_index,
    )
}

pub(super) fn exact_cylinder_arc_trim_curve_between(
    brep: &BRep,
    cylinder_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
) -> Result<CurveOnEdge, ApexBuildError> {
    let start_point = brep.vertices[start].position;
    let end_point = brep.vertices[end].position;
    let Some(SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    }) = brep
        .faces
        .get(cylinder_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                 cylinder arc {label}: target face {cylinder_face:?} is not a cylinder"
            ),
        });
    };
    let Some(axis_unit) = axis.normalized() else {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                 cylinder arc {label}: degenerate cylinder axis"
            ),
        });
    };
    let start_axis = (start_point - *origin).dot(&axis_unit);
    let end_axis = (end_point - *origin).dot(&axis_unit);
    let tol = brep.tolerance.linear().max(1.0e-9);
    if !point_lies_on_cylinder_face(brep, cylinder_face, start_point)
        || !point_lies_on_cylinder_face(brep, cylinder_face, end_point)
    {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                 cylinder arc {label} {start:?}->{end:?}: endpoint is off cylinder"
            ),
        });
    }
    if (start_axis - end_axis).abs() > tol {
        return cylinder_uv_line_trim_curve_between(
            brep,
            *origin,
            axis_unit,
            *radius,
            start_axis,
            end_axis,
            start_point,
            end_point,
            label,
            fillet_index,
        );
    }
    let center = *origin + axis_unit * start_axis;
    let (u_axis, v_axis) = circle_basis_for_normal(axis_unit);
    let start_radial = start_point - center;
    let end_radial = end_point - center;
    let start_angle = start_radial.dot(&v_axis).atan2(start_radial.dot(&u_axis));
    let mut end_angle = end_radial.dot(&v_axis).atan2(end_radial.dot(&u_axis));
    while end_angle - start_angle > std::f64::consts::PI {
        end_angle -= std::f64::consts::TAU;
    }
    while end_angle - start_angle < -std::f64::consts::PI {
        end_angle += std::f64::consts::TAU;
    }
    Ok(CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center,
            normal: axis_unit,
            radius: *radius,
            start_angle,
            end_angle,
        },
        t_start: start_angle,
        t_end: end_angle,
    })
}

#[allow(clippy::too_many_arguments)]
fn cylinder_uv_line_trim_curve_between(
    brep: &BRep,
    origin: Point3,
    axis_unit: Vec3,
    radius: f64,
    start_axis: f64,
    end_axis: f64,
    start_point: Point3,
    end_point: Point3,
    label: &'static str,
    fillet_index: usize,
) -> Result<CurveOnEdge, ApexBuildError> {
    let (u_axis, v_axis) = circle_basis_for_normal(axis_unit);
    let start_center = origin + axis_unit * start_axis;
    let end_center = origin + axis_unit * end_axis;
    let start_radial = start_point - start_center;
    let end_radial = end_point - end_center;
    let start_angle = start_radial.dot(&v_axis).atan2(start_radial.dot(&u_axis));
    let mut end_angle = end_radial.dot(&v_axis).atan2(end_radial.dot(&u_axis));
    while end_angle - start_angle > std::f64::consts::PI {
        end_angle -= std::f64::consts::TAU;
    }
    while end_angle - start_angle < -std::f64::consts::PI {
        end_angle += std::f64::consts::TAU;
    }

    let sample_count = 17;
    let mut points = Vec::with_capacity(sample_count);
    let mut params = Vec::with_capacity(sample_count);
    for i in 0..sample_count {
        let t = i as f64 / (sample_count as f64 - 1.0);
        let angle = start_angle + t * (end_angle - start_angle);
        let station = start_axis + t * (end_axis - start_axis);
        let center = origin + axis_unit * station;
        points.push(center + u_axis * (radius * angle.cos()) + v_axis * (radius * angle.sin()));
        params.push(t);
    }
    let nurbs =
        interpolate_curve(&points, &params, 5).map_err(|err| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                 cylinder UV-line trim {label}: interpolation failed: {err:?}"
            ),
        })?;
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs),
        t_start: 0.0,
        t_end: 1.0,
    };
    let tolerance = brep.tolerance.linear().max(1.0e-7);
    for i in 0..=32 {
        let t = i as f64 / 32.0;
        let point = curve.geometry.point_at(t);
        let station = (point - origin).dot(&axis_unit);
        let center = origin + axis_unit * station;
        let radial = point.distance_to(&center);
        if (radial - radius).abs() > tolerance {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cylinder UV-line \
                     trim {label} exceeded cylinder tolerance at sample {i}: radial deviation \
                     {:.3e}",
                    (radial - radius).abs()
                ),
            });
        }
    }
    Ok(curve)
}

fn circle_basis_for_normal(normal: Vec3) -> (Vec3, Vec3) {
    let axis = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = normal
        .cross(&axis)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = normal
        .cross(&u)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn exact_support_cylinder_intersection_curve_between(
    brep: &BRep,
    cylinder_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    start_point: Point3,
    end_point: Point3,
    label: &'static str,
    fillet_index: usize,
) -> Result<CurveOnEdge, ApexBuildError> {
    let cylinder_surface = brep
        .faces
        .get(cylinder_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                 cylinder trim {label}: target cylinder face {cylinder_face:?} has no surface"
            ),
        })?;
    let mut blockers = Vec::new();
    for (support_face, face_data) in brep.faces.iter() {
        if support_face == cylinder_face {
            continue;
        }
        let Some(support_surface) = face_data.surface.as_ref().map(|surface| &surface.geometry)
        else {
            continue;
        };
        let SurfaceGeometry::Plane { origin, normal, .. } = support_surface else {
            continue;
        };
        let tol = brep.tolerance.linear().max(1.0e-9);
        let start_plane_distance = (start_point - *origin).dot(normal).abs();
        let end_plane_distance = (end_point - *origin).dot(normal).abs();
        if start_plane_distance > tol || end_plane_distance > tol {
            blockers.push(format!(
                "{support_face:?}: plane distances {:.3e}/{:.3e}",
                start_plane_distance, end_plane_distance
            ));
            continue;
        }
        if line_segment_lies_on_cylinder_and_plane(
            brep,
            cylinder_face,
            start_point,
            end_point,
            *origin,
            *normal,
        ) {
            return CurveOnEdge::line_segment(start_point, end_point, brep.tolerance.linear())
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                         line support/cylinder trim for {label}: degenerate endpoints \
                         {start:?}->{end:?}"
                    ),
                });
        }
        let anchors = EdgeAnchors::Open {
            start: start_point,
            end: end_point,
        };
        if let Ok(triple) = compute_edge_geometry(
            cylinder_surface,
            support_surface,
            anchors,
            None,
            &brep.tolerance,
            &EdgeTripleOptions::default(),
        ) {
            let trim_start = triple.curve_3d.point_at(triple.t_start);
            let trim_end = triple.curve_3d.point_at(triple.t_end);
            if trim_start
                .distance_to(&start_point)
                .max(trim_end.distance_to(&end_point))
                <= brep.tolerance.linear().max(1.0e-9)
            {
                return Ok(CurveOnEdge {
                    geometry: triple.curve_3d,
                    t_start: triple.t_start,
                    t_end: triple.t_end,
                });
            }
            blockers.push(format!(
                "{support_face:?}: edge-triple endpoint mismatch {:.3e}/{:.3e}",
                trim_start.distance_to(&start_point),
                trim_end.distance_to(&end_point)
            ));
            continue;
        }
        let result = intersect_surfaces(
            cylinder_surface,
            (0.0, std::f64::consts::TAU),
            (-10.0, 10.0),
            support_surface,
            (-10.0, 10.0),
            (-10.0, 10.0),
            &brep.tolerance,
        )
        .map_err(|err| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] failed analytic \
                 support/cylinder intersection for {label} on {support_face:?}: {err:?}"
            ),
        })?;
        for intersection in result.intersections {
            let SurfaceSurfaceIntersection::Curve(curve_on_surface) = intersection else {
                continue;
            };
            let range = curve_on_surface.curve_3d.parameter_range();
            let candidate = CurveOnEdge {
                geometry: CurveGeometry::Nurbs(curve_on_surface.curve_3d.clone()),
                t_start: range.0,
                t_end: range.1,
            };
            let Ok(t_start) = invert_curve_at_point(&candidate, start_point, &brep.tolerance)
            else {
                continue;
            };
            let Ok(t_end) = invert_curve_at_point(&candidate, end_point, &brep.tolerance) else {
                continue;
            };
            if (t_start - t_end).abs() <= brep.tolerance.linear().max(1.0e-9) {
                continue;
            }
            return Ok(CurveOnEdge {
                geometry: candidate.geometry,
                t_start,
                t_end,
            });
        }
        blockers.push(format!(
            "{support_face:?}: analytic intersection did not invert both endpoints {:?}->{:?}",
            start, end
        ));
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author exact \
             support/cylinder curve for {label} {start:?}->{end:?}; tried support planes [{}]",
            blockers.join("; ")
        ),
    })
}
