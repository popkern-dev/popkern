//! Read-only analytic three-face tangency probe for corner diagnostics.

#[path = "corner_apex_diagnostics_tangency_lines.rs"]
mod tangency_lines;

use super::live_helpers::*;
use crate::euler::DiagnosticTolerance;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::topology::{BRep, FaceHandle, VertexHandle};
use tangency_lines::{
    ConeEquation, CylinderEquation, PlaneEquation, line_cone_intersections,
    line_cylinder_intersections, line_sphere_intersections, plane_plane_offset_line,
};

#[derive(Clone)]
pub(super) struct TangencySupport {
    surface: SurfaceGeometry,
}

impl TangencySupport {
    fn new(surface: SurfaceGeometry) -> Option<Self> {
        matches!(
            surface,
            SurfaceGeometry::Plane { .. }
                | SurfaceGeometry::Cylinder { .. }
                | SurfaceGeometry::Cone { .. }
                | SurfaceGeometry::Sphere { .. }
        )
        .then_some(Self { surface })
    }

    fn signed_distance_and_normal(&self, point: Point3) -> Option<(f64, Vec3)> {
        match self.surface {
            SurfaceGeometry::Plane { origin, normal, .. } => {
                let normal = normal.normalized()?;
                Some(((point - origin).dot(&normal), normal))
            }
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } => {
                let axis = axis.normalized()?;
                let axial = axis * ((point - origin).dot(&axis));
                let radial = (point - origin) - axial;
                let normal = radial.normalized()?;
                Some((radial.length() - radius, normal))
            }
            SurfaceGeometry::Sphere { center, radius } => {
                let radial = point - center;
                let normal = radial.normalized()?;
                Some((radial.length() - radius, normal))
            }
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            } => {
                let axis = axis.normalized()?;
                let from_apex = point - apex;
                let axial_distance = from_apex.dot(&axis);
                let radial = from_apex - axis * axial_distance;
                let radial_direction = radial.normalized()?;
                let normal = radial_direction * half_angle.cos() - axis * half_angle.sin();
                Some((
                    radial.length() * half_angle.cos() - axial_distance * half_angle.sin(),
                    normal,
                ))
            }
            _ => None,
        }
    }

    fn normal_near(&self, point: Point3) -> Option<Vec3> {
        self.signed_distance_and_normal(point)
            .and_then(|(_, normal)| normal.normalized())
    }
}

pub(super) fn tangency_supports(
    brep: &BRep,
    vertex: VertexHandle,
) -> Result<Vec<TangencySupport>, &'static str> {
    let faces = incident_faces_at_vertex(brep, vertex);
    if faces.len() < 3 {
        return Err("three-face-tangency-support-count");
    }
    faces
        .into_iter()
        .take(3)
        .map(|face| tangency_support(brep, face))
        .collect()
}

fn tangency_support(brep: &BRep, face: FaceHandle) -> Result<TangencySupport, &'static str> {
    let Some(surface) = brep.faces.get(face).and_then(|face| face.surface.as_ref()) else {
        return Err("three-face-tangency-nonplanar-support");
    };
    TangencySupport::new(surface.geometry.clone()).ok_or("three-face-tangency-nonplanar-support")
}

pub(super) fn solve_three_face_tangency(
    supports: &[TangencySupport],
    vertex: Point3,
    radius: f64,
    tolerance: DiagnosticTolerance,
) -> Option<Point3> {
    let tolerance_m = tolerance_meters(tolerance);
    for signs in sign_sets() {
        for seed in tangency_seeds(supports, vertex, radius, signs) {
            let Some(center) = newton_solve_tangency(supports, seed, radius, signs, tolerance_m)
            else {
                continue;
            };
            if center_is_tangent(supports, center, radius, signs, tolerance_m)
                && center_has_stable_tangency_normals(supports, center)
            {
                return Some(center);
            }
        }
    }
    None
}

fn center_has_stable_tangency_normals(supports: &[TangencySupport], center: Point3) -> bool {
    supports
        .iter()
        .all(|support| support.normal_near(center).is_some())
}

fn tangency_seeds(
    supports: &[TangencySupport],
    vertex: Point3,
    radius: f64,
    signs: [f64; 3],
) -> Vec<Point3> {
    let normals = supports
        .iter()
        .filter_map(|support| support.normal_near(vertex))
        .collect::<Vec<_>>();
    if normals.len() != 3 {
        return vec![vertex];
    }
    let mut seeds = Vec::new();
    seeds.extend(exact_sphere_two_plane_seeds(supports, radius, signs));
    seeds.extend(exact_cylinder_two_plane_seeds(supports, radius, signs));
    seeds.extend(exact_cone_two_plane_seeds(supports, radius, signs));
    if let Some(seed) = linearized_tangent_seed(vertex, &normals, radius, signs) {
        seeds.push(seed);
    }
    let combined = normals
        .iter()
        .zip(signs)
        .fold(Vec3::zero(), |sum, (normal, sign)| sum + *normal * sign);
    if let Some(direction) = combined.normalized() {
        seeds.push(vertex + direction * radius);
        seeds.push(vertex + direction * (radius * 2.0));
    }
    seeds.push(
        vertex
            + normals[0] * (signs[0] * radius)
            + normals[1] * (signs[1] * radius)
            + normals[2] * (signs[2] * radius),
    );
    seeds.push(vertex);
    seeds
}

fn exact_cylinder_two_plane_seeds(
    supports: &[TangencySupport],
    radius: f64,
    signs: [f64; 3],
) -> Vec<Point3> {
    let mut seeds = Vec::new();
    for cylinder_index in 0..3 {
        let Some(cylinder) = supports[cylinder_index].cylinder() else {
            continue;
        };
        let planes = (0..3)
            .filter(|index| *index != cylinder_index)
            .collect::<Vec<_>>();
        let Some(first_plane) = supports[planes[0]].plane() else {
            continue;
        };
        let Some(second_plane) = supports[planes[1]].plane() else {
            continue;
        };
        let Some(line) = plane_plane_offset_line(
            first_plane,
            signs[planes[0]] * radius,
            second_plane,
            signs[planes[1]] * radius,
            cylinder.origin,
        ) else {
            continue;
        };
        let target_radius = cylinder.radius + signs[cylinder_index] * radius;
        if target_radius > 0.0 {
            seeds.extend(line_cylinder_intersections(line, cylinder, target_radius));
        }
    }
    seeds
}

fn exact_sphere_two_plane_seeds(
    supports: &[TangencySupport],
    radius: f64,
    signs: [f64; 3],
) -> Vec<Point3> {
    let mut seeds = Vec::new();
    for sphere_index in 0..3 {
        let Some((sphere_center, sphere_radius)) = supports[sphere_index].sphere() else {
            continue;
        };
        let planes = (0..3)
            .filter(|index| *index != sphere_index)
            .collect::<Vec<_>>();
        let Some(first_plane) = supports[planes[0]].plane() else {
            continue;
        };
        let Some(second_plane) = supports[planes[1]].plane() else {
            continue;
        };
        let Some(line) = plane_plane_offset_line(
            first_plane,
            signs[planes[0]] * radius,
            second_plane,
            signs[planes[1]] * radius,
            sphere_center,
        ) else {
            continue;
        };
        let target_radius = sphere_radius + signs[sphere_index] * radius;
        if target_radius > 0.0 {
            seeds.extend(line_sphere_intersections(
                line,
                sphere_center,
                target_radius,
            ));
        }
    }
    seeds
}

fn exact_cone_two_plane_seeds(
    supports: &[TangencySupport],
    radius: f64,
    signs: [f64; 3],
) -> Vec<Point3> {
    let mut seeds = Vec::new();
    for cone_index in 0..3 {
        let Some(cone) = supports[cone_index].cone() else {
            continue;
        };
        let planes = (0..3)
            .filter(|index| *index != cone_index)
            .collect::<Vec<_>>();
        let Some(first_plane) = supports[planes[0]].plane() else {
            continue;
        };
        let Some(second_plane) = supports[planes[1]].plane() else {
            continue;
        };
        let Some(line) = plane_plane_offset_line(
            first_plane,
            signs[planes[0]] * radius,
            second_plane,
            signs[planes[1]] * radius,
            cone.apex,
        ) else {
            continue;
        };
        seeds.extend(line_cone_intersections(
            line,
            cone,
            signs[cone_index] * radius,
        ));
    }
    seeds
}

fn linearized_tangent_seed(
    vertex: Point3,
    normals: &[Vec3],
    radius: f64,
    signs: [f64; 3],
) -> Option<Point3> {
    let vertex_vec = Vec3::new(vertex.x, vertex.y, vertex.z);
    let rhs = Vec3::new(
        normals[0].dot(&vertex_vec) + signs[0] * radius,
        normals[1].dot(&vertex_vec) + signs[1] * radius,
        normals[2].dot(&vertex_vec) + signs[2] * radius,
    );
    solve_3x3([normals[0], normals[1], normals[2]], rhs)
        .map(|seed| Point3::new(seed.x, seed.y, seed.z))
}

impl TangencySupport {
    fn plane(&self) -> Option<PlaneEquation> {
        let SurfaceGeometry::Plane { origin, normal, .. } = self.surface else {
            return None;
        };
        let normal = normal.normalized()?;
        let origin = Vec3::new(origin.x, origin.y, origin.z);
        Some(PlaneEquation {
            normal,
            offset: normal.dot(&origin),
        })
    }

    fn sphere(&self) -> Option<(Point3, f64)> {
        let SurfaceGeometry::Sphere { center, radius } = self.surface else {
            return None;
        };
        Some((center, radius))
    }

    fn cylinder(&self) -> Option<CylinderEquation> {
        let SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } = self.surface
        else {
            return None;
        };
        Some(CylinderEquation {
            origin,
            axis: axis.normalized()?,
            radius,
        })
    }

    fn cone(&self) -> Option<ConeEquation> {
        let SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } = self.surface
        else {
            return None;
        };
        Some(ConeEquation {
            apex,
            axis: axis.normalized()?,
            half_angle,
        })
    }
}

fn newton_solve_tangency(
    supports: &[TangencySupport],
    mut center: Point3,
    radius: f64,
    signs: [f64; 3],
    tolerance: f64,
) -> Option<Point3> {
    let step_limit = (radius * 32.0).max(1.0).max(tolerance * 16.0);
    for _ in 0..64 {
        let residual = tangency_residual(supports, center, radius, signs)?;
        if max_abs(residual) <= tolerance.max(radius * 1e-8) {
            return Some(center);
        }
        let jacobian = finite_difference_jacobian(supports, center, radius, signs)?;
        let mut delta = solve_3x3(
            jacobian,
            Vec3::new(-residual[0], -residual[1], -residual[2]),
        )?;
        let length = delta.length();
        if !length.is_finite() {
            return None;
        }
        if length > step_limit {
            delta = delta * (step_limit / length);
        }
        center = center + delta;
    }
    None
}

fn tangency_residual(
    supports: &[TangencySupport],
    center: Point3,
    radius: f64,
    signs: [f64; 3],
) -> Option<[f64; 3]> {
    if supports.len() != 3 {
        return None;
    }
    Some([
        supports[0].signed_distance_and_normal(center)?.0 - signs[0] * radius,
        supports[1].signed_distance_and_normal(center)?.0 - signs[1] * radius,
        supports[2].signed_distance_and_normal(center)?.0 - signs[2] * radius,
    ])
}

fn finite_difference_jacobian(
    supports: &[TangencySupport],
    center: Point3,
    radius: f64,
    signs: [f64; 3],
) -> Option<[Vec3; 3]> {
    let h = radius.max(1.0) * 1e-6;
    let x = finite_difference_column(supports, center, radius, signs, Vec3::new(h, 0.0, 0.0), h)?;
    let y = finite_difference_column(supports, center, radius, signs, Vec3::new(0.0, h, 0.0), h)?;
    let z = finite_difference_column(supports, center, radius, signs, Vec3::new(0.0, 0.0, h), h)?;
    Some([
        Vec3::new(x[0], y[0], z[0]),
        Vec3::new(x[1], y[1], z[1]),
        Vec3::new(x[2], y[2], z[2]),
    ])
}

fn finite_difference_column(
    supports: &[TangencySupport],
    center: Point3,
    radius: f64,
    signs: [f64; 3],
    step: Vec3,
    h: f64,
) -> Option<[f64; 3]> {
    let forward = tangency_residual(supports, center + step, radius, signs)?;
    let backward = tangency_residual(supports, center - step, radius, signs)?;
    Some([
        (forward[0] - backward[0]) / (2.0 * h),
        (forward[1] - backward[1]) / (2.0 * h),
        (forward[2] - backward[2]) / (2.0 * h),
    ])
}

fn center_is_tangent(
    supports: &[TangencySupport],
    center: Point3,
    radius: f64,
    signs: [f64; 3],
    tolerance: f64,
) -> bool {
    tangency_residual(supports, center, radius, signs)
        .is_some_and(|residual| max_abs(residual) <= tolerance.max(radius * 1e-8))
}

fn solve_3x3(rows: [Vec3; 3], rhs: Vec3) -> Option<Vec3> {
    let det = rows[0].dot(&rows[1].cross(&rows[2]));
    if det.abs() <= 1e-12 {
        return None;
    }
    let x = rhs.dot(&rows[1].cross(&rows[2])) / det;
    let y = rows[0].dot(&rhs.cross(&rows[2])) / det;
    let z = rows[0].dot(&rows[1].cross(&rhs)) / det;
    (x.is_finite() && y.is_finite() && z.is_finite()).then_some(Vec3::new(x, y, z))
}

fn max_abs(values: [f64; 3]) -> f64 {
    values
        .into_iter()
        .map(f64::abs)
        .fold(0.0, |max, value| max.max(value))
}

fn sign_sets() -> [[f64; 3]; 8] {
    [
        [-1.0, -1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, -1.0],
        [1.0, 1.0, 1.0],
    ]
}

fn tolerance_meters(tolerance: DiagnosticTolerance) -> f64 {
    tolerance.nanometers.unwrap_or(10).max(1) as f64 * 1e-9
}
