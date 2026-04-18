//! Face measurement functions.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::tight_bounded_parameter_range;
use popkern_core::geometry::{CurveEval, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_mesh::TessellationOptions;

use crate::MeasureError;
use crate::face_helpers::{loop_positions, pcurve_uv_range};

// ---------------------------------------------------------------------------
// Private geometry utilities
// ---------------------------------------------------------------------------

/// Compute the area of a planar polygon given its 3-D vertex positions using
/// the generalised Shoelace formula (cross-product method).
fn planar_polygon_area(positions: &[Point3], normal: Vec3) -> f64 {
    if positions.len() < 3 {
        return 0.0;
    }
    let origin = positions[0];
    let mut sum = Vec3::zero();
    for i in 1..positions.len() - 1 {
        let a = positions[i] - origin;
        let b = positions[i + 1] - origin;
        sum = Vec3::new(
            sum.x + a.cross(&b).x,
            sum.y + a.cross(&b).y,
            sum.z + a.cross(&b).z,
        );
    }
    0.5 * sum.dot(&normal).abs()
}

/// Compute face area by tessellating and summing triangle areas.
fn tessellate_area(
    brep: &BRep,
    face: FaceHandle,
    tolerance: &Tolerance,
) -> Result<f64, MeasureError> {
    let opts = TessellationOptions {
        chord_height: tolerance.linear().max(1e-4),
        ..TessellationOptions::default()
    };
    let face_mesh = popkern_mesh::tessellate_face(brep, face, &opts)
        .map_err(|e| MeasureError::TessellationFailed(e.to_string()))?;
    let area: f64 = face_mesh
        .triangles
        .iter()
        .map(|tri| {
            triangle_area(
                face_mesh.vertices[tri[0]].position,
                face_mesh.vertices[tri[1]].position,
                face_mesh.vertices[tri[2]].position,
            )
        })
        .sum();
    Ok(area)
}

/// Triangle area from three vertices.
fn triangle_area(a: Point3, b: Point3, c: Point3) -> f64 {
    0.5 * (b - a).cross(&(c - a)).length()
}

// ---------------------------------------------------------------------------
// Public measurement functions
// ---------------------------------------------------------------------------

/// Compute the surface area of a face.
///
/// Exact analytic formulas are used for `Plane`, `Cylinder`, `Sphere`,
/// `Cone`, and `Torus` faces that have no inner-loop holes. For `Nurbs` and
/// any curved face that contains holes, the face is tessellated and triangle
/// areas are summed.
pub fn face_area(
    brep: &BRep,
    face: FaceHandle,
    tolerance: &Tolerance,
) -> Result<f64, MeasureError> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| MeasureError::InvalidHandle("face".into()))?;
    let surface_on_face = face_data
        .surface
        .as_ref()
        .ok_or(MeasureError::MissingSurface)?;

    match &surface_on_face.geometry {
        SurfaceGeometry::Plane { normal, .. } => {
            let n = *normal;
            let outer_positions = loop_positions(brep, face_data.outer_loop);
            let outer_area = planar_polygon_area(&outer_positions, n);
            let hole_area: f64 = face_data
                .inner_loops
                .iter()
                .map(|&lh| planar_polygon_area(&loop_positions(brep, lh), n))
                .sum();
            Ok(outer_area - hole_area)
        }
        // Cylinder: area element = radius * du * dv.
        SurfaceGeometry::Cylinder { radius, .. } if face_data.inner_loops.is_empty() => {
            if let Some(((u0, u1), (v0, v1))) = pcurve_uv_range(brep, face) {
                Ok(radius * (u1 - u0).abs() * (v1 - v0).abs())
            } else {
                tessellate_area(brep, face, tolerance)
            }
        }
        // Sphere: area element = r^2*sin(v) du dv -> r^2*|du|*|cos(v0)-cos(v1)|.
        SurfaceGeometry::Sphere { radius, .. } if face_data.inner_loops.is_empty() => {
            if let Some(((u0, u1), (v0, v1))) = pcurve_uv_range(brep, face) {
                Ok(radius * radius * (u1 - u0).abs() * (v0.cos() - v1.cos()).abs())
            } else {
                tessellate_area(brep, face, tolerance)
            }
        }
        // Cone: area element = v*tan(ha)/cos(ha) du dv.
        // Integrated: tan(ha)/cos(ha) * |du| * (v1^2-v0^2)/2.
        SurfaceGeometry::Cone { half_angle, .. } if face_data.inner_loops.is_empty() => {
            if let Some(((u0, u1), (v0, v1))) = pcurve_uv_range(brep, face) {
                let delta_v_sq = (v1 * v1 - v0 * v0).abs() / 2.0;
                Ok(half_angle.tan() / half_angle.cos() * (u1 - u0).abs() * delta_v_sq)
            } else {
                tessellate_area(brep, face, tolerance)
            }
        }
        // Torus: area element = r*(R + r*cos(v)) du dv.
        // Integrated: r * |du| * |R*(v1-v0) + r*(sin(v1)-sin(v0))|.
        // Only valid for rectangular UV patches (|du| <= pi, |dv| <= pi);
        // seam-closing faces that span > half the period use tessellation.
        SurfaceGeometry::Torus {
            major_radius,
            minor_radius,
            ..
        } if face_data.inner_loops.is_empty() => {
            if let Some(((u0, u1), (v0, v1))) = pcurve_uv_range(brep, face) {
                let du = (u1 - u0).abs();
                let dv = (v1 - v0).abs();
                // Guard: seam faces get inflated UV bbox (> half-period).
                if du <= std::f64::consts::PI && dv <= std::f64::consts::PI {
                    let r = *minor_radius;
                    let big_r = *major_radius;
                    let integral = (big_r * dv + r * (v1.sin() - v0.sin())).abs();
                    Ok(r * du * integral)
                } else {
                    tessellate_area(brep, face, tolerance)
                }
            } else {
                tessellate_area(brep, face, tolerance)
            }
        }
        // §5.3 item #9 — Gregory-native surface-integral quadrature.
        // 8×8 Gauss-Legendre on the pcurve UV rectangle integrates
        // `|∂u × ∂v|` to ~1e-6 relative on smooth Gregory patches —
        // orders of magnitude tighter than tessellation. Falls back
        // to tessellation when the face has inner loops (trim-aware
        // integration is a separate rabbit hole) or the pcurve UV
        // range can't be determined.
        SurfaceGeometry::Gregory(_) if face_data.inner_loops.is_empty() => {
            if let Some(((u0, u1), (v0, v1))) = pcurve_uv_range(brep, face) {
                Ok(gauss_legendre_surface_area(
                    &surface_on_face.geometry,
                    (u0, u1),
                    (v0, v1),
                ))
            } else {
                tessellate_area(brep, face, tolerance)
            }
        }
        _ => tessellate_area(brep, face, tolerance),
    }
}

/// 8-point tensor-product Gauss-Legendre quadrature of
/// `∫∫ |∂u × ∂v| du dv` over a UV rectangle. Polynomial-exact to
/// degree 15 per axis; on the Gregory sqrt integrand the 8×8 grid
/// gives ~1e-6 relative accuracy — far tighter than tessellation-
/// based area.
fn gauss_legendre_surface_area(
    geom: &SurfaceGeometry,
    (u0, u1): (f64, f64),
    (v0, v1): (f64, f64),
) -> f64 {
    // 8-point Gauss-Legendre nodes + weights on [-1, 1].
    const GL_X: [f64; 8] = [
        -0.960_289_856_497_536_2,
        -0.796_666_477_413_626_7,
        -0.525_532_409_916_329,
        -0.183_434_642_495_649_8,
        0.183_434_642_495_649_8,
        0.525_532_409_916_329,
        0.796_666_477_413_626_7,
        0.960_289_856_497_536_2,
    ];
    const GL_W: [f64; 8] = [
        0.101_228_536_290_376_26,
        0.222_381_034_453_374_4,
        0.313_706_645_877_887_3,
        0.362_683_783_378_362,
        0.362_683_783_378_362,
        0.313_706_645_877_887_3,
        0.222_381_034_453_374_4,
        0.101_228_536_290_376_26,
    ];
    let du = 0.5 * (u1 - u0);
    let dv = 0.5 * (v1 - v0);
    let mu = 0.5 * (u1 + u0);
    let mv = 0.5 * (v1 + v0);
    let mut sum = 0.0_f64;
    for i in 0..8 {
        let u = mu + du * GL_X[i];
        for j in 0..8 {
            let v = mv + dv * GL_X[j];
            let su = geom.partial_u(u, v);
            let sv = geom.partial_v(u, v);
            sum += GL_W[i] * GL_W[j] * su.cross(&sv).length();
        }
    }
    du.abs() * dv.abs() * sum
}

/// Compute the perimeter of a face's outer boundary loop.
///
/// Sums the arc lengths of every edge in the outer loop.
pub fn face_perimeter(brep: &BRep, face: FaceHandle) -> Result<f64, MeasureError> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| MeasureError::InvalidHandle("face".into()))?;
    let loop_ = brep
        .loops
        .get(face_data.outer_loop)
        .ok_or_else(|| MeasureError::InvalidHandle("outer loop".into()))?;

    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Err(MeasureError::InvalidHandle("outer loop half-edge".into()));
    }

    let mut total = 0.0_f64;
    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        let edge = brep
            .edges
            .get(he.edge)
            .ok_or_else(|| MeasureError::InvalidHandle("edge".into()))?;
        let curve_on_edge = edge.curve.as_ref().ok_or(MeasureError::MissingCurve)?;
        total += curve_on_edge
            .geometry
            .length_between(curve_on_edge.t_start, curve_on_edge.t_end);
        current = he.next;
        if current == start_he {
            break;
        }
    }
    Ok(total)
}

/// Compute the outward unit normal of a face.
///
/// For planar faces the normal is read directly from the surface definition and
/// flipped if the surface orientation is reversed. For all other surface types
/// the normal is evaluated at the UV centroid of the face's parameter domain.
pub fn face_normal(brep: &BRep, face: FaceHandle) -> Result<Vec3, MeasureError> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| MeasureError::InvalidHandle("face".into()))?;
    let surface_on_face = face_data
        .surface
        .as_ref()
        .ok_or(MeasureError::MissingSurface)?;

    let raw_normal = match &surface_on_face.geometry {
        SurfaceGeometry::Plane { normal, .. } => *normal,
        geom => {
            let positions = brep.face_outer_loop_positions(face);
            let ((u0, u1), (v0, v1)) = tight_bounded_parameter_range(&positions, geom);
            let u_mid = 0.5 * (u0 + u1);
            let v_mid = 0.5 * (v0 + v1);
            geom.normal_at(u_mid, v_mid)
        }
    };

    let oriented = if surface_on_face.orientation {
        raw_normal
    } else {
        Vec3::new(-raw_normal.x, -raw_normal.y, -raw_normal.z)
    };

    oriented
        .normalized()
        .ok_or_else(|| MeasureError::DegenerateGeometry("zero-length face normal".into()))
}

/// Compute the area-weighted centroid of a face.
///
/// Tessellates the BRep, filters triangles belonging to this face, then
/// returns the area-weighted average of triangle centroids.
pub fn face_centroid(
    brep: &BRep,
    face: FaceHandle,
    tolerance: &Tolerance,
) -> Result<Point3, MeasureError> {
    brep.faces
        .get(face)
        .ok_or_else(|| MeasureError::InvalidHandle("face".into()))?;

    let opts = TessellationOptions {
        chord_height: tolerance.linear().max(1e-4),
        ..TessellationOptions::default()
    };
    let face_mesh = popkern_mesh::tessellate_face(brep, face, &opts)
        .map_err(|e| MeasureError::TessellationFailed(e.to_string()))?;

    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    let mut cz = 0.0_f64;
    let mut total_area = 0.0_f64;

    for tri in &face_mesh.triangles {
        let a = face_mesh.vertices[tri[0]].position;
        let b = face_mesh.vertices[tri[1]].position;
        let c = face_mesh.vertices[tri[2]].position;
        let area = triangle_area(a, b, c);
        cx += (a.x + b.x + c.x) / 3.0 * area;
        cy += (a.y + b.y + c.y) / 3.0 * area;
        cz += (a.z + b.z + c.z) / 3.0 * area;
        total_area += area;
    }

    if total_area < f64::EPSILON {
        return Err(MeasureError::DegenerateGeometry(
            "face has zero tessellated area".into(),
        ));
    }

    Ok(Point3::new(
        cx / total_area,
        cy / total_area,
        cz / total_area,
    ))
}

/// Return the characteristic radius of a face's surface.
///
/// - `Cylinder` -> cylinder radius
/// - `Sphere` -> sphere radius
/// - `Torus` -> minor (tube) radius
/// - All other surface types -> [`MeasureError::IncompatibleSurface`]
pub fn face_radius(brep: &BRep, face: FaceHandle) -> Result<f64, MeasureError> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| MeasureError::InvalidHandle("face".into()))?;
    let surface_on_face = face_data
        .surface
        .as_ref()
        .ok_or(MeasureError::MissingSurface)?;

    match &surface_on_face.geometry {
        SurfaceGeometry::Cylinder { radius, .. } => Ok(*radius),
        SurfaceGeometry::Sphere { radius, .. } => Ok(*radius),
        SurfaceGeometry::Torus { minor_radius, .. } => Ok(*minor_radius),
        other => Err(MeasureError::IncompatibleSurface {
            expected: "Cylinder, Sphere, or Torus",
            got: surface_type_name(other).to_string(),
        }),
    }
}

/// Return the principal curvatures `(k_min, k_max)` of a face's surface.
///
/// Analytic values are used for `Plane`, `Cylinder`, and `Sphere`. For all
/// other surfaces the curvatures are approximated from the first and second
/// fundamental form evaluated at the UV centroid.
pub fn face_curvature(brep: &BRep, face: FaceHandle) -> Result<(f64, f64), MeasureError> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| MeasureError::InvalidHandle("face".into()))?;
    let surface_on_face = face_data
        .surface
        .as_ref()
        .ok_or(MeasureError::MissingSurface)?;

    match &surface_on_face.geometry {
        SurfaceGeometry::Plane { .. } => Ok((0.0, 0.0)),
        SurfaceGeometry::Cylinder { radius, .. } => Ok((0.0, 1.0 / radius)),
        SurfaceGeometry::Sphere { radius, .. } => {
            let k = 1.0 / radius;
            Ok((k, k))
        }
        SurfaceGeometry::Cone { .. } => Ok((0.0, 0.0)),
        geom => {
            let positions = brep.face_outer_loop_positions(face);
            let ((u0, u1), (v0, v1)) = tight_bounded_parameter_range(&positions, geom);
            let u_mid = 0.5 * (u0 + u1);
            let v_mid = 0.5 * (v0 + v1);

            let su = geom.partial_u(u_mid, v_mid);
            let sv = geom.partial_v(u_mid, v_mid);

            // First fundamental form coefficients.
            let e = su.dot(&su);
            let f = su.dot(&sv);
            let g = sv.dot(&sv);
            let denom = e * g - f * f;
            if denom.abs() < f64::EPSILON {
                return Ok((0.0, 0.0));
            }

            // Second fundamental form via finite-difference of the normal field.
            let h = (u1 - u0).max(v1 - v0) * 1e-5;
            let normal_center = geom.normal_at(u_mid, v_mid);
            let normal_u = geom.normal_at(u_mid + h, v_mid);
            let normal_v = geom.normal_at(u_mid, v_mid + h);
            let dn_u = Vec3::new(
                (normal_u.x - normal_center.x) / h,
                (normal_u.y - normal_center.y) / h,
                (normal_u.z - normal_center.z) / h,
            );
            let dn_v = Vec3::new(
                (normal_v.x - normal_center.x) / h,
                (normal_v.y - normal_center.y) / h,
                (normal_v.z - normal_center.z) / h,
            );

            let l = -su.dot(&dn_u);
            let m = -su.dot(&dn_v);
            let n_coeff = -sv.dot(&dn_v);

            // Mean and Gaussian curvatures -> principal curvatures.
            let mean_k = (e * n_coeff - 2.0 * f * m + g * l) / (2.0 * denom);
            let gauss_k = (l * n_coeff - m * m) / denom;
            let sqrt_d = (mean_k * mean_k - gauss_k).max(0.0).sqrt();
            Ok((mean_k - sqrt_d, mean_k + sqrt_d))
        }
    }
}

// ---------------------------------------------------------------------------
// Private utility
// ---------------------------------------------------------------------------

/// Return a human-readable name for a surface geometry variant.
fn surface_type_name(geom: &SurfaceGeometry) -> &'static str {
    match geom {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cylinder",
        SurfaceGeometry::Sphere { .. } => "Sphere",
        SurfaceGeometry::Cone { .. } => "Cone",
        SurfaceGeometry::Torus { .. } => "Torus",
        SurfaceGeometry::Nurbs(_) => "Nurbs",
        SurfaceGeometry::Gregory(_) => "Gregory",
    }
}

#[cfg(test)]
#[path = "face_tests.rs"]
mod tests;
