//! Face measurement functions.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::tight_bounded_parameter_range;
use popkern_core::geometry::{CurveEval, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_mesh::TessellationOptions;

use crate::MeasureError;

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Collect 3-D vertex positions for every half-edge in a loop.
fn loop_positions(brep: &BRep, loop_handle: popkern_core::topology::LoopHandle) -> Vec<Point3> {
    let Some(loop_) = brep.loops.get(loop_handle) else {
        return Vec::new();
    };
    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }
    let mut positions = Vec::new();
    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        positions.push(brep.vertices[he.origin].position);
        current = he.next;
        if current == start_he {
            break;
        }
    }
    positions
}

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
    // Project onto the face normal to get the signed area.
    0.5 * sum.dot(&normal).abs()
}

/// Compute the area contribution of a triangulated face by summing triangle areas.
fn triangle_area(a: Point3, b: Point3, c: Point3) -> f64 {
    let ab = b - a;
    let ac = c - a;
    0.5 * ab.cross(&ac).length()
}

// ---------------------------------------------------------------------------
// Public measurement functions
// ---------------------------------------------------------------------------

/// Compute the surface area of a face.
///
/// For planar faces the exact Shoelace formula is used (outer loop minus inner
/// loop holes). For all other surface types the face is tessellated and the
/// triangle areas are summed.
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
                .map(|&lh| {
                    let positions = loop_positions(brep, lh);
                    planar_polygon_area(&positions, n)
                })
                .sum();
            Ok(outer_area - hole_area)
        }
        _ => {
            let opts = TessellationOptions {
                chord_height: tolerance.linear.max(1e-4),
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
    }
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
        chord_height: tolerance.linear.max(1e-4),
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
/// - `Cylinder` → cylinder radius
/// - `Sphere` → sphere radius
/// - `Torus` → minor (tube) radius
/// - All other surface types → [`MeasureError::IncompatibleSurface`]
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
            // Approximate via shape operator: use surface partials and their
            // derivatives at the UV centroid.
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

            // Second fundamental form: finite-difference the normal field.
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

            // Mean and Gaussian curvatures → principal curvatures.
            let mean_k = (e * n_coeff - 2.0 * f * m + g * l) / (2.0 * denom);
            let gauss_k = (l * n_coeff - m * m) / denom;

            let discriminant = (mean_k * mean_k - gauss_k).max(0.0);
            let sqrt_d = discriminant.sqrt();
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
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn make_cube() -> BRep {
        BRep::unit_cube(Tolerance::default()).expect("cube")
    }

    #[test]
    fn unit_cube_face_area_is_one() {
        let cube = make_cube();
        let tol = Tolerance::default();
        let face = cube.faces.keys().next().unwrap();
        let area = face_area(&cube, face, &tol).unwrap();
        assert!((area - 1.0).abs() < 1e-10, "expected 1.0, got {area}");
    }

    #[test]
    fn unit_cube_face_perimeter_is_four() {
        let cube = make_cube();
        let face = cube.faces.keys().next().unwrap();
        let perim = face_perimeter(&cube, face).unwrap();
        assert!((perim - 4.0).abs() < 1e-10, "expected 4.0, got {perim}");
    }

    #[test]
    fn unit_cube_face_normal_is_unit() {
        let cube = make_cube();
        let face = cube.faces.keys().next().unwrap();
        let n = face_normal(&cube, face).unwrap();
        assert!((n.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn planar_face_radius_is_incompatible() {
        let cube = make_cube();
        let face = cube.faces.keys().next().unwrap();
        assert!(matches!(
            face_radius(&cube, face),
            Err(MeasureError::IncompatibleSurface { .. })
        ));
    }

    #[test]
    fn planar_face_curvature_is_zero() {
        let cube = make_cube();
        let face = cube.faces.keys().next().unwrap();
        let (k_min, k_max) = face_curvature(&cube, face).unwrap();
        assert!(k_min.abs() < 1e-10);
        assert!(k_max.abs() < 1e-10);
    }
}
