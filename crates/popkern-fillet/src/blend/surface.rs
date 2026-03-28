//! Vertex blend surface construction (mitred and smooth modes).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::fillet_surface::FilletSurfaceResult;
use crate::types::VertexMode;

use super::boundary::{BlendBoundary, BoundarySegment};

/// Result of vertex blend computation.
#[derive(Debug, Clone)]
pub struct BlendSurfaceResult {
    /// The constructed blend surface.
    pub surface: SurfaceGeometry,
    /// Boundary curves of the blend region.
    pub boundary_curves: Vec<NurbsCurve>,
}

/// Build a vertex blend surface.
///
/// For [`VertexMode::Mitred`], fits a least-squares plane through the boundary
/// points and returns a planar surface.
///
/// For [`VertexMode::Smooth`] with two fillets, attempts a spherical cap
/// tangent to both fillet surfaces. If G1 continuity cannot be achieved,
/// falls back to NURBS patch fitting. For three or more fillets, degrades
/// gracefully to the mitred approach.
///
/// # Errors
///
/// Returns an error if the boundary has no segments or if the geometry is
/// degenerate (e.g., all boundary points coincide).
pub fn build_blend_surface(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
    vertex_mode: VertexMode,
    tolerance: &Tolerance,
) -> Result<BlendSurfaceResult, String> {
    if boundary.segments.is_empty() {
        return Err("blend boundary has no segments".into());
    }

    // Collect all boundary points.
    let boundary_curves = collect_boundary_curves(boundary);
    let points = sample_boundary_points(&boundary_curves, 10);

    if points.len() < 3 {
        return Err("need at least 3 boundary points for blend surface".into());
    }

    match vertex_mode {
        VertexMode::Mitred => build_mitred_blend(&points, boundary_curves),
        VertexMode::Smooth => build_smooth_blend(
            boundary,
            fillet_surfaces,
            &points,
            boundary_curves,
            tolerance,
            brep,
        ),
    }
}

/// Build a mitred (planar) blend surface.
///
/// Fits a least-squares plane through the boundary sample points.
fn build_mitred_blend(
    points: &[Point3],
    boundary_curves: Vec<NurbsCurve>,
) -> Result<BlendSurfaceResult, String> {
    let (origin, normal) = fit_plane_least_squares(points)?;

    // Construct orthonormal frame.
    let u_axis = pick_perpendicular(&normal);
    let v_axis = normal.cross(&u_axis);
    let v_axis = v_axis
        .normalized()
        .ok_or_else(|| "degenerate v_axis in blend plane".to_string())?;

    let surface = SurfaceGeometry::Plane {
        origin,
        normal,
        u_axis,
        v_axis,
    };

    Ok(BlendSurfaceResult {
        surface,
        boundary_curves,
    })
}

/// Build a smooth blend surface.
///
/// Fallback chain: sphere → polar NURBS → mitred plane.
///
/// For two-fillet corners, the torus axes define a valid sphere bisector,
/// so the sphere is tried first. If it fails, polar NURBS is tried as an
/// intermediate fallback before the mitred plane.
///
/// For three or more fillets (e.g., a cube corner), the geometry is more
/// complex. Polar NURBS is tried first; the mitred plane is the final
/// fallback.
fn build_smooth_blend(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    points: &[Point3],
    boundary_curves: Vec<NurbsCurve>,
    tolerance: &Tolerance,
    brep: &BRep,
) -> Result<BlendSurfaceResult, String> {
    let n_fillets = boundary
        .segments
        .iter()
        .filter(|s| matches!(s, BoundarySegment::FilletEnd { .. }))
        .count();

    // Try spherical cap first.
    if let Some(result) = try_spherical_blend(
        boundary,
        fillet_surfaces,
        points,
        &boundary_curves,
        tolerance,
    ) {
        return Ok(result);
    }

    // For two-fillet corners, try polar NURBS as intermediate fallback.
    if n_fillets == 2
        && let Some(result) =
            build_polar_nurbs_blend(boundary, fillet_surfaces, brep, boundary_curves.clone())
    {
        return Ok(result);
    }

    // Final fallback: mitred plane.
    build_mitred_blend(points, boundary_curves)
}

/// Try to build a spherical cap blend.
///
/// Derives the sphere center from the torus geometry of the fillet surfaces:
/// for each `FilletEnd` segment, computes the inward direction from the vertex
/// toward the torus center and averages the directions to find the bisector.
/// The sphere center is placed at `vertex_position + avg_radius * bisector`.
/// Returns `None` if no torus fillets are found or the geometry is degenerate.
fn try_spherical_blend(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    points: &[Point3],
    boundary_curves: &[NurbsCurve],
    tolerance: &Tolerance,
) -> Option<BlendSurfaceResult> {
    if points.len() < 3 {
        return None;
    }

    let mut inward_sum = Vec3::new(0.0, 0.0, 0.0);
    let mut radius_sum = 0.0;
    let mut count = 0usize;

    for seg in &boundary.segments {
        if let BoundarySegment::FilletEnd { fillet_index, .. } = seg
            && let Some(fillet) = fillet_surfaces.get(*fillet_index)
            && let SurfaceGeometry::Torus {
                center,
                minor_radius,
                ..
            } = &fillet.surface
        {
            let dir = Vec3::new(
                center.x - boundary.vertex_position.x,
                center.y - boundary.vertex_position.y,
                center.z - boundary.vertex_position.z,
            );
            if let Some(n) = dir.normalized() {
                inward_sum = Vec3::new(inward_sum.x + n.x, inward_sum.y + n.y, inward_sum.z + n.z);
            }
            radius_sum += minor_radius;
            count += 1;
        }
    }

    if count == 0 {
        return None;
    }

    let avg_radius = radius_sum / count as f64;
    let bisector = inward_sum.normalized()?;

    let center = Point3::new(
        boundary.vertex_position.x + avg_radius * bisector.x,
        boundary.vertex_position.y + avg_radius * bisector.y,
        boundary.vertex_position.z + avg_radius * bisector.z,
    );

    let max_deviation: f64 = points
        .iter()
        .map(|p| (p.distance_to(&center) - avg_radius).abs())
        .fold(0.0_f64, f64::max);

    if max_deviation > tolerance.linear * 10.0 {
        return None;
    }

    let surface = SurfaceGeometry::Sphere {
        center,
        radius: avg_radius,
    };

    Some(BlendSurfaceResult {
        surface,
        boundary_curves: boundary_curves.to_vec(),
    })
}

/// Build a polar NURBS blend surface.
///
/// Uses a polar parameterization centered at the vertex position.
/// The boundary points form the outer ring, and the vertex position
/// is the degenerate pole. Interior control points are adjusted for
/// G1 continuity with adjacent surfaces.
fn build_polar_nurbs_blend(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
    boundary_curves: Vec<NurbsCurve>,
) -> Option<BlendSurfaceResult> {
    use super::boundary::sample_boundary_by_arc_length;

    let n_samples = 24;
    let samples = sample_boundary_by_arc_length(boundary, fillet_surfaces, brep, n_samples);

    if samples.len() < 3 {
        return None;
    }

    let v_pos = boundary.vertex_position;
    let m = samples.len();

    // Build the num_u × num_v polar control point grid.
    // num_u = m samples around the boundary (u direction).
    // num_v = 3 layers in the radial direction (v direction):
    //   v=0: degenerate pole (all at vertex position)
    //   v=0.5: G1-adjusted interior ring
    //   v=1: boundary sample positions
    //
    // NurbsSurface stores row-major: cp(i_u, i_v) = control_points[i_u * num_v + i_v]
    let num_u = m;
    let num_v = 3usize;
    let degree_u = 3usize.min(m - 1);
    let degree_v = 2usize;

    let mut control_points = Vec::with_capacity(num_u * num_v);
    for sample in &samples {
        // v=0: pole
        control_points.push(v_pos);
        // v=1: G1 interior — step inward from the boundary along the surface normal
        let radial_dist = sample.position.distance_to(&v_pos);
        let d = radial_dist / 3.0;
        control_points.push(Point3::new(
            sample.position.x - d * sample.normal.x,
            sample.position.y - d * sample.normal.y,
            sample.position.z - d * sample.normal.z,
        ));
        // v=2: boundary
        control_points.push(sample.position);
    }

    // Clamped uniform knot vector in u: num_u + degree_u + 1 elements,
    // with degree_u+1 repeated values at each end and uniform interior knots.
    let n_interior_u = num_u - degree_u - 1;
    let mut knots_u: Vec<f64> = Vec::with_capacity(num_u + degree_u + 1);
    knots_u.extend(std::iter::repeat_n(0.0, degree_u + 1));
    for i in 1..=n_interior_u {
        knots_u.push(i as f64 / (n_interior_u + 1) as f64);
    }
    knots_u.extend(std::iter::repeat_n(1.0, degree_u + 1));

    // Clamped knot vector in v (degree 2, 3 control points): [0,0,0,1,1,1]
    let knots_v = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let weights = vec![1.0; num_u * num_v];

    let nurbs = popkern_core::geometry::NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        num_u,
        num_v,
        weights,
    };

    Some(BlendSurfaceResult {
        surface: SurfaceGeometry::Nurbs(nurbs),
        boundary_curves,
    })
}

/// Try to build a NURBS patch blend surface.
///
/// Creates a simple bilinear or Coons-style patch through the boundary
/// points. Returns `None` if the construction fails.
#[allow(dead_code)]
fn try_nurbs_blend(
    points: &[Point3],
    boundary_curves: &[NurbsCurve],
) -> Option<BlendSurfaceResult> {
    if points.len() < 4 {
        return None;
    }

    // Build a simple bilinear patch from the first 4 boundary points.
    let n_u = 2;
    let n_v = 2;
    let control_points = vec![points[0], points[1], points[2], points[3]];
    let weights = vec![1.0; 4];

    let nurbs = NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        num_u: n_u,
        num_v: n_v,
        weights,
    };

    Some(BlendSurfaceResult {
        surface: SurfaceGeometry::Nurbs(nurbs),
        boundary_curves: boundary_curves.to_vec(),
    })
}

/// Collect boundary curves from blend boundary segments.
fn collect_boundary_curves(boundary: &BlendBoundary) -> Vec<NurbsCurve> {
    boundary
        .segments
        .iter()
        .map(|seg| match seg {
            BoundarySegment::FilletEnd { curve, .. } => curve.clone(),
            BoundarySegment::FaceEdge { curve, .. } => curve.clone(),
        })
        .collect()
}

/// Sample points from a set of boundary curves.
fn sample_boundary_points(curves: &[NurbsCurve], samples_per_curve: usize) -> Vec<Point3> {
    let mut points = Vec::with_capacity(curves.len() * samples_per_curve);

    for curve in curves {
        let (t_min, t_max) = curve.parameter_range();
        for i in 0..samples_per_curve {
            let t = t_min + (t_max - t_min) * (i as f64) / ((samples_per_curve - 1).max(1) as f64);
            points.push(curve.point_at(t));
        }
    }

    points
}

/// Compute the centroid of a point cloud.
fn compute_centroid(points: &[Point3]) -> Point3 {
    let n = points.len() as f64;
    let x = points.iter().map(|p| p.x).sum::<f64>() / n;
    let y = points.iter().map(|p| p.y).sum::<f64>() / n;
    let z = points.iter().map(|p| p.z).sum::<f64>() / n;
    Point3::new(x, y, z)
}

/// Fit a least-squares plane through a set of points.
///
/// Returns `(centroid, normal)`. The normal is the eigenvector corresponding
/// to the smallest eigenvalue of the covariance matrix — computed via a
/// simplified power-iteration approach for the 3x3 case.
fn fit_plane_least_squares(points: &[Point3]) -> Result<(Point3, Vec3), String> {
    if points.len() < 3 {
        return Err("need at least 3 points to fit a plane".into());
    }

    let centroid = compute_centroid(points);

    // Build 3x3 covariance matrix.
    let mut cov = [[0.0f64; 3]; 3];
    for p in points {
        let dx = p.x - centroid.x;
        let dy = p.y - centroid.y;
        let dz = p.z - centroid.z;
        let d = [dx, dy, dz];
        for r in 0..3 {
            for c in 0..3 {
                cov[r][c] += d[r] * d[c];
            }
        }
    }

    // Find the eigenvector with the smallest eigenvalue using inverse
    // power iteration. For robustness, we use the cross-product of the
    // two largest-variance directions instead.
    let normal = smallest_eigenvector_3x3(&cov)?;

    Ok((centroid, normal))
}

/// Find the eigenvector corresponding to the smallest eigenvalue of a
/// symmetric 3x3 matrix.
///
/// Uses a simplified approach: compute the normal as the cross product
/// of the two most spread directions.
fn smallest_eigenvector_3x3(cov: &[[f64; 3]; 3]) -> Result<Vec3, String> {
    // The diagonal elements give variance along each axis.
    // The smallest-variance axis is approximately the normal.
    // For a more robust result, try cross products of sample vectors.

    // Power iteration to find the largest eigenvector.
    let mut v = [1.0, 0.0, 0.0];
    for _ in 0..50 {
        let mut w = [0.0f64; 3];
        for r in 0..3 {
            for c in 0..3 {
                w[r] += cov[r][c] * v[c];
            }
        }
        let len = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
        if len < 1e-30 {
            break;
        }
        v = [w[0] / len, w[1] / len, w[2] / len];
    }
    let e1 = Vec3::new(v[0], v[1], v[2]);

    // Deflate and find second largest.
    let lambda1: f64 = {
        let mut w = [0.0f64; 3];
        for r in 0..3 {
            for c in 0..3 {
                w[r] += cov[r][c] * v[c];
            }
        }
        (w[0] * v[0] + w[1] * v[1] + w[2] * v[2]).max(0.0)
    };

    let mut cov2 = *cov;
    for r in 0..3 {
        for c in 0..3 {
            cov2[r][c] -= lambda1 * v[r] * v[c];
        }
    }

    let mut v2 = [0.0, 1.0, 0.0];
    for _ in 0..50 {
        let mut w = [0.0f64; 3];
        for r in 0..3 {
            for c in 0..3 {
                w[r] += cov2[r][c] * v2[c];
            }
        }
        let len = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
        if len < 1e-30 {
            break;
        }
        v2 = [w[0] / len, w[1] / len, w[2] / len];
    }
    let e2 = Vec3::new(v2[0], v2[1], v2[2]);

    // Normal = cross product of the two largest eigenvectors.
    let normal = e1.cross(&e2);
    normal
        .normalized()
        .ok_or_else(|| "degenerate covariance — cannot fit plane".to_string())
}

/// Pick a vector perpendicular to the given vector.
fn pick_perpendicular(v: &Vec3) -> Vec3 {
    let candidate = if v.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let perp = v.cross(&candidate);
    perp.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0))
}

#[cfg(test)]
#[path = "surface_tests.rs"]
mod tests;
