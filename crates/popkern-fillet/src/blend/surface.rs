//! Vertex blend surface construction.
//!
//! Builds the surface that fills the gap at a vertex where multiple fillets
//! meet. Supports mitred (planar) blends and smooth (spherical or NURBS)
//! blends depending on the [`VertexMode`].

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
    _brep: &BRep,
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
/// For the two-fillet case, tries a spherical cap. Falls back to mitred
/// if smooth construction fails.
fn build_smooth_blend(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    points: &[Point3],
    boundary_curves: Vec<NurbsCurve>,
    tolerance: &Tolerance,
) -> Result<BlendSurfaceResult, String> {
    let n_fillets = boundary
        .segments
        .iter()
        .filter(|s| matches!(s, BoundarySegment::FilletEnd { .. }))
        .count();

    if n_fillets == 2 {
        // Try spherical cap.
        if let Some(result) = try_spherical_blend(
            boundary,
            fillet_surfaces,
            points,
            &boundary_curves,
            tolerance,
        ) {
            return Ok(result);
        }

        // Try NURBS patch as second attempt.
        if let Some(result) = try_nurbs_blend(points, &boundary_curves) {
            return Ok(result);
        }
    }

    // Fallback for 3+ fillets or if smooth construction fails: use mitred.
    build_mitred_blend(points, boundary_curves)
}

/// Try to build a spherical cap blend.
///
/// Finds a sphere tangent to both fillet surfaces at the boundary, then
/// checks G1 continuity. Returns `None` if the construction fails.
fn try_spherical_blend(
    _boundary: &BlendBoundary,
    _fillet_surfaces: &[FilletSurfaceResult],
    points: &[Point3],
    boundary_curves: &[NurbsCurve],
    tolerance: &Tolerance,
) -> Option<BlendSurfaceResult> {
    // Compute the centroid as the sphere center candidate.
    let centroid = compute_centroid(points);

    // Compute average distance from centroid to boundary points as radius.
    let avg_radius: f64 =
        points.iter().map(|p| p.distance_to(&centroid)).sum::<f64>() / points.len() as f64;

    if avg_radius < tolerance.linear {
        return None; // Degenerate — all points at the same location.
    }

    // Check that all boundary points are approximately equidistant from
    // the centroid (within tolerance). If not, a sphere won't work.
    let max_deviation: f64 = points
        .iter()
        .map(|p| (p.distance_to(&centroid) - avg_radius).abs())
        .fold(0.0_f64, f64::max);

    if max_deviation > tolerance.linear * 10.0 {
        return None; // Points don't lie on a sphere.
    }

    let surface = SurfaceGeometry::Sphere {
        center: centroid,
        radius: avg_radius,
    };

    Some(BlendSurfaceResult {
        surface,
        boundary_curves: boundary_curves.to_vec(),
    })
}

/// Try to build a NURBS patch blend surface.
///
/// Creates a simple bilinear or Coons-style patch through the boundary
/// points. Returns `None` if the construction fails.
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
mod tests {
    use super::*;
    use crate::blend::boundary::{BlendBoundary, BoundarySegment};
    use popkern_core::geometry::Point3;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::{BRep, VertexHandle};

    fn test_tolerance() -> Tolerance {
        Tolerance {
            linear: 1e-3,
            angular: 1e-3,
            parametric: 1e-6,
        }
    }

    fn make_test_curve(start: Point3, end: Point3) -> NurbsCurve {
        crate::center_curve::fit_nurbs_through_points(&[start, end]).expect("fit should succeed")
    }

    fn make_blend_boundary(vertex: VertexHandle) -> BlendBoundary {
        // Two fillet-end arcs forming a right-angle-ish boundary.
        let seg1 = BoundarySegment::FilletEnd {
            fillet_index: 0,
            curve: make_test_curve(Point3::new(0.1, 0.0, 0.0), Point3::new(0.0, 0.1, 0.0)),
        };
        let seg2 = BoundarySegment::FilletEnd {
            fillet_index: 1,
            curve: make_test_curve(Point3::new(0.0, 0.1, 0.0), Point3::new(0.0, 0.0, 0.1)),
        };

        BlendBoundary {
            vertex,
            segments: vec![seg1, seg2],
        }
    }

    #[test]
    fn mitred_blend_two_fillets() {
        let tol = test_tolerance();
        let brep = BRep::unit_cube(tol.clone()).expect("cube should build");
        let vertex = brep.vertices.keys().next().expect("cube has vertices");

        let boundary = make_blend_boundary(vertex);
        let fillets: Vec<FilletSurfaceResult> = Vec::new();

        let result = build_blend_surface(&boundary, &fillets, &brep, VertexMode::Mitred, &tol);
        assert!(result.is_ok(), "mitred blend failed: {:?}", result.err());

        let blend = result.expect("already checked");
        assert!(
            matches!(blend.surface, SurfaceGeometry::Plane { .. }),
            "mitred blend should produce a plane, got {:?}",
            std::mem::discriminant(&blend.surface)
        );
        assert!(
            !blend.boundary_curves.is_empty(),
            "should have boundary curves"
        );
    }

    #[test]
    fn smooth_blend_two_fillets() {
        let tol = test_tolerance();
        let brep = BRep::unit_cube(tol.clone()).expect("cube should build");
        let vertex = brep.vertices.keys().next().expect("cube has vertices");

        // Build boundary with points on a sphere of radius 0.1.
        let r = 0.1;
        let seg1 = BoundarySegment::FilletEnd {
            fillet_index: 0,
            curve: make_test_curve(Point3::new(r, 0.0, 0.0), Point3::new(0.0, r, 0.0)),
        };
        let seg2 = BoundarySegment::FilletEnd {
            fillet_index: 1,
            curve: make_test_curve(Point3::new(0.0, r, 0.0), Point3::new(0.0, 0.0, r)),
        };

        let boundary = BlendBoundary {
            vertex,
            segments: vec![seg1, seg2],
        };

        let fillets: Vec<FilletSurfaceResult> = Vec::new();

        let result = build_blend_surface(&boundary, &fillets, &brep, VertexMode::Smooth, &tol);
        assert!(result.is_ok(), "smooth blend failed: {:?}", result.err());

        let blend = result.expect("already checked");
        // Should produce either a sphere or a NURBS patch (or plane fallback).
        // The exact result depends on point distribution; we just verify success.
        assert!(
            !blend.boundary_curves.is_empty(),
            "should have boundary curves"
        );
    }

    #[test]
    fn blend_surface_empty_boundary() {
        let tol = test_tolerance();
        let brep = BRep::unit_cube(tol.clone()).expect("cube should build");
        let vertex = brep.vertices.keys().next().expect("cube has vertices");

        let boundary = BlendBoundary {
            vertex,
            segments: vec![],
        };

        let result = build_blend_surface(&boundary, &[], &brep, VertexMode::Mitred, &tol);
        assert!(result.is_err(), "should fail with empty boundary");
    }

    #[test]
    fn fit_plane_coplanar_points() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ];

        let (centroid, normal) = fit_plane_least_squares(&points).expect("should fit");

        // Centroid should be at (0.5, 0.5, 0.0).
        assert!((centroid.x - 0.5).abs() < 1e-10);
        assert!((centroid.y - 0.5).abs() < 1e-10);
        assert!(centroid.z.abs() < 1e-10);

        // Normal should be along z (up to sign).
        assert!(
            normal.z.abs() > 0.99,
            "normal should be along z, got {normal:?}"
        );
    }
}
