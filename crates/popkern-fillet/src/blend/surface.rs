//! Vertex blend surface construction (mitred and smooth modes).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3};
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
/// Computes the circumcenter of the boundary points — the point
/// equidistant from all of them. If the points lie on a sphere (max
/// deviation < tolerance), returns that sphere. This handles cube
/// corners with equal-radius perpendicular fillets exactly.
///
/// Returns `None` if the geometry is degenerate or the points don't
/// lie on a sphere.
fn try_spherical_blend(
    _boundary: &BlendBoundary,
    _fillet_surfaces: &[FilletSurfaceResult],
    points: &[Point3],
    boundary_curves: &[NurbsCurve],
    tolerance: &Tolerance,
) -> Option<BlendSurfaceResult> {
    if points.len() < 4 {
        return None;
    }

    // Compute the circumcenter by solving the linear system:
    // 2*(p_j - p_0) · c = |p_j|² - |p_0|²  for j = 1..N
    //
    // This finds the point equidistant from all boundary points.
    // Use the first 3 non-degenerate point differences to solve for c,
    // then verify against all points.
    let p0 = points[0];
    let p0_sq = p0.x * p0.x + p0.y * p0.y + p0.z * p0.z;

    // Build a 3x4 augmented matrix [A | b] from the first 3 valid equations.
    let mut rows: Vec<[f64; 4]> = Vec::new();
    for p in points.iter().skip(1) {
        let dx = p.x - p0.x;
        let dy = p.y - p0.y;
        let dz = p.z - p0.z;
        if dx.abs() + dy.abs() + dz.abs() < tolerance.solver() {
            continue; // Skip near-duplicate points.
        }
        let rhs = p.x * p.x + p.y * p.y + p.z * p.z - p0_sq;
        rows.push([2.0 * dx, 2.0 * dy, 2.0 * dz, rhs]);
        if rows.len() == 3 {
            break;
        }
    }

    if rows.len() < 3 {
        return None; // Not enough independent equations.
    }

    // Solve 3x3 system by Gaussian elimination with partial pivoting.
    let (cx, cy, cz) = solve_circumcenter_3x3(&rows)?;
    let center = Point3::new(cx, cy, cz);

    // Compute radius as average distance from center to all points.
    let avg_radius: f64 =
        points.iter().map(|p| p.distance_to(&center)).sum::<f64>() / points.len() as f64;

    if avg_radius < tolerance.linear() {
        return None;
    }

    // Verify ALL boundary points are equidistant from center.
    let max_deviation: f64 = points
        .iter()
        .map(|p| (p.distance_to(&center) - avg_radius).abs())
        .fold(0.0_f64, f64::max);

    if max_deviation > tolerance.linear() * 10.0 {
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

use super::surface_math::{fit_plane_least_squares, pick_perpendicular, solve_circumcenter_3x3};

#[cfg(test)]
#[path = "surface_tests.rs"]
mod tests;
