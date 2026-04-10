//! Chamfer curve and surface computation.
//!
//! For each edge to be chamfered, computes offset curves on each adjacent face
//! and builds a ruled surface (flat plane) between them.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::analysis::EdgeAnalysis;
use crate::center_curve::fit_nurbs_through_points;
use crate::types::ChamferKind;

/// Result of chamfer geometry computation for one edge.
#[derive(Debug)]
pub struct ChamferResult {
    /// The ruled surface connecting the two trim curves.
    pub surface: SurfaceGeometry,
    /// 3D trim curve on the first face.
    pub trim_curve_1: NurbsCurve,
    /// 3D trim curve on the second face.
    pub trim_curve_2: NurbsCurve,
    /// Parameter-space trim curve on the first face (Point3(u, v, 0) convention).
    pub pcurve_1: NurbsCurve,
    /// Parameter-space trim curve on the second face (Point3(u, v, 0) convention).
    pub pcurve_2: NurbsCurve,
}

/// Compute chamfer geometry for an edge.
///
/// Based on the `ChamferKind`, computes the offset distance on each face,
/// finds trim curves at the offset distance, builds a ruled surface between
/// them, and computes parameter-space curves on each face.
pub fn compute_chamfer(
    analysis: &EdgeAnalysis,
    brep: &BRep,
    kind: &ChamferKind,
    tolerance: &Tolerance,
) -> Result<ChamferResult, String> {
    // 1. Determine offset distances on each face.
    let (dist_1, dist_2) = compute_offset_distances(analysis, kind)?;

    // 2. Get the edge curve.
    let edge = brep
        .edges
        .get(analysis.edge)
        .ok_or_else(|| format!("edge {:?} not found", analysis.edge))?;
    let curve = edge
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {:?} has no curve geometry", analysis.edge))?;

    // 3. Sample points along the edge.
    let n_samples = 20usize;
    let (t_start, t_end) = (curve.t_start, curve.t_end);

    let mut edge_points = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let t = t_start + (t_end - t_start) * (i as f64) / ((n_samples - 1).max(1) as f64);
        edge_points.push(curve.geometry.point_at(t));
    }

    // 4. Compute offset points on each face.
    let offset_points_1 =
        compute_offset_points(&edge_points, &analysis.surface_1, dist_1, tolerance)?;
    let offset_points_2 =
        compute_offset_points(&edge_points, &analysis.surface_2, dist_2, tolerance)?;

    // 5. Fit trim curves through the offset points.
    let trim_curve_1 = fit_nurbs_through_points(&offset_points_1)?;
    let trim_curve_2 = fit_nurbs_through_points(&offset_points_2)?;

    // 6. Compute pcurves by projecting onto face UV spaces.
    let pcurve_1 = compute_pcurve(&offset_points_1, &analysis.surface_1)?;
    let pcurve_2 = compute_pcurve(&offset_points_2, &analysis.surface_2)?;

    // 7. Build the ruled surface between the two trim curves.
    let surface = build_ruled_surface(&offset_points_1, &offset_points_2)?;

    Ok(ChamferResult {
        surface,
        trim_curve_1,
        trim_curve_2,
        pcurve_1,
        pcurve_2,
    })
}

/// Compute offset distances on each face based on chamfer kind.
fn compute_offset_distances(
    analysis: &EdgeAnalysis,
    kind: &ChamferKind,
) -> Result<(f64, f64), String> {
    match kind {
        ChamferKind::Equal(d) => Ok((*d, *d)),
        ChamferKind::TwoDistance {
            face_1: ref_face,
            dist_1,
            dist_2,
        } => {
            if *ref_face == analysis.face_1 {
                Ok((*dist_1, *dist_2))
            } else if *ref_face == analysis.face_2 {
                Ok((*dist_2, *dist_1))
            } else {
                Err("reference face does not match either adjacent face".into())
            }
        }
        ChamferKind::DistanceAngle {
            face: ref_face,
            distance,
            angle,
        } => {
            let other_dist = distance * angle.tan();
            if *ref_face == analysis.face_1 {
                Ok((*distance, other_dist))
            } else if *ref_face == analysis.face_2 {
                Ok((other_dist, *distance))
            } else {
                Err("reference face does not match either adjacent face".into())
            }
        }
    }
}

/// Compute offset points on a surface at a given distance from edge points.
///
/// For each edge point, computes the direction on the surface perpendicular
/// to the edge and offsets along that direction by the given distance.
fn compute_offset_points(
    edge_points: &[Point3],
    surface: &SurfaceGeometry,
    distance: f64,
    tolerance: &Tolerance,
) -> Result<Vec<Point3>, String> {
    if edge_points.len() < 2 {
        return Err("need at least 2 edge points".into());
    }

    let mut offset_points = Vec::with_capacity(edge_points.len());

    for (i, pt) in edge_points.iter().enumerate() {
        // Compute approximate edge tangent from neighboring points.
        let tangent = if i == 0 {
            edge_points[1] - edge_points[0]
        } else if i == edge_points.len() - 1 {
            edge_points[i] - edge_points[i - 1]
        } else {
            edge_points[i + 1] - edge_points[i - 1]
        };
        let tangent_len = tangent.length();
        if tangent_len < tolerance.numerical_floor() {
            offset_points.push(*pt);
            continue;
        }
        let tangent_unit = Vec3::new(
            tangent.x / tangent_len,
            tangent.y / tangent_len,
            tangent.z / tangent_len,
        );

        // Find the surface normal at this point.
        let (u, v) = find_uv_for_point(surface, pt);
        let normal = surface.normal_at(u, v);
        let normal_len = normal.length();
        if normal_len < tolerance.numerical_floor() {
            offset_points.push(*pt);
            continue;
        }

        // Offset direction: perpendicular to the edge tangent, lying on the surface.
        // offset_dir = normalize(normal x tangent)
        let offset_raw = normal.cross(&tangent_unit);
        let offset_len = offset_raw.length();
        if offset_len < tolerance.numerical_floor() {
            offset_points.push(*pt);
            continue;
        }
        let offset_dir = Vec3::new(
            offset_raw.x / offset_len,
            offset_raw.y / offset_len,
            offset_raw.z / offset_len,
        );

        let offset_pt = *pt + offset_dir * distance;

        // Project offset point back onto the surface for accuracy.
        let projected = project_onto_surface(surface, &offset_pt);
        offset_points.push(projected);
    }

    Ok(offset_points)
}

/// Project a 3D point onto a surface (closest point).
fn project_onto_surface(surface: &SurfaceGeometry, point: &Point3) -> Point3 {
    let (u, v) = find_uv_for_point(surface, point);
    surface.point_at(u, v)
}

/// Find (u, v) parameters on a surface closest to a 3D point.
fn find_uv_for_point(surface: &SurfaceGeometry, point: &Point3) -> (f64, f64) {
    if let Some(uv) = popkern_core::geometry::uv_projection::project_point_to_uv(surface, point) {
        return uv;
    }

    // Fallback: grid search.
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let u_lo = if u_min.is_finite() { u_min } else { -10.0 };
    let u_hi = if u_max.is_finite() { u_max } else { 10.0 };
    let v_lo = if v_min.is_finite() { v_min } else { -10.0 };
    let v_hi = if v_max.is_finite() { v_max } else { 10.0 };

    let n = 20;
    let mut best_u = (u_lo + u_hi) * 0.5;
    let mut best_v = (v_lo + v_hi) * 0.5;
    let mut best_dist = f64::INFINITY;

    for iu in 0..n {
        let u = u_lo + (u_hi - u_lo) * (iu as f64) / ((n - 1).max(1) as f64);
        for iv in 0..n {
            let v = v_lo + (v_hi - v_lo) * (iv as f64) / ((n - 1).max(1) as f64);
            let sp = surface.point_at(u, v);
            let d = sp.distance_to(point);
            if d < best_dist {
                best_dist = d;
                best_u = u;
                best_v = v;
            }
        }
    }

    (best_u, best_v)
}

/// Compute a parameter-space (pcurve) by projecting 3D points into surface UV.
fn compute_pcurve(points: &[Point3], surface: &SurfaceGeometry) -> Result<NurbsCurve, String> {
    let mut uv_points = Vec::with_capacity(points.len());
    for pt in points {
        let (u, v) = find_uv_for_point(surface, pt);
        uv_points.push(Point3::new(u, v, 0.0));
    }
    fit_nurbs_through_points(&uv_points)
}

/// Build a ruled (degree-1 in v) NURBS surface between two sets of points.
///
/// The surface is degree-1 in both u (along edge) and v (across chamfer).
/// Row 0 = points on face 1, row 1 = points on face 2.
fn build_ruled_surface(
    points_1: &[Point3],
    points_2: &[Point3],
) -> Result<SurfaceGeometry, String> {
    let n = points_1.len();
    if n < 2 || points_2.len() != n {
        return Err("ruled surface requires matching point arrays with >= 2 points".into());
    }

    // Control points: row-major, [row0 = points_1, row1 = points_2].
    // U direction = along the edge (n points), V direction = across chamfer (2 points).
    let mut control_points = Vec::with_capacity(2 * n);
    control_points.extend_from_slice(points_1);
    control_points.extend_from_slice(points_2);

    let degree_u = 1.min(n - 1);
    let degree_v = 1;

    let knots_u = build_clamped_uniform_knots(n, degree_u);
    let knots_v = vec![0.0, 0.0, 1.0, 1.0]; // degree-1 with 2 control points

    let weights = vec![1.0; 2 * n];

    let surface = NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        num_u: n,
        num_v: 2,
        weights,
    };

    Ok(SurfaceGeometry::Nurbs(surface))
}

/// Build clamped uniform knot vector for `n` control points at given degree.
fn build_clamped_uniform_knots(n: usize, degree: usize) -> Vec<f64> {
    let n_knots = n + degree + 1;
    let n_interior = n.saturating_sub(degree);
    let mut knots = Vec::with_capacity(n_knots);

    knots.extend(std::iter::repeat_n(0.0, degree + 1));
    for i in 1..n_interior {
        knots.push(i as f64 / n_interior as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));

    knots
}

#[cfg(test)]
#[path = "chamfer_curve_tests.rs"]
mod tests;
