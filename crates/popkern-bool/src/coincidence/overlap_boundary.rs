//! Overlap boundary computation for coincident face pairs.
//!
//! When two faces lie on the same plane (CoincidentRegion from
//! `intersect_surfaces`), compute the boundary of their overlap region
//! as intersection curves that can be used by the graph-merge pipeline.
//!
//! For coplanar planar faces, this is a 2D polygon intersection problem
//! in the shared UV parameter space.

use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip::face_vertex_positions;

/// Compute the overlap boundary curves between two coplanar faces.
///
/// Returns intersection curves along the boundary of the overlap region.
/// These curves can be fed into the graph-merge pipeline the same way
/// that surface-surface intersection curves are.
///
/// Currently handles planar-planar coplanar faces only.
pub fn compute_overlap_boundary(
    a: &BRep,
    face_a: FaceHandle,
    b: &BRep,
    face_b: FaceHandle,
    tolerance: &Tolerance,
) -> Vec<CurveOnSurface> {
    let surf_a = match a.faces.get(face_a).and_then(|f| f.surface.as_ref()) {
        Some(s) => &s.geometry,
        None => return Vec::new(),
    };
    let surf_b = match b.faces.get(face_b).and_then(|f| f.surface.as_ref()) {
        Some(s) => &s.geometry,
        None => return Vec::new(),
    };

    // Only handle planar-planar for now
    let (origin, u_axis, v_axis) = match surf_a {
        SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        } => (*origin, *u_axis, *v_axis),
        _ => return Vec::new(),
    };

    if !matches!(surf_b, SurfaceGeometry::Plane { .. }) {
        return Vec::new();
    }

    // Project both face boundaries into the shared UV space
    let poly_a = project_face_to_uv(a, face_a, &origin, &u_axis, &v_axis);
    let poly_b = project_face_to_uv(b, face_b, &origin, &u_axis, &v_axis);

    if poly_a.len() < 3 || poly_b.len() < 3 {
        return Vec::new();
    }

    // Find edges of polygon B that cross inside polygon A, and vice versa.
    // Each crossing edge produces a boundary curve segment.
    let mut curves = Vec::new();

    // Check edges of A against boundary of B
    for i in 0..poly_a.len() {
        let j = (i + 1) % poly_a.len();
        let p0 = poly_a[i];
        let p1 = poly_a[j];

        // If this edge of A lies inside B, it's part of the overlap boundary
        let mid_u = (p0.0 + p1.0) * 0.5;
        let mid_v = (p0.1 + p1.1) * 0.5;
        if point_in_polygon_2d(mid_u, mid_v, &poly_b) {
            // This edge segment is inside both polygons — it's an overlap boundary
            let start_3d = origin + u_axis * p0.0 + v_axis * p0.1;
            let end_3d = origin + u_axis * p1.0 + v_axis * p1.1;
            if let Some(cos) =
                make_line_curve_on_surface(start_3d, end_3d, surf_a, surf_b, tolerance)
            {
                curves.push(cos);
            }
        }
    }

    // Check edges of B against boundary of A
    for i in 0..poly_b.len() {
        let j = (i + 1) % poly_b.len();
        let p0 = poly_b[i];
        let p1 = poly_b[j];

        let mid_u = (p0.0 + p1.0) * 0.5;
        let mid_v = (p0.1 + p1.1) * 0.5;
        if point_in_polygon_2d(mid_u, mid_v, &poly_a) {
            let start_3d = origin + u_axis * p0.0 + v_axis * p0.1;
            let end_3d = origin + u_axis * p1.0 + v_axis * p1.1;
            if let Some(cos) =
                make_line_curve_on_surface(start_3d, end_3d, surf_a, surf_b, tolerance)
            {
                curves.push(cos);
            }
        }
    }

    curves
}

/// Project a face's vertex positions into a 2D UV space defined by (origin, u_axis, v_axis).
fn project_face_to_uv(
    brep: &BRep,
    face: FaceHandle,
    origin: &Point3,
    u_axis: &Vec3,
    v_axis: &Vec3,
) -> Vec<(f64, f64)> {
    let positions = face_vertex_positions(brep, face);
    positions
        .iter()
        .map(|p| {
            let d = *p - *origin;
            (d.dot(u_axis), d.dot(v_axis))
        })
        .collect()
}

/// 2D point-in-polygon test using ray casting (crossing number).
pub(crate) fn point_in_polygon_2d(u: f64, v: f64, poly: &[(f64, f64)]) -> bool {
    let n = poly.len();
    let mut crossings = 0;
    for i in 0..n {
        let j = (i + 1) % n;
        let (ui, vi) = poly[i];
        let (uj, vj) = poly[j];

        if (vi <= v && vj > v) || (vj <= v && vi > v) {
            let t = (v - vi) / (vj - vi);
            let u_cross = ui + t * (uj - ui);
            if u < u_cross {
                crossings += 1;
            }
        }
    }
    crossings % 2 == 1
}

/// Create a CurveOnSurface for a line segment between two 3D points.
fn make_line_curve_on_surface(
    start: Point3,
    end: Point3,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> Option<CurveOnSurface> {
    use popkern_core::geometry::surface::SurfaceEval;

    let dir = end - start;
    let len = dir.length();
    if len < tolerance.linear {
        return None;
    }

    let curve_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, len, len],
        control_points: vec![start, end],
        weights: vec![1.0, 1.0],
    };

    // Project endpoints onto both surfaces for pcurves
    let (u1_start, v1_start) = surf_a.closest_point(&start)?;
    let (u1_end, v1_end) = surf_a.closest_point(&end)?;
    let pcurve_a = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, len, len],
        control_points: vec![
            Point3::new(u1_start, v1_start, 0.0),
            Point3::new(u1_end, v1_end, 0.0),
        ],
        weights: vec![1.0, 1.0],
    };

    let (u2_start, v2_start) = surf_b.closest_point(&start)?;
    let (u2_end, v2_end) = surf_b.closest_point(&end)?;
    let pcurve_b = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, len, len],
        control_points: vec![
            Point3::new(u2_start, v2_start, 0.0),
            Point3::new(u2_end, v2_end, 0.0),
        ],
        weights: vec![1.0, 1.0],
    };

    Some(CurveOnSurface {
        curve_3d,
        curve_uv1: pcurve_a,
        curve_uv2: pcurve_b,
        is_closed: false,
    })
}

#[cfg(test)]
#[path = "overlap_boundary_tests.rs"]
mod tests;
