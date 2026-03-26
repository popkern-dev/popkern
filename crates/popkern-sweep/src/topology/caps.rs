//! Cap face creation for swept and lofted solids.
//!
//! Creates planar cap faces at the start and end of a sweep/loft path
//! by wiring boundary half-edges into a loop.

use std::f64::consts::FRAC_1_SQRT_2;

use popkern_core::constructors::builder::BRepBuilder;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{FaceHandle, HalfEdgeHandle, ShellHandle};

use crate::error::SweepError;

/// Create a planar cap face from boundary half-edges.
///
/// The `half_edges` slice contains the twin half-edges of the side face
/// boundary edges that form the cap boundary. These are wired into a loop
/// and a planar surface is fitted through the boundary points.
///
/// `is_start_cap` controls the orientation of the cap normal: start caps
/// face opposite to the path direction, end caps face along it.
pub fn create_cap_face(
    builder: &mut BRepBuilder,
    shell: ShellHandle,
    half_edges: &[HalfEdgeHandle],
    cap_curves: &[NurbsCurve],
    is_start_cap: bool,
) -> Result<FaceHandle, SweepError> {
    if half_edges.is_empty() {
        return Err(SweepError::Topology("cap has no boundary edges".into()));
    }

    let face = builder.add_face(shell);

    // Compute a planar surface from the cap boundary.
    let cap_surface = compute_cap_surface(builder, half_edges, is_start_cap);
    builder.set_face_surface(face, cap_surface);

    builder
        .add_loop(half_edges, face)
        .map_err(|e| SweepError::Topology(format!("cap loop: {e}")))?;

    // Attach cap curve geometry to each boundary edge.
    for (he_handle, curve) in half_edges.iter().zip(cap_curves.iter()) {
        let edge_handle = builder.brep().half_edges[*he_handle].edge;
        let (t_start, t_end) = curve.parameter_range();
        builder.set_edge_curve(
            edge_handle,
            CurveOnEdge {
                geometry: CurveGeometry::Nurbs(curve.clone()),
                t_start,
                t_end,
            },
        );
    }

    // Attach pcurves to cap half-edges.
    attach_cap_pcurves(builder, face, half_edges, cap_curves);

    Ok(face)
}

/// Attach pcurves to cap face half-edges.
///
/// For multi-edge caps, projects the edge curve endpoints into the
/// face's UV space and creates degree-1 line pcurves.
///
/// For self-loop caps (1 half-edge with a full-circle NURBS curve),
/// creates a degree-2 rational NURBS circle pcurve matching the
/// pattern used in the cylinder constructor.
fn attach_cap_pcurves(
    builder: &mut BRepBuilder,
    face: FaceHandle,
    half_edges: &[HalfEdgeHandle],
    cap_curves: &[NurbsCurve],
) {
    let surface = match &builder.brep().faces[face].surface {
        Some(s) => s.geometry.clone(),
        None => return,
    };

    if half_edges.len() == 1 {
        // Self-loop cap: full circle edge. Use a degree-2 rational circle pcurve.
        let he = half_edges[0];
        let edge = builder.brep().half_edges[he].edge;
        if let Some(ref curve_on_edge) = builder.brep().edges[edge].curve {
            let start_3d = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
            if let Some(center_uv) = compute_cap_uv_center(builder, half_edges, &surface)
                && let Some((u_start, v_start)) = project_point_to_uv(&surface, &start_3d)
            {
                let du = u_start - center_uv.0;
                let dv = v_start - center_uv.1;
                let radius = (du * du + dv * dv).sqrt();
                if radius > 1e-15 {
                    let pcurve = cap_circle_pcurve(radius, center_uv.0, center_uv.1);
                    builder.set_half_edge_pcurve(he, pcurve);
                }
            }
        }
    } else {
        // Multi-edge cap: project endpoints to UV and create line pcurves.
        for (&he, curve) in half_edges.iter().zip(cap_curves.iter()) {
            let (t_start, t_end) = curve.parameter_range();
            let p_start = curve.point_at(t_start);
            let p_end = curve.point_at(t_end);

            // Determine direction: compare half-edge origin with curve start.
            let origin_vh = builder.brep().half_edges[he].origin;
            let origin_pos = builder.brep().vertices[origin_vh].position;
            let forward = origin_pos.distance_to(&p_start) < origin_pos.distance_to(&p_end);

            let (pt_a, pt_b) = if forward {
                (p_start, p_end)
            } else {
                (p_end, p_start)
            };

            if let (Some((u0, v0)), Some((u1, v1))) = (
                project_point_to_uv(&surface, &pt_a),
                project_point_to_uv(&surface, &pt_b),
            ) {
                let pcurve = NurbsCurve {
                    degree: 1,
                    knots: vec![0.0, 0.0, 1.0, 1.0],
                    control_points: vec![Point3::new(u0, v0, 0.0), Point3::new(u1, v1, 0.0)],
                    weights: vec![1.0, 1.0],
                };
                builder.set_half_edge_pcurve(he, pcurve);
            }
        }
    }
}

/// Compute the UV center of a self-loop cap face.
///
/// Samples points along the edge curve at several parameter values,
/// projects them into UV space, and returns their centroid. This avoids
/// the degenerate case where the only vertex (seam point) projects to
/// the UV origin.
fn compute_cap_uv_center(
    builder: &BRepBuilder,
    half_edges: &[HalfEdgeHandle],
    surface: &SurfaceGeometry,
) -> Option<(f64, f64)> {
    let brep = builder.brep();

    // Sample 3D points from the edge curve to find the true UV center.
    let he = half_edges[0];
    let edge = brep.half_edges[he].edge;
    let curve_on_edge = brep.edges[edge].curve.as_ref()?;

    let t0 = curve_on_edge.t_start;
    let t1 = curve_on_edge.t_end;
    let n_samples = 8;

    let mut u_sum = 0.0;
    let mut v_sum = 0.0;
    let mut count = 0usize;

    for i in 0..n_samples {
        let t = t0 + (t1 - t0) * (i as f64 / n_samples as f64);
        let pt = curve_on_edge.geometry.point_at(t);
        if let Some((u, v)) = project_point_to_uv(surface, &pt) {
            u_sum += u;
            v_sum += v;
            count += 1;
        }
    }

    if count == 0 {
        return None;
    }
    Some((u_sum / count as f64, v_sum / count as f64))
}

/// Build a degree-2 rational NURBS circle pcurve in the cap's UV plane.
///
/// Creates a full-circle pcurve centered at `(cu, cv)` with the given
/// `radius`, starting at `(cu + radius, cv)`.
fn cap_circle_pcurve(radius: f64, cu: f64, cv: f64) -> NurbsCurve {
    let r = radius;
    let w = FRAC_1_SQRT_2;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            Point3::new(cu + r, cv, 0.0),
            Point3::new(cu + r, cv + r, 0.0),
            Point3::new(cu, cv + r, 0.0),
            Point3::new(cu - r, cv + r, 0.0),
            Point3::new(cu - r, cv, 0.0),
            Point3::new(cu - r, cv - r, 0.0),
            Point3::new(cu, cv - r, 0.0),
            Point3::new(cu + r, cv - r, 0.0),
            Point3::new(cu + r, cv, 0.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Compute a planar surface for a cap face.
///
/// Estimates the plane from the half-edge origin vertices using
/// a centroid and Newell's method for the normal.
fn compute_cap_surface(
    builder: &BRepBuilder,
    half_edges: &[HalfEdgeHandle],
    is_start_cap: bool,
) -> SurfaceOnFace {
    let brep = builder.brep();

    // Collect vertex positions from half-edge origins.
    let points: Vec<Point3> = half_edges
        .iter()
        .map(|&he| brep.vertices[brep.half_edges[he].origin].position)
        .collect();

    let centroid = compute_centroid(&points);
    let normal = compute_newell_normal(&points);

    // Choose u/v axes perpendicular to the normal.
    let (u_axis, v_axis) = orthonormal_basis(&normal);

    // Flip normal for start cap (faces away from path start).
    let oriented_normal = if is_start_cap { normal * -1.0 } else { normal };

    SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: centroid,
            normal: oriented_normal,
            u_axis,
            v_axis,
        },
        orientation: true,
    }
}

/// Compute the centroid of a set of points.
fn compute_centroid(points: &[Point3]) -> Point3 {
    let n = points.len().max(1) as f64;
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sz = 0.0;
    for p in points {
        sx += p.x;
        sy += p.y;
        sz += p.z;
    }
    Point3::new(sx / n, sy / n, sz / n)
}

/// Compute a normal using Newell's method.
fn compute_newell_normal(points: &[Point3]) -> Vec3 {
    let n = points.len();
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;
    for i in 0..n {
        let curr = &points[i];
        let next = &points[(i + 1) % n];
        nx += (curr.y - next.y) * (curr.z + next.z);
        ny += (curr.z - next.z) * (curr.x + next.x);
        nz += (curr.x - next.x) * (curr.y + next.y);
    }
    let len = (nx * nx + ny * ny + nz * nz).sqrt();
    if len < 1e-30 {
        // Degenerate — fall back to Z-up.
        Vec3::new(0.0, 0.0, 1.0)
    } else {
        Vec3::new(nx / len, ny / len, nz / len)
    }
}

/// Build an orthonormal basis perpendicular to a given normal.
fn orthonormal_basis(normal: &Vec3) -> (Vec3, Vec3) {
    let candidate = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = normal
        .cross(&candidate)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = normal
        .cross(&u)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn centroid_of_triangle() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(0.0, 3.0, 0.0),
        ];
        let c = compute_centroid(&points);
        assert!((c.x - 1.0).abs() < 1e-10);
        assert!((c.y - 1.0).abs() < 1e-10);
        assert!(c.z.abs() < 1e-10);
    }

    #[test]
    fn newell_normal_of_xy_triangle() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let n = compute_newell_normal(&points);
        // Should point in +Z direction.
        assert!(n.z.abs() > 0.9, "expected Z-normal, got {:?}", n);
    }

    #[test]
    fn orthonormal_basis_perpendicular() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let (u, v) = orthonormal_basis(&normal);
        // u and v should be perpendicular to normal and to each other.
        assert!(normal.dot(&u).abs() < 1e-10);
        assert!(normal.dot(&v).abs() < 1e-10);
        assert!(u.dot(&v).abs() < 1e-10);
    }

    #[test]
    fn degenerate_newell_falls_back() {
        // All points are the same.
        let points = vec![Point3::new(1.0, 1.0, 1.0), Point3::new(1.0, 1.0, 1.0)];
        let n = compute_newell_normal(&points);
        // Should fall back to Z-up.
        assert!((n.z - 1.0).abs() < 1e-10);
    }
}
