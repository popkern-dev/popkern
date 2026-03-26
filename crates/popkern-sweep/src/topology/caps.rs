//! Cap face creation for swept and lofted solids.
//!
//! Creates planar cap faces at the start and end of a sweep/loft path
//! by wiring boundary half-edges into a loop.

use popkern_core::constructors::builder::BRepBuilder;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
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

    Ok(face)
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
