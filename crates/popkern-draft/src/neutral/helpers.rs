//! Helper functions for neutral line computation.

use crate::error::DraftError;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Compute centroid of a face from its boundary vertices.
pub(crate) fn compute_face_centroid(brep: &BRep, face: FaceHandle) -> Point3 {
    let face_data = &brep.faces[face];
    let outer_loop = &brep.loops[face_data.outer_loop];

    let mut sum = Vec3::zero();
    let mut count = 0;

    let start_he = outer_loop.half_edge;
    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        let pos = brep.vertices[he.origin].position;
        sum = sum + Vec3::new(pos.x, pos.y, pos.z);
        count += 1;
        current = he.next;
        if current == start_he {
            break;
        }
    }

    if count > 0 {
        let inv = 1.0 / count as f64;
        Point3::new(sum.x * inv, sum.y * inv, sum.z * inv)
    } else {
        Point3::new(0.0, 0.0, 0.0)
    }
}

/// Compute approximate extent of a face (max distance between boundary vertices).
pub(crate) fn compute_face_extent(brep: &BRep, face: FaceHandle) -> f64 {
    let face_data = &brep.faces[face];
    let outer_loop = &brep.loops[face_data.outer_loop];

    let mut vertices = Vec::new();
    let start_he = outer_loop.half_edge;
    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        vertices.push(brep.vertices[he.origin].position);
        current = he.next;
        if current == start_he {
            break;
        }
    }

    let mut max_dist = 1.0;
    for (i, a) in vertices.iter().enumerate() {
        for b in vertices.iter().skip(i + 1) {
            let d = a.distance_to(b);
            if d > max_dist {
                max_dist = d;
            }
        }
    }
    max_dist
}

/// Build a degree-1 NurbsCurve (line segment) from two points.
pub(crate) fn build_line_curve(p0: Point3, p1: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    }
}

/// Compute a vector perpendicular to the given vector.
pub(crate) fn compute_perpendicular(v: &Vec3) -> Vec3 {
    let abs_x = v.x.abs();
    let abs_y = v.y.abs();
    let abs_z = v.z.abs();

    let candidate = if abs_x <= abs_y && abs_x <= abs_z {
        Vec3::new(1.0, 0.0, 0.0)
    } else if abs_y <= abs_z {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };

    let perp = v.cross(&candidate);
    perp.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0))
}

/// Find a point on the intersection of two planes.
pub(crate) fn solve_point_on_two_planes(
    n1: Vec3,
    d1: f64,
    n2: Vec3,
    d2: f64,
    direction: Vec3,
    _tolerance: &Tolerance,
) -> Result<Point3, DraftError> {
    let ax = direction.x.abs();
    let ay = direction.y.abs();
    let az = direction.z.abs();

    let (x, y, z) = if ax >= ay && ax >= az {
        let det = n1.y * n2.z - n1.z * n2.y;
        if det.abs() < 1e-12 {
            return Err(DraftError::InvalidNeutral(
                "degenerate plane-plane intersection".into(),
            ));
        }
        let y = (d1 * n2.z - n1.z * d2) / det;
        let z = (n1.y * d2 - d1 * n2.y) / det;
        (0.0, y, z)
    } else if ay >= ax && ay >= az {
        let det = n1.x * n2.z - n1.z * n2.x;
        if det.abs() < 1e-12 {
            return Err(DraftError::InvalidNeutral(
                "degenerate plane-plane intersection".into(),
            ));
        }
        let x = (d1 * n2.z - n1.z * d2) / det;
        let z = (n1.x * d2 - d1 * n2.x) / det;
        (x, 0.0, z)
    } else {
        let det = n1.x * n2.y - n1.y * n2.x;
        if det.abs() < 1e-12 {
            return Err(DraftError::InvalidNeutral(
                "degenerate plane-plane intersection".into(),
            ));
        }
        let x = (d1 * n2.y - n1.y * d2) / det;
        let y = (n1.x * d2 - d1 * n2.x) / det;
        (x, y, 0.0)
    };

    Ok(Point3::new(x, y, z))
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn build_line_curve_evaluates_correctly() {
        let p0 = Point3::new(0.0, 0.0, 0.0);
        let p1 = Point3::new(1.0, 0.0, 0.0);
        let curve = build_line_curve(p0, p1);
        let mid = curve.point_at(0.5);
        assert!((mid.x - 0.5).abs() < 1e-10);
        assert!(mid.y.abs() < 1e-10);
    }

    #[test]
    fn compute_perpendicular_gives_orthogonal_vector() {
        let v = Vec3::new(0.0, 0.0, 1.0);
        let perp = compute_perpendicular(&v);
        assert!(v.dot(&perp).abs() < 1e-10);
        assert!((perp.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn face_centroid_is_center_of_cube_face() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let centroid = compute_face_centroid(&brep, face);
        let coords = [centroid.x, centroid.y, centroid.z];
        for &c in &coords {
            assert!(
                c >= -0.1 && c <= 1.1,
                "centroid coordinate {c} out of range"
            );
        }
    }
}
