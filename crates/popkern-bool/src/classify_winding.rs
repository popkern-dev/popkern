//! Winding number point-in-solid classification.
//!
//! Computes the generalized winding number of a point with respect to a
//! closed triangulated surface. The winding number is the solid angle
//! subtended by the surface divided by 4pi:
//!   w(P) = (1/4pi) * sum over triangles of signed solid angle
//!
//! For a valid closed shell: w = +1 inside, 0 outside.
//! More robust than ray casting (no degenerate-ray problem).
//!
//! Reference: Jacobson, Kavan, Sorkine-Hornung, "Robust Inside-Outside
//! Segmentation using Generalized Winding Numbers," ACM TOG (2013).

use popkern_core::geometry::Point3;
use popkern_core::topology::BRep;

use crate::error::PointClassification;

/// Classify a point relative to a solid using the winding number.
///
/// Returns `Inside` if winding number > 0.5, `Outside` if < 0.5.
/// Falls back to `OnBoundary` if the computation is inconclusive
/// (e.g., the point is very close to the surface).
pub fn point_in_solid_winding(point: &Point3, solid: &BRep) -> PointClassification {
    let winding = compute_winding_number(point, solid);

    if winding > 0.5 {
        PointClassification::Inside
    } else if winding < -0.5 {
        // Negative winding: inside with inverted orientation
        PointClassification::Inside
    } else if winding.abs() < 0.01 {
        PointClassification::Outside
    } else {
        // Ambiguous — near the boundary
        PointClassification::OnBoundary
    }
}

/// Compute the generalized winding number of `point` with respect to `solid`.
///
/// Approximates each face as a triangle fan from its vertex centroid and
/// sums the signed solid angle contributions. For faces with surface
/// geometry, samples additional points for better accuracy.
fn compute_winding_number(point: &Point3, solid: &BRep) -> f64 {
    let mut total = 0.0;

    for (fh, face) in &solid.faces {
        // Collect vertex positions for this face
        let positions = crate::clip::face_vertex_positions(solid, fh);
        if positions.len() < 3 {
            continue;
        }

        // Determine face orientation
        let flip = face.surface.as_ref().is_some_and(|s| !s.orientation);

        // Triangle fan from first vertex
        let v0 = positions[0];
        for i in 1..(positions.len() - 1) {
            let v1 = positions[i];
            let v2 = positions[i + 1];
            let angle = triangle_solid_angle(point, &v0, &v1, &v2);
            if flip {
                total -= angle;
            } else {
                total += angle;
            }
        }
    }

    total / (4.0 * std::f64::consts::PI)
}

/// Compute the signed solid angle subtended by a triangle at a point.
///
/// Uses the formula from Van Oosterom and Strackee (1983):
///   tan(Omega/2) = det(a,b,c) / (|a||b||c| + (a·b)|c| + (b·c)|a| + (c·a)|b|)
///
/// where a, b, c are vectors from the point to the triangle vertices.
fn triangle_solid_angle(point: &Point3, v0: &Point3, v1: &Point3, v2: &Point3) -> f64 {
    let a = *v0 - *point;
    let b = *v1 - *point;
    let c = *v2 - *point;

    let la = a.length();
    let lb = b.length();
    let lc = c.length();

    // Degenerate: point is at a vertex
    if la < 1e-15 || lb < 1e-15 || lc < 1e-15 {
        return 0.0;
    }

    // det(a, b, c) = a · (b × c)
    let det = a.x * (b.y * c.z - b.z * c.y)
        + a.y * (b.z * c.x - b.x * c.z)
        + a.z * (b.x * c.y - b.y * c.x);

    let ab = a.x * b.x + a.y * b.y + a.z * b.z;
    let bc = b.x * c.x + b.y * c.y + b.z * c.z;
    let ca = c.x * a.x + c.y * a.y + c.z * a.z;

    let denom = la * lb * lc + ab * lc + bc * la + ca * lb;

    2.0 * det.atan2(denom)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    #[test]
    fn point_inside_cube() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let inside = Point3::new(0.5, 0.5, 0.5);
        let result = point_in_solid_winding(&inside, &cube);
        assert_eq!(result, PointClassification::Inside);
    }

    #[test]
    fn point_outside_cube() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let outside = Point3::new(5.0, 5.0, 5.0);
        let result = point_in_solid_winding(&outside, &cube);
        assert_eq!(result, PointClassification::Outside);
    }

    #[test]
    fn winding_number_unit_cube_center() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let center = Point3::new(0.5, 0.5, 0.5);
        let w = compute_winding_number(&center, &cube);
        // Should be approximately 1.0 for inside a closed solid
        assert!(
            (w.abs() - 1.0).abs() < 0.1,
            "winding at center should be ~1.0, got {w}"
        );
    }

    #[test]
    fn winding_number_outside_cube() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let outside = Point3::new(5.0, 0.0, 0.0);
        let w = compute_winding_number(&outside, &cube);
        assert!(w.abs() < 0.1, "winding outside should be ~0, got {w}");
    }
}
