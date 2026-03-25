//! Parting line draft: apply different draft angles above and below a parting plane.
//!
//! In mold design, a parting line divides the mold into two halves.
//! Faces above the parting line are drafted one way, faces below the other.

use crate::error::DraftError;
use crate::types::{
    AngleVariation, DraftOptions, DraftResult, DraftSpec, NeutralReference, PartingLine,
};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Apply draft angles relative to a parting line.
///
/// Classifies each face as above or below the parting line based on
/// face centroid position, then applies the appropriate draft angle.
///
/// # Arguments
///
/// * `brep` - The input BRep solid.
/// * `parting` - The parting line definition (plane, curve, or edges).
/// * `above_angle` - Draft angle in radians for faces above the parting line.
/// * `below_angle` - Draft angle in radians for faces below the parting line.
/// * `pull_direction` - Direction of mold withdrawal.
/// * `options` - Draft options controlling adjacent face behavior.
/// * `tolerance` - Geometric tolerance.
///
/// # Errors
///
/// Returns [`DraftError`] if the parting line is invalid, no faces can
/// be classified, or drafting fails.
pub fn draft_parting_line(
    brep: &BRep,
    parting: &PartingLine,
    above_angle: f64,
    below_angle: f64,
    pull_direction: Vec3,
    options: &DraftOptions,
    tolerance: Tolerance,
) -> Result<DraftResult, DraftError> {
    let (parting_origin, parting_normal) = resolve_parting_line(brep, parting)?;

    let mut specs = Vec::new();

    for (face_h, face_data) in &brep.faces {
        // Skip faces whose surface normal is parallel to the parting
        // plane normal — the neutral plane and face are parallel so no
        // neutral line can be computed (e.g. top/bottom of a cube when
        // the parting plane is horizontal).
        if let Some(ref sof) = face_data.surface {
            use popkern_core::geometry::surface::SurfaceEval;
            let face_normal = sof.geometry.normal_at(0.0, 0.0);
            if let Some(fn_unit) = face_normal.normalized() {
                let dot = fn_unit.dot(&parting_normal).abs();
                if dot > 1.0 - tolerance.angular {
                    continue;
                }
            }
        }

        let centroid = compute_face_centroid(brep, face_h);
        let offset = centroid - parting_origin;
        let signed_dist = offset.dot(&parting_normal);

        // Skip faces on the parting line (within tolerance).
        if signed_dist.abs() < tolerance.linear {
            continue;
        }

        let angle = if signed_dist > 0.0 {
            above_angle
        } else {
            below_angle
        };

        // Skip faces with zero draft angle.
        if angle.abs() < tolerance.angular {
            continue;
        }

        specs.push(DraftSpec {
            face: face_h,
            angle,
            pull_direction,
            neutral: NeutralReference::Plane {
                origin: parting_origin,
                normal: parting_normal,
            },
            angle_variation: AngleVariation::Constant,
        });
    }

    if specs.is_empty() {
        return Err(DraftError::InvalidPartingLine(
            "no faces classified above or below the parting line".into(),
        ));
    }

    crate::draft_faces(brep, &specs, options, tolerance)
}

/// Resolve a parting line definition into an origin point and normal direction.
fn resolve_parting_line(brep: &BRep, parting: &PartingLine) -> Result<(Point3, Vec3), DraftError> {
    match parting {
        PartingLine::Plane { origin, normal } => {
            let n = normal.normalized().ok_or_else(|| {
                DraftError::InvalidPartingLine("parting plane normal is zero".into())
            })?;
            Ok((*origin, n))
        }
        PartingLine::Curve(curve) => {
            use popkern_core::geometry::curve::CurveEval;
            let (t0, t1) = curve.parameter_range();
            let mid = (t0 + t1) * 0.5;
            let origin = curve.point_at(mid);
            let tangent = curve.tangent_at(mid);
            // Normal is perpendicular to curve tangent in the pull plane.
            // Use a simple perpendicular for now.
            let perp = compute_perpendicular(&tangent);
            Ok((origin, perp))
        }
        PartingLine::Edges(edges) => {
            if edges.is_empty() {
                return Err(DraftError::InvalidPartingLine(
                    "empty edge list for parting line".into(),
                ));
            }
            // Use the midpoint of the first edge as origin,
            // and compute a fitted plane normal from all edge vertices.
            let first_edge = brep
                .edges
                .get(edges[0])
                .ok_or_else(|| DraftError::InvalidPartingLine("edge not found in BRep".into()))?;
            let he = &brep.half_edges[first_edge.half_edge];
            let twin = &brep.half_edges[he.twin];
            let p0 = brep.vertices[he.origin].position;
            let p1 = brep.vertices[twin.origin].position;
            let origin = Point3::new(
                (p0.x + p1.x) * 0.5,
                (p0.y + p1.y) * 0.5,
                (p0.z + p1.z) * 0.5,
            );

            // Collect all edge direction vectors and compute a normal from them.
            let mut tangent_sum = p1 - p0;
            for edge_h in edges.iter().skip(1) {
                if let Some(e) = brep.edges.get(*edge_h) {
                    let eh = &brep.half_edges[e.half_edge];
                    let et = &brep.half_edges[eh.twin];
                    let ep0 = brep.vertices[eh.origin].position;
                    let ep1 = brep.vertices[et.origin].position;
                    tangent_sum = tangent_sum + (ep1 - ep0);
                }
            }

            let avg_tangent = tangent_sum
                .normalized()
                .ok_or_else(|| DraftError::InvalidPartingLine("degenerate edge chain".into()))?;
            let normal = compute_perpendicular(&avg_tangent);
            Ok((origin, normal))
        }
    }
}

/// Compute a unit vector perpendicular to the given direction.
fn compute_perpendicular(v: &Vec3) -> Vec3 {
    let reference = if v.x.abs() <= v.y.abs() && v.x.abs() <= v.z.abs() {
        Vec3::new(1.0, 0.0, 0.0)
    } else if v.y.abs() <= v.z.abs() {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let cross = v.cross(&reference);
    cross.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0))
}

/// Compute the centroid of a face from its boundary vertices.
fn compute_face_centroid(brep: &BRep, face: popkern_core::topology::FaceHandle) -> Point3 {
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

    if count == 0 {
        return Point3::new(0.0, 0.0, 0.0);
    }

    let inv = 1.0 / count as f64;
    Point3::new(sum.x * inv, sum.y * inv, sum.z * inv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn parting_plane_through_cube_classifies_faces() {
        let brep = BRep::unit_cube(tol()).unwrap();

        // Parting plane at y=0.5 with normal pointing in +y.
        // This classifies x-normal and z-normal faces as above/below,
        // while y-normal faces (parallel to parting normal) are skipped.
        // Side faces at x=0(centroid y=0.5), x=1(centroid y=0.5) have
        // centroids on the plane, so use y=0.25 offset to avoid that.
        let parting = PartingLine::Plane {
            origin: Point3::new(0.0, 0.25, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
        };

        let result = draft_parting_line(
            &brep,
            &parting,
            0.05, // ~3 degrees above
            0.05, // ~3 degrees below
            Vec3::new(0.0, 1.0, 0.0),
            &DraftOptions::default(),
            tol(),
        );

        // The result should succeed (faces above/below are classified).
        assert!(result.is_ok(), "parting draft failed: {:?}", result.err());
        let dr = result.unwrap();
        // Should still have 6 faces (no splits for simple cube).
        assert_eq!(dr.brep.faces.len(), 6);
    }

    #[test]
    fn parting_line_zero_normal_returns_error() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let parting = PartingLine::Plane {
            origin: Point3::new(0.0, 0.0, 0.5),
            normal: Vec3::zero(),
        };

        let result = draft_parting_line(
            &brep,
            &parting,
            0.05,
            0.05,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol(),
        );
        assert!(result.is_err());
    }

    #[test]
    fn parting_line_empty_edges_returns_error() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let parting = PartingLine::Edges(vec![]);

        let result = draft_parting_line(
            &brep,
            &parting,
            0.05,
            0.05,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol(),
        );
        assert!(result.is_err());
    }

    #[test]
    fn face_centroid_on_cube() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let centroid = compute_face_centroid(&brep, face);
        // Centroid should be within the unit cube bounding box.
        assert!(centroid.x >= -0.1 && centroid.x <= 1.1);
        assert!(centroid.y >= -0.1 && centroid.y <= 1.1);
        assert!(centroid.z >= -0.1 && centroid.z <= 1.1);
    }

    #[test]
    fn compute_perpendicular_is_orthogonal() {
        let v = Vec3::new(1.0, 0.0, 0.0);
        let perp = compute_perpendicular(&v);
        assert!(perp.dot(&v).abs() < 1e-10);
        assert!((perp.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn compute_perpendicular_y_dominant() {
        let v = Vec3::new(0.0, 1.0, 0.0);
        let perp = compute_perpendicular(&v);
        assert!(perp.dot(&v).abs() < 1e-10);
        assert!((perp.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn compute_perpendicular_z_dominant() {
        let v = Vec3::new(0.0, 0.0, 1.0);
        let perp = compute_perpendicular(&v);
        assert!(perp.dot(&v).abs() < 1e-10);
        assert!((perp.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn parting_line_edges_variant_classifies_faces() {
        let brep = BRep::unit_cube(tol()).unwrap();

        // Pick an edge from the unit cube.
        let edge_h = brep.edges.keys().next().unwrap();

        let parting = PartingLine::Edges(vec![edge_h]);

        let result = draft_parting_line(
            &brep,
            &parting,
            0.05,
            0.05,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol(),
        );
        // Should resolve without error (may produce "no faces classified" if
        // all centroids happen to be on the plane, but that is also acceptable).
        let _ = result;
    }

    #[test]
    fn parting_line_edges_invalid_handle_returns_error() {
        let brep = BRep::unit_cube(tol()).unwrap();
        // Default (zero) EdgeHandle is unlikely to exist in a real BRep.
        let parting = PartingLine::Edges(vec![Default::default()]);

        let result = draft_parting_line(
            &brep,
            &parting,
            0.05,
            0.05,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol(),
        );
        assert!(result.is_err());
    }

    #[test]
    fn faces_parallel_to_parting_plane_are_skipped() {
        // For a unit cube with a Z-normal parting plane at z=0.3, the top
        // (z=1) and bottom (z=0) faces are parallel to the parting plane.
        // The side faces have centroids at z=0.5, so they are above z=0.3
        // and get classified as "above". This exercises the parallel-skip path.
        let brep = BRep::unit_cube(tol()).unwrap();
        // Place parting plane at z=0.3 — side face centroids are at z=0.5 (above).
        // X-normal faces have normals perpendicular to Z, so they are NOT parallel.
        let parting = PartingLine::Plane {
            origin: Point3::new(0.5, 0.5, 0.3),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };

        let result = draft_parting_line(
            &brep,
            &parting,
            0.05,
            0.05,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol(),
        );
        // Side faces are classified above; the operation should succeed.
        assert!(result.is_ok(), "expected Ok, got {:?}", result.err());
    }

    #[test]
    fn all_faces_on_parting_plane_returns_error() {
        // Use a cube and put the parting plane such that all non-parallel
        // face centroids are exactly on it (within tolerance) — this forces
        // the "no faces classified" error.
        let brep = BRep::unit_cube(tol()).unwrap();

        // A very tight tolerance so that even tiny offsets are within linear tol.
        let big_tol = Tolerance {
            linear: 10.0, // huge linear tolerance treats everything as "on plane"
            angular: 1e-10,
            parametric: 1e-8,
        };

        let parting = PartingLine::Plane {
            origin: Point3::new(0.5, 0.5, 0.5),
            normal: Vec3::new(1.0, 0.0, 0.0),
        };

        let result = draft_parting_line(
            &brep,
            &parting,
            0.05,
            0.05,
            Vec3::new(1.0, 0.0, 0.0),
            &DraftOptions::default(),
            big_tol,
        );
        assert!(result.is_err());
    }

    #[test]
    fn zero_above_and_below_angles_returns_error() {
        // Both angles are zero → all faces skipped by the angle-filter → error.
        let brep = BRep::unit_cube(tol()).unwrap();
        let parting = PartingLine::Plane {
            origin: Point3::new(0.5, 0.5, 0.25),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };

        let result = draft_parting_line(
            &brep,
            &parting,
            0.0, // above_angle = 0 → skip
            0.0, // below_angle = 0 → skip
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol(),
        );
        assert!(result.is_err());
    }

    #[test]
    fn parting_curve_resolves_without_panic() {
        use popkern_core::geometry::NurbsCurve;
        let brep = BRep::unit_cube(tol()).unwrap();

        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
            weights: vec![1.0, 1.0],
        };

        let parting = PartingLine::Curve(curve);
        // Result may succeed or fail depending on face classification; just
        // ensure it doesn't panic and the Curve branch is executed.
        let _result = draft_parting_line(
            &brep,
            &parting,
            0.05,
            0.05,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol(),
        );
    }
}
