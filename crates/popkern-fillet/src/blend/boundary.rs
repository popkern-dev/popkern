//! Blend region boundary identification.
//!
//! At a vertex where multiple fillets meet, this module computes the closed
//! boundary loop around the vertex. The boundary alternates between fillet-end
//! arcs (the transverse edge of each fillet surface near the vertex) and
//! original face edge segments between adjacent fillets.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::fillet_surface::FilletSurfaceResult;

/// Boundary of a vertex blend region.
#[derive(Debug)]
pub struct BlendBoundary {
    /// The vertex being blended.
    pub vertex: VertexHandle,
    /// Position of the vertex (saved before vertex removal by blend insertion).
    pub vertex_position: Point3,
    /// Boundary curve segments (alternating fillet-end arcs and face edge segments).
    pub segments: Vec<BoundarySegment>,
}

/// A single segment of a blend boundary.
#[derive(Debug, Clone)]
pub enum BoundarySegment {
    /// End arc of a fillet surface at the vertex.
    FilletEnd {
        /// Index of the fillet in the global fillet list.
        fillet_index: usize,
        /// The boundary curve segment.
        curve: NurbsCurve,
    },
    /// Edge segment on an original face between two adjacent fillets.
    FaceEdge {
        /// The face this segment lies on.
        face: FaceHandle,
        /// The boundary curve segment.
        curve: NurbsCurve,
    },
}

/// Compute the blend boundary at a vertex where multiple fillets meet.
///
/// Collects the end arcs of each fillet surface near the vertex, plus the
/// original face edge segments between fillets. These form a closed loop
/// around the vertex that defines the blend region.
///
/// # Errors
///
/// Returns an error if any fillet index is out of range or if the vertex
/// position cannot be determined.
pub fn compute_blend_boundary(
    vertex: VertexHandle,
    fillet_indices: &[usize],
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
    tolerance: &Tolerance,
) -> Result<BlendBoundary, String> {
    if fillet_indices.is_empty() {
        return Err("no fillet indices provided for blend boundary".into());
    }

    // Resolve the vertex position.
    let vertex_data = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} not found in BRep"))?;
    let vertex_pos = vertex_data.position;

    let mut segments = Vec::with_capacity(fillet_indices.len() * 2);

    for &fi in fillet_indices {
        let fillet = fillet_surfaces
            .get(fi)
            .ok_or_else(|| format!("fillet index {fi} out of range"))?;

        // Extract the fillet end arc near the vertex.
        let end_curve = extract_fillet_end_near_vertex(fillet, &vertex_pos, tolerance)?;

        segments.push(BoundarySegment::FilletEnd {
            fillet_index: fi,
            curve: end_curve,
        });
    }

    // For a complete boundary we would also insert face-edge segments between
    // consecutive fillet ends. This requires walking the topology around the
    // vertex. For now, the fillet-end segments alone define the blend region;
    // face-edge segments can be added when topology walking is available.

    Ok(BlendBoundary {
        vertex,
        vertex_position: vertex_pos,
        segments,
    })
}

/// Extract the end curve of a fillet surface closest to a vertex.
///
/// Takes the trim curves of the fillet and finds the endpoint nearest to the
/// vertex position, then builds a short transverse arc between the two trim
/// curve endpoints.
fn extract_fillet_end_near_vertex(
    fillet: &FilletSurfaceResult,
    vertex_pos: &Point3,
    _tolerance: &Tolerance,
) -> Result<NurbsCurve, String> {
    // Find which end of trim_curve_1 is nearest to the vertex.
    let (t1_min, t1_max) = fillet.trim_curve_1.parameter_range();
    let p1_start = fillet.trim_curve_1.point_at(t1_min);
    let p1_end = fillet.trim_curve_1.point_at(t1_max);

    let (t2_min, t2_max) = fillet.trim_curve_2.parameter_range();
    let p2_start = fillet.trim_curve_2.point_at(t2_min);
    let p2_end = fillet.trim_curve_2.point_at(t2_max);

    // Pick the pair of endpoints nearest to the vertex.
    let d_start = p1_start.distance_to(vertex_pos) + p2_start.distance_to(vertex_pos);
    let d_end = p1_end.distance_to(vertex_pos) + p2_end.distance_to(vertex_pos);

    let (arc_start, arc_end) = if d_start < d_end {
        (p1_start, p2_start)
    } else {
        (p1_end, p2_end)
    };

    // Build a simple linear NURBS between the two endpoints.
    // A true arc could be constructed using the fillet surface geometry.
    let mid = Point3::new(
        (arc_start.x + arc_end.x) * 0.5,
        (arc_start.y + arc_end.y) * 0.5,
        (arc_start.z + arc_end.z) * 0.5,
    );

    crate::center_curve::fit_nurbs_through_points(&[arc_start, mid, arc_end])
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::geometry::{Point3, Vec3};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn test_tolerance() -> Tolerance {
        Tolerance {
            linear: 1e-3,
            angular: 1e-3,
            parametric: 1e-6,
        }
    }

    fn make_fillet_result(
        start1: Point3,
        end1: Point3,
        start2: Point3,
        end2: Point3,
    ) -> FilletSurfaceResult {
        let trim_1 = crate::center_curve::fit_nurbs_through_points(&[start1, end1])
            .expect("fit should succeed");
        let trim_2 = crate::center_curve::fit_nurbs_through_points(&[start2, end2])
            .expect("fit should succeed");

        FilletSurfaceResult {
            surface: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
            trim_curve_1: trim_1,
            trim_curve_2: trim_2,
        }
    }

    #[test]
    fn blend_boundary_two_fillets() {
        let tol = test_tolerance();
        let brep = BRep::unit_cube(tol).expect("cube should build");

        // Pick the first vertex.
        let vertex = brep.vertices.keys().next().expect("cube has vertices");
        let vertex_pos = brep.vertices.get(vertex).expect("vertex exists").position;

        // Create two mock fillet results near the vertex.
        let f1 = make_fillet_result(
            vertex_pos,
            Point3::new(vertex_pos.x + 1.0, vertex_pos.y, vertex_pos.z),
            Point3::new(vertex_pos.x, vertex_pos.y + 0.1, vertex_pos.z),
            Point3::new(vertex_pos.x + 1.0, vertex_pos.y + 0.1, vertex_pos.z),
        );
        let f2 = make_fillet_result(
            vertex_pos,
            Point3::new(vertex_pos.x, vertex_pos.y + 1.0, vertex_pos.z),
            Point3::new(vertex_pos.x + 0.1, vertex_pos.y, vertex_pos.z),
            Point3::new(vertex_pos.x + 0.1, vertex_pos.y + 1.0, vertex_pos.z),
        );

        let fillets = vec![f1, f2];
        let result = compute_blend_boundary(vertex, &[0, 1], &fillets, &brep, &tol);
        assert!(result.is_ok(), "blend boundary failed: {:?}", result.err());

        let boundary = result.expect("already checked");
        assert_eq!(
            boundary.segments.len(),
            2,
            "should have 2 fillet-end segments"
        );

        // Verify the boundary is around the right vertex.
        assert_eq!(boundary.vertex, vertex);

        // Each segment should be a FilletEnd.
        for seg in &boundary.segments {
            assert!(
                matches!(seg, BoundarySegment::FilletEnd { .. }),
                "expected FilletEnd segment"
            );
        }
    }

    #[test]
    fn blend_boundary_empty_fillets() {
        let tol = test_tolerance();
        let brep = BRep::unit_cube(tol).expect("cube should build");
        let vertex = brep.vertices.keys().next().expect("cube has vertices");

        let result = compute_blend_boundary(vertex, &[], &[], &brep, &tol);
        assert!(result.is_err(), "should fail with no fillet indices");
    }

    #[test]
    fn blend_boundary_invalid_fillet_index() {
        let tol = test_tolerance();
        let brep = BRep::unit_cube(tol).expect("cube should build");
        let vertex = brep.vertices.keys().next().expect("cube has vertices");

        let result = compute_blend_boundary(vertex, &[99], &[], &brep, &tol);
        assert!(result.is_err(), "should fail with out-of-range index");
    }
}
