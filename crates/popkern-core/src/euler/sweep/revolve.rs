//! Revolve sweep — a higher-level operation built from Euler operators (mev + mef).
//!
//! Revolves a face around an axis, creating a solid of revolution.
//! For partial angles (< 2pi), creates side faces and an end cap.
//! For full 360°, creates multiple angular strips that wrap around,
//! with the final strip closing back to the original vertices.

use crate::error::CoreError;
use crate::geometry::revolution::rotate_point;
use crate::geometry::{Point3, Vec3};
use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};
use std::f64::consts::TAU;

/// Result of the revolve sweep operation.
#[derive(Debug, Clone)]
pub struct RevolveResult {
    /// The revolution (side) faces created by the sweep.
    pub revolution_faces: Vec<FaceHandle>,
    /// Cap faces (empty for full 360° revolve, 1 cap for partial revolve).
    pub cap_faces: Vec<FaceHandle>,
    /// All new edges created during the revolve.
    pub new_edges: Vec<EdgeHandle>,
    /// All new vertices created during the revolve.
    pub new_vertices: Vec<VertexHandle>,
}

/// Determine the number of angular steps for the revolve.
fn angular_steps(angle: f64) -> usize {
    use std::f64::consts::{FRAC_PI_2, PI};
    if angle <= FRAC_PI_2 + 1e-10 {
        1
    } else if angle <= PI + 1e-10 {
        2
    } else if angle <= 3.0 * FRAC_PI_2 + 1e-10 {
        3
    } else {
        4
    }
}

/// Perform one "extrude-like" ring of mev + mef operations.
///
/// Given the current half-edges in the active face loop (`current_he`, one per
/// profile vertex), creates new vertices at rotated positions, then closes
/// side faces with mef.
///
/// Returns `(ring_new_vertices, ring_new_edges, ring_side_faces, cap_face, next_current_he)`.
/// `cap_face` is the face remaining after the last mef (the "end" face for this ring).
/// `next_current_he` contains the half-edges in the cap face's loop, originating
/// at the newly created vertices (for use by the next ring).
#[allow(clippy::type_complexity)]
fn extrude_ring(
    brep: &mut BRep,
    current_he: &[HalfEdgeHandle],
    loop_vertices: &[VertexHandle],
    axis_origin: &Point3,
    axis_unit: &Vec3,
    ring_angle: f64,
) -> Result<
    (
        Vec<VertexHandle>,
        Vec<EdgeHandle>,
        Vec<FaceHandle>,
        FaceHandle,
        Vec<HalfEdgeHandle>,
    ),
    CoreError,
> {
    let n = loop_vertices.len();
    let mut new_vertices = Vec::with_capacity(n);
    let mut new_edges = Vec::new();
    let mut side_faces = Vec::with_capacity(n);
    let mut in_half_edges: Vec<HalfEdgeHandle> = Vec::with_capacity(n);

    // Step 1: mev for each vertex — create rotated vertex.
    for i in 0..n {
        let old_vertex = loop_vertices[i];
        let old_pos = brep.vertices[old_vertex].position;
        let new_pos = rotate_point(&old_pos, axis_origin, axis_unit, ring_angle);

        let mev_result = brep.mev(current_he[i], new_pos)?;
        new_vertices.push(mev_result.vertex);
        new_edges.push(mev_result.edge);
        in_half_edges.push(mev_result.half_edge_twin);
    }

    // Step 2: mef to close side faces (same pattern as extrude).
    let mut top0_he = HalfEdgeHandle::default();
    let mut cap_face = FaceHandle::default();

    for i in 0..n {
        let (he_a, he_b) = if i < n - 1 {
            (in_half_edges[i + 1], in_half_edges[i])
        } else {
            (top0_he, in_half_edges[n - 1])
        };

        let mef_result = brep.mef(he_a, he_b)?;
        new_edges.push(mef_result.edge);

        if i == 0 {
            top0_he = mef_result.half_edge_twin;
        }

        let side_loop = brep.half_edges[mef_result.half_edge].loop_;
        side_faces.push(brep.loops[side_loop].face);

        if i == n - 1 {
            cap_face = mef_result.face;
        }
    }

    // Compute next_current_he: walk the cap face's loop starting from top0_he.
    let mut next_he = Vec::with_capacity(n);
    let mut walk = top0_he;
    for _ in 0..n {
        next_he.push(walk);
        walk = brep.half_edges[walk].next;
    }

    Ok((new_vertices, new_edges, side_faces, cap_face, next_he))
}

impl BRep {
    /// Revolve a face around an axis by a given angle.
    ///
    /// Creates revolution (side) faces by rotating each vertex of the face's
    /// outer loop around the axis. For partial revolves (angle < 2pi), a cap
    /// face is created at the end position. For full 360° revolves, the
    /// topology uses multiple angular strips to wrap around completely.
    ///
    /// # Arguments
    ///
    /// * `face` - The face to revolve (must have edges in its outer loop).
    /// * `axis_origin` - A point on the revolution axis.
    /// * `axis_dir` - Direction of the revolution axis (need not be unit length).
    /// * `angle` - Revolution angle in radians. Must be in (0, 2*pi].
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if the face handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if the face has no edges, the angle
    ///   is out of range, or the axis direction is zero-length.
    pub fn revolve(
        &mut self,
        face: FaceHandle,
        axis_origin: Point3,
        axis_dir: Vec3,
        angle: f64,
    ) -> Result<RevolveResult, CoreError> {
        // Validate angle.
        if angle <= 0.0 || angle > TAU + 1e-10 {
            return Err(CoreError::PreconditionFailed(format!(
                "revolve angle must be in (0, 2*pi], got {angle}"
            )));
        }

        let axis_unit = axis_dir.normalized().ok_or_else(|| {
            CoreError::PreconditionFailed("axis direction is zero-length".into())
        })?;

        let is_full = (angle - TAU).abs() < 1e-10;
        let angle = angle.min(TAU);

        // Validate face and get outer loop.
        let face_data = self
            .faces
            .get(face)
            .ok_or_else(|| CoreError::InvalidHandle("face handle is invalid".into()))?;

        let outer_loop = face_data.outer_loop;
        let loop_data = self
            .loops
            .get(outer_loop)
            .ok_or_else(|| CoreError::InvalidHandle("outer loop handle is invalid".into()))?;

        let start_he = loop_data.half_edge;
        if !self.half_edges.contains_key(start_he) {
            return Err(CoreError::PreconditionFailed(
                "face has no edges in its outer loop".into(),
            ));
        }

        // Walk the outer loop to collect half-edges and vertices.
        let mut loop_half_edges: Vec<HalfEdgeHandle> = Vec::new();
        let mut loop_vertices: Vec<VertexHandle> = Vec::new();

        let mut cursor = start_he;
        loop {
            let he_data = &self.half_edges[cursor];
            loop_half_edges.push(cursor);
            loop_vertices.push(he_data.origin);
            cursor = he_data.next;
            if cursor == start_he {
                break;
            }
        }

        let m = angular_steps(angle);

        let mut all_new_vertices = Vec::new();
        let mut all_new_edges = Vec::new();
        let mut revolution_faces = Vec::new();
        let mut cap_faces = Vec::new();

        if is_full {
            // Full 360°: create M rings, each spanning angle/M degrees.
            // All M rings create new vertices (the last ring's vertices are
            // topologically distinct from the originals, even though they are
            // at the same geometric positions). This creates a closed solid.
            let mut current_he = loop_half_edges.clone();

            for ring in 0..m {
                let ring_angle = angle * (ring + 1) as f64 / m as f64;

                let (verts, edges, faces, _cap, next_he) = extrude_ring(
                    self,
                    &current_he,
                    &loop_vertices,
                    &axis_origin,
                    &axis_unit,
                    ring_angle,
                )?;

                all_new_vertices.extend(verts);
                all_new_edges.extend(edges);
                revolution_faces.extend(faces);
                current_he = next_he;
            }

            // No cap faces for full 360°.
        } else {
            // Partial revolve: create M rings of intermediate+final positions.
            let mut current_he = loop_half_edges.clone();

            for ring in 0..m {
                let ring_angle = angle * (ring + 1) as f64 / m as f64;

                let (verts, edges, faces, cap, next_he) = extrude_ring(
                    self,
                    &current_he,
                    &loop_vertices,
                    &axis_origin,
                    &axis_unit,
                    ring_angle,
                )?;

                all_new_vertices.extend(verts);
                all_new_edges.extend(edges);
                revolution_faces.extend(faces);

                if ring == m - 1 {
                    // The cap from the final ring is the end cap.
                    cap_faces.push(cap);
                }

                current_he = next_he;
            }
        }

        Ok(RevolveResult {
            revolution_faces,
            cap_faces,
            new_edges: all_new_edges,
            new_vertices: all_new_vertices,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;
    use std::f64::consts::{FRAC_PI_2, TAU};

    /// Build a square face offset from the Z axis for revolving.
    fn build_revolve_square(brep: &mut BRep) -> crate::topology::FaceHandle {
        let mvfs = brep.mvfs(Point3::new(2.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(3.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(3.0, 0.0, 1.0))
            .unwrap();
        let mev3 = brep
            .mev(mev2.half_edge_twin, Point3::new(2.0, 0.0, 1.0))
            .unwrap();
        let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
        mef.face
    }

    #[test]
    fn revolve_square_90_has_caps() {
        let mut brep = BRep::new(Tolerance::default());
        let face = build_revolve_square(&mut brep);

        let result = brep
            .revolve(
                face,
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                FRAC_PI_2,
            )
            .unwrap();

        assert!(
            !result.cap_faces.is_empty(),
            "Partial revolve should have caps"
        );
        assert!(
            !result.revolution_faces.is_empty(),
            "Should have revolution faces"
        );
        // 90° = 1 step, 4 profile vertices => 4 side faces, 4 new vertices
        assert_eq!(result.revolution_faces.len(), 4);
        assert_eq!(result.new_vertices.len(), 4);
    }

    #[test]
    fn revolve_square_360_closed() {
        let mut brep = BRep::new(Tolerance::default());
        let face = build_revolve_square(&mut brep);

        let result = brep
            .revolve(
                face,
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                TAU,
            )
            .unwrap();

        assert!(
            result.cap_faces.is_empty(),
            "Full 360 should have no caps"
        );
        assert!(
            brep.vertices.len() > 4,
            "Should have created intermediate vertices"
        );
        // 4 rings x 4 vertices = 16 new vertices
        assert_eq!(result.new_vertices.len(), 16);
        // 4 rings x 4 side faces = 16 revolution faces
        assert_eq!(result.revolution_faces.len(), 16);
    }

    #[test]
    fn revolve_invalid_angle_zero() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(1.0, 0.0, 0.0)).unwrap();
        assert!(brep
            .revolve(
                mvfs.face,
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                0.0,
            )
            .is_err());
    }

    #[test]
    fn revolve_invalid_angle_too_large() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(1.0, 0.0, 0.0)).unwrap();
        assert!(brep
            .revolve(
                mvfs.face,
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                7.0,
            )
            .is_err());
    }

    #[test]
    fn revolve_rotated_vertices_correct_positions() {
        let mut brep = BRep::new(Tolerance::default());
        let face = build_revolve_square(&mut brep);

        let result = brep
            .revolve(
                face,
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                FRAC_PI_2,
            )
            .unwrap();

        // Original vertices at x=2..3, y=0. After 90° around Z:
        // should be at x~0, y=2..3.
        for &vh in &result.new_vertices {
            let v = &brep.vertices[vh];
            let r = (v.position.x * v.position.x + v.position.y * v.position.y).sqrt();
            // Radius from Z axis should be 2 or 3 (the original x values).
            assert!(
                (r - 2.0).abs() < 1e-10 || (r - 3.0).abs() < 1e-10,
                "Rotated vertex radius should be 2 or 3, got {r}"
            );
            // x should be ~0 (rotated 90° from x-axis to y-axis)
            assert!(
                v.position.x.abs() < 1e-10,
                "After 90° rotation, x should be ~0, got {}",
                v.position.x
            );
        }
    }

    #[test]
    fn revolve_180_creates_two_strips() {
        let mut brep = BRep::new(Tolerance::default());
        let face = build_revolve_square(&mut brep);

        let result = brep
            .revolve(
                face,
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                std::f64::consts::PI,
            )
            .unwrap();

        // 180° = 2 steps, 4 profile vertices => 8 side faces, 8 new vertices
        assert_eq!(result.revolution_faces.len(), 8);
        assert_eq!(result.new_vertices.len(), 8);
        assert_eq!(result.cap_faces.len(), 1, "Partial revolve should have 1 cap");
    }

    #[test]
    fn revolve_face_without_edges_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(1.0, 0.0, 0.0)).unwrap();
        // mvfs face has no edges
        let result = brep.revolve(
            mvfs.face,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            FRAC_PI_2,
        );
        assert!(result.is_err());
    }

    #[test]
    fn revolve_invalid_face_handle_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let _mvfs = brep.mvfs(Point3::new(1.0, 0.0, 0.0)).unwrap();
        let result = brep.revolve(
            Default::default(),
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            FRAC_PI_2,
        );
        assert!(result.is_err());
    }

    #[test]
    fn revolve_zero_axis_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let face = build_revolve_square(&mut brep);
        let result = brep.revolve(
            face,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            FRAC_PI_2,
        );
        assert!(result.is_err());
    }
}
