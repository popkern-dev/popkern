//! Self-intersection validation check (#11).
//!
//! Detects when non-adjacent faces in a solid intersect each other.

use crate::error::ValidationError;
use crate::geometry::intersection::degeneracy::DefaultDegeneracyHandler;
use crate::geometry::intersection::intersect_surfaces_with_handler;
use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::surface::SurfaceEval;
use crate::tolerance::Tolerance;
use crate::topology::BRep;
use crate::validation::bvh::Bvh;
use std::collections::HashSet;

impl BRep {
    /// Check #11: no face intersects another non-adjacent face.
    ///
    /// The check proceeds in three stages:
    ///
    /// 1. Build a set of adjacent face-index pairs by walking shared half-edges.
    /// 2. Sample each face surface on a 5×5 grid to compute an AABB; build a
    ///    BVH from those AABBs for efficient spatial queries.
    /// 3. For each non-adjacent pair whose AABBs overlap (from BVH query), call
    ///    `intersect_surfaces_with_handler` with `DefaultDegeneracyHandler` and
    ///    a coarse tolerance (10× linear) to test for actual surface intersection.
    ///    Any intersections found are reported as a `SelfIntersection` error.
    pub fn validate_self_intersection(&self) -> Result<(), ValidationError> {
        let face_keys: Vec<_> = self.faces.keys().collect();
        let n = face_keys.len();
        if n < 2 {
            return Ok(());
        }

        // ------------------------------------------------------------------ //
        // 1. Build adjacency set (by face index pairs, not handles)           //
        // ------------------------------------------------------------------ //
        let mut adjacent: HashSet<(usize, usize)> = HashSet::new();

        for (i, &fk) in face_keys.iter().enumerate() {
            let face = &self.faces[fk];
            let start_he_handle = self.loops[face.outer_loop].half_edge;
            // Empty loop — skip
            if !self.half_edges.contains_key(start_he_handle) {
                continue;
            }
            let mut he_handle = start_he_handle;
            // Walk the outer loop; guard against cycles larger than the mesh.
            let bound = self.half_edges.len().max(1);
            for _ in 0..=bound {
                let he = match self.half_edges.get(he_handle) {
                    Some(h) => h,
                    None => break,
                };
                let twin_handle = he.twin;
                if let Some(twin) = self.half_edges.get(twin_handle) {
                    let twin_loop_handle = twin.loop_;
                    if let Some(twin_loop) = self.loops.get(twin_loop_handle) {
                        let twin_face = twin_loop.face;
                        if let Some(j) = face_keys.iter().position(|&k| k == twin_face)
                            && i != j
                        {
                            let pair = if i < j { (i, j) } else { (j, i) };
                            adjacent.insert(pair);
                        }
                    }
                }
                he_handle = he.next;
                if he_handle == start_he_handle {
                    break;
                }
            }
        }

        // ------------------------------------------------------------------ //
        // 2. Compute face AABBs via 5×5 sampling; build BVH                  //
        // ------------------------------------------------------------------ //
        let sample_params = [0.0, 0.25, 0.5, 0.75, 1.0];

        let aabbs: Vec<Option<Aabb>> = face_keys
            .iter()
            .map(|&fk| {
                let face = &self.faces[fk];
                face.surface.as_ref().and_then(|s| {
                    let ((u0, u1), (v0, v1)) = s.geometry.parameter_range();
                    let mut pts = Vec::with_capacity(25);
                    for &tp in &sample_params {
                        let u = u0 + tp * (u1 - u0);
                        for &tq in &sample_params {
                            let v = v0 + tq * (v1 - v0);
                            pts.push(s.geometry.point_at(u, v));
                        }
                    }
                    Aabb::from_points(&pts)
                })
            })
            .collect();

        // Build BVH items: only include faces that have valid AABBs.
        let bvh_items: Vec<(usize, Aabb)> = aabbs
            .iter()
            .enumerate()
            .filter_map(|(i, aabb)| aabb.as_ref().map(|a| (i, a.clone())))
            .collect();

        if bvh_items.len() < 2 {
            return Ok(());
        }

        let bvh = Bvh::build(&bvh_items);

        // Use coarse tolerance for the BVH query padding.
        let linear_tol = self.tolerance.linear;
        let coarse_tol = Tolerance::default()
            .with_linear(linear_tol * 10.0)
            .with_angular(self.tolerance.angular)
            .with_parametric(self.tolerance.parametric);

        // ------------------------------------------------------------------ //
        // 3. Test non-adjacent, AABB-overlapping pairs with intersect_surfaces //
        // ------------------------------------------------------------------ //
        let overlapping = bvh.find_overlapping_pairs(linear_tol);

        for (i, j) in overlapping {
            let pair = if i < j { (i, j) } else { (j, i) };
            if adjacent.contains(&pair) {
                continue;
            }

            let surf_i = match self.faces[face_keys[i]].surface.as_ref() {
                Some(s) => s,
                None => continue,
            };
            let surf_j = match self.faces[face_keys[j]].surface.as_ref() {
                Some(s) => s,
                None => continue,
            };

            let ((u0i, u1i), (v0i, v1i)) = surf_i.geometry.parameter_range();
            let ((u0j, u1j), (v0j, v1j)) = surf_j.geometry.parameter_range();

            // Skip faces with degenerate (zero-area) parameter ranges.
            if (u1i - u0i).abs() < f64::EPSILON || (v1i - v0i).abs() < f64::EPSILON {
                continue;
            }
            if (u1j - u0j).abs() < f64::EPSILON || (v1j - v0j).abs() < f64::EPSILON {
                continue;
            }

            let result = intersect_surfaces_with_handler(
                &surf_i.geometry,
                (u0i, u1i),
                (v0i, v1i),
                &surf_j.geometry,
                (u0j, u1j),
                (v0j, v1j),
                &coarse_tol,
                &DefaultDegeneracyHandler,
            );

            match result {
                Ok(ref r) if !r.intersections.is_empty() => {
                    return Err(ValidationError::SelfIntersection {
                        face_a: face_keys[i],
                        face_b: face_keys[j],
                    });
                }
                _ => {}
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::error::ValidationError;
    use crate::geometry::nurbs_surface::NurbsSurface;
    use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    /// Build a degree-1 bilinear NURBS patch (a flat quad) spanning a finite
    /// parameter region.  The four corner control points are given in the order
    /// `(u0,v0), (u0,v1), (u1,v0), (u1,v1)` (row-major, u varies in each row).
    fn bilinear_patch(p00: Point3, p01: Point3, p10: Point3, p11: Point3) -> SurfaceGeometry {
        SurfaceGeometry::Nurbs(NurbsSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![p00, p01, p10, p11],
            num_u: 2,
            num_v: 2,
            weights: vec![1.0, 1.0, 1.0, 1.0],
        })
    }

    #[test]
    fn cube_passes_full_self_intersection() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        assert!(cube.validate_self_intersection().is_ok());
    }

    #[test]
    fn prism_passes_full_self_intersection() {
        let profile = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
        ];
        let prism = BRep::prism(
            &profile,
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
            Tolerance::default(),
        )
        .unwrap();
        assert!(prism.validate_self_intersection().is_ok());
    }

    #[test]
    fn empty_brep_passes() {
        let brep = BRep::new(Tolerance::default());
        assert!(brep.validate_self_intersection().is_ok());
    }

    #[test]
    fn single_face_brep_passes() {
        // A brep with only one face cannot self-intersect.
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev3 = brep
            .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let _mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
        assert!(brep.validate_self_intersection().is_ok());
    }

    /// Faces without surface geometry are skipped gracefully.
    #[test]
    fn face_without_surface_skips_gracefully() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face_keys: Vec<_> = brep.faces.keys().collect();

        // Remove surface from one face.
        brep.faces[face_keys[0]].surface = None;

        // Should still pass — faces without geometry are skipped.
        assert!(brep.validate_self_intersection().is_ok());
    }

    /// Test that adjacent faces in a cube are not reported as self-intersecting.
    #[test]
    fn adjacent_faces_not_reported() {
        // The cube has geometry on all faces; no non-adjacent faces share a plane.
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        assert!(cube.validate_self_intersection().is_ok());
    }

    /// Test the adjacency walking with a valid cube.
    #[test]
    fn unit_cube_with_geometry_no_self_intersection() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        assert!(cube.validate_self_intersection().is_ok());
    }

    /// Faces whose surfaces have no geometry (all None) cannot yield any BVH pairs.
    #[test]
    fn all_faces_without_surfaces_pass() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        // Remove all surface geometry.
        for fk in brep.faces.keys().collect::<Vec<_>>() {
            brep.faces[fk].surface = None;
        }
        assert!(brep.validate_self_intersection().is_ok());
    }

    /// Build a minimal two-face BRep (no shared edges) and attach crossing NURBS
    /// bilinear patches so the self-intersection check fires the error path.
    ///
    /// Topology: two completely disconnected faces created via mvfs twice —
    /// they share no edges, so the adjacency set is empty.  `SurfaceGeometry::Nurbs`
    /// is used because it has finite parameter ranges (from the knot vectors),
    /// which means the AABB sampling in `validate_self_intersection` produces real
    /// numbers, enabling the BVH to pair the two faces.
    ///
    /// Patch A: flat quad in the XY plane (z=0), u∈[0,1], v∈[0,1].
    ///   corners: (0,0,0), (0,1,0), (1,0,0), (1,1,0).
    /// Patch B: flat quad in the XZ plane (y=0), u∈[0,1], v∈[0,1].
    ///   corners: (0,0,0), (0,0,1), (1,0,0), (1,0,1).
    ///
    /// Both patches share the segment x∈[0,1], y=0, z=0, so their AABBs overlap
    /// and `intersect_surfaces_with_handler` finds an intersection line.
    #[test]
    fn self_intersecting_faces_detected() {
        let mut brep = BRep::new(Tolerance::default());

        // Face A — isolated topology node; carries patch A geometry.
        let mvfs_a = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev_a1 = brep
            .mev_at_vertex(mvfs_a.vertex, mvfs_a.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev_a2 = brep
            .mev(mev_a1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev_a3 = brep
            .mev(mev_a2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let mef_a = brep.mef(mev_a3.half_edge_twin, mev_a1.half_edge).unwrap();

        // Face B — isolated topology node; carries patch B geometry.
        let mvfs_b = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev_b1 = brep
            .mev_at_vertex(mvfs_b.vertex, mvfs_b.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev_b2 = brep
            .mev(mev_b1.half_edge_twin, Point3::new(1.0, 0.0, 1.0))
            .unwrap();
        let mev_b3 = brep
            .mev(mev_b2.half_edge_twin, Point3::new(0.0, 0.0, 1.0))
            .unwrap();
        let mef_b = brep.mef(mev_b3.half_edge_twin, mev_b1.half_edge).unwrap();

        // Patch A: flat quad in the XY plane (z = 0).
        //   row-major order: (u0,v0), (u0,v1), (u1,v0), (u1,v1)
        brep.faces[mef_a.face].surface = Some(SurfaceOnFace {
            geometry: bilinear_patch(
                Point3::new(0.0, 0.0, 0.0), // u=0, v=0
                Point3::new(0.0, 1.0, 0.0), // u=0, v=1
                Point3::new(1.0, 0.0, 0.0), // u=1, v=0
                Point3::new(1.0, 1.0, 0.0), // u=1, v=1
            ),
            orientation: true,
        });
        // Patch B: flat quad in the XZ plane (y = 0).
        //   intersects patch A along x∈[0,1], y=0, z=0.
        brep.faces[mef_b.face].surface = Some(SurfaceOnFace {
            geometry: bilinear_patch(
                Point3::new(0.0, 0.0, 0.0), // u=0, v=0
                Point3::new(0.0, 0.0, 1.0), // u=0, v=1
                Point3::new(1.0, 0.0, 0.0), // u=1, v=0
                Point3::new(1.0, 0.0, 1.0), // u=1, v=1
            ),
            orientation: true,
        });

        let result = brep.validate_self_intersection();
        assert!(
            matches!(result, Err(ValidationError::SelfIntersection { .. })),
            "expected SelfIntersection error, got: {result:?}"
        );
    }

    /// When one face in a non-adjacent pair has `surface = None` the pair is
    /// silently skipped (the `None => continue` branch in the intersection loop).
    /// Build the same two-disconnected-face topology but leave face B's surface
    /// as None — the check must pass.
    #[test]
    fn non_adjacent_pair_with_one_missing_surface_passes() {
        let mut brep = BRep::new(Tolerance::default());

        // Face A — bilinear NURBS patch in XY plane.
        let mvfs_a = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev_a1 = brep
            .mev_at_vertex(mvfs_a.vertex, mvfs_a.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev_a2 = brep
            .mev(mev_a1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev_a3 = brep
            .mev(mev_a2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let mef_a = brep.mef(mev_a3.half_edge_twin, mev_a1.half_edge).unwrap();

        // Face B — disconnected, surface stays None.
        let mvfs_b = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev_b1 = brep
            .mev_at_vertex(mvfs_b.vertex, mvfs_b.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev_b2 = brep
            .mev(mev_b1.half_edge_twin, Point3::new(1.0, 0.0, 1.0))
            .unwrap();
        let mev_b3 = brep
            .mev(mev_b2.half_edge_twin, Point3::new(0.0, 0.0, 1.0))
            .unwrap();
        let mef_b = brep.mef(mev_b3.half_edge_twin, mev_b1.half_edge).unwrap();

        // Only face A gets a surface; face B stays None.
        brep.faces[mef_a.face].surface = Some(SurfaceOnFace {
            geometry: bilinear_patch(
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ),
            orientation: true,
        });
        let _ = mef_b; // face B surface stays None

        assert!(
            brep.validate_self_intersection().is_ok(),
            "pair with one missing surface should not be reported as self-intersecting"
        );
    }
}
