//! Torus constructor — revolve-based with seam-closing post-process.

use crate::error::CoreError;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

#[path = "torus_pcurves.rs"]
mod torus_pcurves;
use torus_pcurves::{attach_torus_edge_curves, attach_torus_pcurves};

impl BRep {
    /// Build a torus by revolving a circular profile 360° around the Z axis,
    /// then closing the 0°/360° seam to produce a topologically correct
    /// genus-1 surface.
    ///
    /// The torus is centered at the origin. `major_radius` is the distance
    /// from the Z axis to the centre of the tube; `minor_radius` is the
    /// radius of the tube itself.
    ///
    /// # Topology
    ///
    /// V=16, E=32, F=16. Euler characteristic: V − E + F = 0 (genus 1).
    ///
    /// # Errors
    ///
    /// Returns [`CoreError::PreconditionFailed`] if either radius is
    /// non-positive, or if `minor_radius >= major_radius`.
    pub fn torus(
        major_radius: f64,
        minor_radius: f64,
        tolerance: Tolerance,
    ) -> Result<Self, CoreError> {
        if major_radius <= 0.0 || minor_radius <= 0.0 {
            return Err(CoreError::PreconditionFailed(
                "Radii must be positive".into(),
            ));
        }
        if minor_radius >= major_radius {
            return Err(CoreError::PreconditionFailed(
                "Minor radius must be less than major radius".into(),
            ));
        }

        let mut brep = BRep::new(tolerance);

        // Build a 4-vertex circular profile in the XZ plane centred at
        // (major_radius, 0, 0) with radius minor_radius.  The polygon loop
        // will visit vertices in the order [v0, v3, v2, v1] due to the
        // half-edge wiring produced by mev_at_vertex / mev / mef.
        let mvfs = brep.mvfs(Point3::new(major_radius + minor_radius, 0.0, 0.0))?;
        let mev1 = brep.mev_at_vertex(
            mvfs.vertex,
            mvfs.loop_,
            Point3::new(major_radius, 0.0, minor_radius),
        )?;
        let mev2 = brep.mev(
            mev1.half_edge_twin,
            Point3::new(major_radius - minor_radius, 0.0, 0.0),
        )?;
        let mev3 = brep.mev(
            mev2.half_edge_twin,
            Point3::new(major_radius, 0.0, -minor_radius),
        )?;
        let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge)?;

        // The polygon loop visits vertices in order: [v0, v3, v2, v1].
        // The revolve will create new_vertices[12..16] as the 360° copies of
        // [v0, v3, v2, v1] respectively.
        let orig_in_loop_order = [
            mvfs.vertex, // v0 — loop position 0
            mev3.vertex, // v3 — loop position 1
            mev2.vertex, // v2 — loop position 2
            mev1.vertex, // v1 — loop position 3
        ];

        // Revolve the profile face 360° around the Z axis.
        // Produces V=20, E=36, F=18 with χ=2 before seam-closing.
        let result = brep.revolve(
            mef.face,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            std::f64::consts::TAU,
        )?;

        // The last 4 new_vertices are the 360°-copies of orig_in_loop_order.
        if result.new_vertices.len() < 16 {
            return Err(CoreError::PreconditionFailed(
                "revolve did not produce expected vertex count".into(),
            ));
        }
        let vd: [VertexHandle; 4] = [
            result.new_vertices[12],
            result.new_vertices[13],
            result.new_vertices[14],
            result.new_vertices[15],
        ];

        // Close the 0°/360° seam to turn the open cylinder into a torus.
        close_torus_seam(&mut brep, mvfs.face, &orig_in_loop_order, &vd)?;

        // Assign torus surface geometry to every remaining face.
        let torus_surface = SurfaceGeometry::Torus {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            major_radius,
            minor_radius,
        };
        let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
        for fh in face_handles {
            brep.faces[fh].surface = Some(SurfaceOnFace {
                geometry: torus_surface.clone(),
                orientation: true,
            });
        }

        // Author exact arc edge curves, then pcurves (which anchor
        // their branch to the edge curve's midpoint).
        attach_torus_edge_curves(&mut brep, major_radius);
        attach_torus_pcurves(&mut brep)?;

        Ok(brep)
    }
}

/// Close the seam of a 360° revolve to produce a genus-1 (torus) topology.
///
/// After a full 360° revolve of a 4-vertex profile, the BRep is topologically
/// an open cylinder with two boundary loops:
/// - The "0° boundary": edges shared between ring-1 side faces and the mvfs
///   outer face.
/// - The "360° boundary": edges shared between ring-4 side faces and the
///   360° cap face (created but not returned by revolve).
///
/// This function stitches those two boundaries together:
/// 1. Merges each 360° duplicate vertex with its corresponding 0° original.
/// 2. Finds the matching boundary edges in the mvfs and cap face loops.
/// 3. Makes the ring-1 and ring-4 half-edges of each pair into twins.
/// 4. Removes the redundant half-edges, edges, faces, loops and vertices.
fn close_torus_seam(
    brep: &mut BRep,
    mvfs_face: FaceHandle,
    orig_verts: &[VertexHandle; 4],
    vd_verts: &[VertexHandle; 4],
) -> Result<(), CoreError> {
    // -----------------------------------------------------------------------
    // Step 1: redirect all half-edges that originate at vd_verts[i] to use
    //         orig_verts[i] instead, then remove vd_verts[i].
    // -----------------------------------------------------------------------
    for i in 0..4 {
        let vd = vd_verts[i];
        let vo = orig_verts[i];
        // Collect handles of half-edges that start at vd.
        let he_handles: Vec<HalfEdgeHandle> = brep
            .half_edges
            .iter()
            .filter(|(_, he)| he.origin == vd)
            .map(|(h, _)| h)
            .collect();
        for h in he_handles {
            brep.half_edges[h].origin = vo;
        }
        // Remove the duplicate vertex.
        brep.vertices.remove(vd);
    }

    // -----------------------------------------------------------------------
    // Step 2: find the 360° cap face.
    //
    // The cap face is the unique face (other than mvfs_face) whose outer loop
    // contains only half-edges originating at orig_verts (which now includes
    // the former vd_verts, after the merge above).  Since all vertices in the
    // BRep originate at orig_verts or the 90°/180°/270° positions, we
    // identify the cap face as the non-revolution, non-mvfs face whose outer
    // loop has exactly 4 half-edges, all with origins in orig_verts.
    //
    // In practice it is the face created as the cap of the last ring, which
    // is NOT included in result.revolution_faces (which has 16 entries) and
    // NOT the mvfs_face.  We find it by exhaustive search.
    // -----------------------------------------------------------------------
    let orig_set: std::collections::HashSet<VertexHandle> = orig_verts.iter().copied().collect();

    let cap_face = find_cap_face(brep, mvfs_face, &orig_set)
        .ok_or_else(|| CoreError::PreconditionFailed("could not locate 360° cap face".into()))?;

    // -----------------------------------------------------------------------
    // Step 3: collect the boundary half-edges from both faces.
    //
    // mvfs_boundary: half-edges in the mvfs face's outer loop.
    //   Each such he has its twin in some ring-1 side face.
    //
    // cap_boundary: half-edges in the cap face's outer loop.
    //   Each such he has its twin in some ring-4 side face.
    // -----------------------------------------------------------------------
    let mvfs_loop = brep
        .faces
        .get(mvfs_face)
        .ok_or_else(|| CoreError::InvalidHandle("mvfs face missing".into()))?
        .outer_loop;

    let cap_loop = brep
        .faces
        .get(cap_face)
        .ok_or_else(|| CoreError::InvalidHandle("cap face missing".into()))?
        .outer_loop;

    let mvfs_boundary = collect_loop_half_edges(brep, mvfs_loop)?;
    let cap_boundary = collect_loop_half_edges(brep, cap_loop)?;

    // -----------------------------------------------------------------------
    // Step 4: pair each mvfs boundary half-edge with the corresponding cap
    //         boundary half-edge (they connect the same vertex pair in
    //         opposite directions after the vertex merge), then stitch:
    //
    //   ring1_he  (twin of mvfs_he)  } become twins of each other,
    //   ring4_he  (twin of cap_he)   } replacing mvfs_he / cap_he.
    // -----------------------------------------------------------------------
    //
    // We collect (ring1_he, ring4_he, mvfs_he, cap_he, edge_to_remove) tuples
    // so that we can perform the stitching after building the full mapping.
    let mut stitches: Vec<(
        HalfEdgeHandle,
        HalfEdgeHandle,
        HalfEdgeHandle,
        HalfEdgeHandle,
    )> = Vec::with_capacity(4);

    for &he_m in &mvfs_boundary {
        let he_m_data = brep
            .half_edges
            .get(he_m)
            .ok_or_else(|| CoreError::InvalidHandle("mvfs boundary half-edge missing".into()))?;
        let a = he_m_data.origin;
        let b = brep.half_edges[brep.half_edges[he_m].next].origin;
        let ring1_he = he_m_data.twin;

        // Find the cap half-edge that goes from b→a (opposite direction).
        let he_c = cap_boundary
            .iter()
            .copied()
            .find(|&h| {
                let hd = &brep.half_edges[h];
                let c_next_origin = brep.half_edges[brep.half_edges[h].next].origin;
                hd.origin == b && c_next_origin == a
            })
            .ok_or_else(|| {
                CoreError::PreconditionFailed("no matching cap half-edge for mvfs boundary".into())
            })?;

        let ring4_he = brep.half_edges[he_c].twin;

        stitches.push((ring1_he, ring4_he, he_m, he_c));
    }

    // Apply stitches.
    for &(ring1_he, ring4_he, he_m, he_c) in &stitches {
        // Rewire twins.
        brep.half_edges[ring1_he].twin = ring4_he;
        brep.half_edges[ring4_he].twin = ring1_he;

        // The ring-4 half-edge inherits ring-1's edge handle so that one
        // canonical Edge entity survives for this seam edge.
        let surviving_edge = brep.half_edges[ring1_he].edge;
        let removed_edge = brep.half_edges[he_c].edge;
        brep.half_edges[ring4_he].edge = surviving_edge;

        // Keep the surviving edge's half_edge pointer valid.
        brep.edges[surviving_edge].half_edge = ring1_he;

        // Remove the stale boundary half-edges and their now-redundant edge.
        brep.remove_half_edge_reanchored(he_m);
        brep.remove_half_edge_reanchored(he_c);
        brep.edges.remove(removed_edge);
    }

    // -----------------------------------------------------------------------
    // Step 4b: the removed boundary half-edges may have been some merged
    // vertex's `outgoing` pointer (measured: v0 kept a dangling handle,
    // which silently exempted it from downstream per-vertex passes like
    // the shell vertex recompute). Re-anchor every seam vertex to a
    // SURVIVING half-edge that originates at it.
    // -----------------------------------------------------------------------
    for &vo in orig_verts {
        let outgoing = brep.vertices[vo].outgoing;
        if !brep.half_edges.contains_key(outgoing)
            && let Some((h, _)) = brep.half_edges.iter().find(|(_, he)| he.origin == vo)
        {
            brep.vertices[vo].outgoing = h;
        }
    }

    // -----------------------------------------------------------------------
    // Step 5: remove the two redundant loops and faces.
    // -----------------------------------------------------------------------
    let mvfs_loop_handle = walk_and_return_loop(brep, mvfs_face)?;
    let cap_loop_handle = walk_and_return_loop(brep, cap_face)?;

    brep.loops.remove(mvfs_loop_handle);
    brep.loops.remove(cap_loop_handle);
    brep.faces.remove(mvfs_face);
    brep.faces.remove(cap_face);

    Ok(())
}

/// Collect all half-edge handles in the loop in traversal order.
fn collect_loop_half_edges(
    brep: &BRep,
    loop_handle: LoopHandle,
) -> Result<Vec<HalfEdgeHandle>, CoreError> {
    let loop_data = brep
        .loops
        .get(loop_handle)
        .ok_or_else(|| CoreError::InvalidHandle("loop handle invalid".into()))?;

    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return Err(CoreError::PreconditionFailed(
            "loop has no half-edges".into(),
        ));
    }

    let mut result = Vec::new();
    let mut cursor = start;
    loop {
        result.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            break;
        }
        if result.len() > brep.half_edges.len() {
            return Err(CoreError::PreconditionFailed(
                "loop traversal did not close".into(),
            ));
        }
    }
    Ok(result)
}

/// Find the 360° cap face: the face (other than `mvfs_face`) whose outer loop
/// contains only half-edges originating at vertices in `orig_set`, and which
/// has exactly 4 half-edges.
fn find_cap_face(
    brep: &BRep,
    mvfs_face: FaceHandle,
    orig_set: &std::collections::HashSet<VertexHandle>,
) -> Option<FaceHandle> {
    for (fh, face_data) in &brep.faces {
        if fh == mvfs_face {
            continue;
        }
        let loop_he = brep.loops.get(face_data.outer_loop)?.half_edge;
        if !brep.half_edges.contains_key(loop_he) {
            continue;
        }

        // Walk the loop, check that all origins are in orig_set and count = 4.
        let mut count = 0usize;
        let mut all_orig = true;
        let mut cursor = loop_he;
        loop {
            let he = &brep.half_edges[cursor];
            if !orig_set.contains(&he.origin) {
                all_orig = false;
                break;
            }
            count += 1;
            cursor = he.next;
            if cursor == loop_he {
                break;
            }
            if count > 4 {
                break;
            }
        }

        if all_orig && count == 4 {
            return Some(fh);
        }
    }
    None
}

/// Return the outer loop handle of the given face (helper for removal).
fn walk_and_return_loop(brep: &BRep, face: FaceHandle) -> Result<LoopHandle, CoreError> {
    brep.faces.get(face).map(|f| f.outer_loop).ok_or_else(|| {
        CoreError::InvalidHandle("face handle invalid in walk_and_return_loop".into())
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;

    #[test]
    fn torus_creates_valid_topology() {
        let torus = BRep::torus(3.0, 1.0, Tolerance::default()).unwrap();
        assert!(!torus.vertices.is_empty());
        assert!(!torus.edges.is_empty());
        assert!(!torus.faces.is_empty());
        // Euler characteristic for genus-1: V - E + F = 0
        let v = torus.vertices.len() as i64;
        let e = torus.edges.len() as i64;
        let f = torus.faces.len() as i64;
        assert_eq!(
            v - e + f,
            0,
            "Torus should have Euler characteristic 0 (genus 1): V={v} E={e} F={f}"
        );
    }

    #[test]
    fn torus_vertices_on_surface() {
        let torus = BRep::torus(3.0, 1.0, Tolerance::default()).unwrap();
        for v in torus.vertices.values() {
            let p = &v.position;
            let dist_from_axis = (p.x * p.x + p.y * p.y).sqrt();
            let dist_from_major = ((dist_from_axis - 3.0).powi(2) + p.z * p.z).sqrt();
            assert!(
                (dist_from_major - 1.0).abs() < 0.1,
                "Vertex not on torus: dist={dist_from_major}"
            );
        }
    }

    #[test]
    fn torus_half_edges_have_pcurves() {
        let torus = BRep::torus(3.0, 1.0, Tolerance::default()).unwrap();
        for he_handle in torus.half_edges.keys() {
            assert!(
                torus.half_edges[he_handle].pcurve.is_some(),
                "torus half-edge {:?} missing pcurve",
                he_handle
            );
        }
    }

    #[test]
    fn torus_invalid_zero_radius() {
        assert!(BRep::torus(0.0, 1.0, Tolerance::default()).is_err());
        assert!(BRep::torus(1.0, 0.0, Tolerance::default()).is_err());
    }

    #[test]
    fn torus_minor_greater_than_major() {
        assert!(BRep::torus(1.0, 2.0, Tolerance::default()).is_err());
    }
}
