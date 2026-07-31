//! Boundary-ring fan capping — close a hole with a fan of triangular
//! faces meeting at a new pole vertex.
//!
//! Composition of sanctioned operators executed on a candidate clone:
//! `cap_boundary_ring` + `mev` (pole spur) + (N−1) × `mef` (sector
//! peel), verified as one aggregate delta and swapped in atomically.
//! The residual cap loop after the last `mef` is the final sector.
//!
//! Delta: V +1, E +N, HE +3N, L +N, F +N, S 0, Sol 0,
//! boundary_edges −N. χ effect: ΔV − ΔE + ΔF = +1.

use crate::error::CoreError;
use crate::euler::contract::{EulerDelta, EulerTopologyCounts, verify_loop_contains_half_edges};
use crate::geometry::Point3;
use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

/// Result of capping a closed boundary ring with a pole fan.
///
/// Ring positions are the committed cap-cycle order: rim half-edge k
/// runs from ring vertex `w_k` to `w_{k+1}` on the hole side.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundaryRingFanResult {
    /// The new pole vertex all sectors meet at.
    pub pole_vertex: VertexHandle,
    /// Sector faces in ring order: face k's rim is rim half-edge k.
    pub sector_faces: Vec<FaceHandle>,
    /// Spoke edges in ring order: spoke k runs ring vertex `w_k` ↔
    /// pole and is shared by sectors k−1 and k.
    pub spoke_edges: Vec<EdgeHandle>,
    /// Hole-side rim half-edge per sector, in ring order.
    pub rim_half_edges: Vec<HalfEdgeHandle>,
}

impl BRep {
    /// Cap a closed ring of free boundary edges with a fan of N
    /// triangular faces meeting at a new pole vertex at `pole`.
    ///
    /// Ring validation is identical to [`BRep::cap_boundary_ring`].
    /// All new faces are surfaceless and all new edges curveless —
    /// attaching sector surfaces, spoke curves, and pcurves is the
    /// caller's next step. The pole position is taken as given; a
    /// geometric degeneracy gate (pole coincident with a ring vertex)
    /// belongs to the caller, which owns the tolerance.
    ///
    /// # Errors
    ///
    /// - Everything [`BRep::cap_boundary_ring`] declines.
    /// - [`CoreError::PreconditionFailed`] if the aggregate Euler
    ///   delta or per-sector loop membership fails verification.
    pub fn cap_boundary_ring_fan(
        &mut self,
        ring: &[EdgeHandle],
        pole: Point3,
    ) -> Result<BoundaryRingFanResult, CoreError> {
        let before = EulerTopologyCounts::from_brep(self);
        let mut candidate = self.clone();
        let result = candidate.apply_boundary_ring_fan(ring, pole)?;
        candidate.verify_boundary_ring_fan_postconditions(before, &result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_boundary_ring_fan(
        &mut self,
        ring: &[EdgeHandle],
        pole: Point3,
    ) -> Result<BoundaryRingFanResult, CoreError> {
        let cap = self.cap_boundary_ring(ring)?;
        let rim = cap.half_edges;
        let count = rim.len();

        // Plant the pole spur at rim vertex w_0: he_out runs w_0 → pole,
        // he_in runs pole → w_0 with he_in.next == rim[0].
        let spur = self.mev(rim[0], pole)?;
        let pole_vertex = spur.vertex;
        let mut spoke_edges = Vec::with_capacity(count);
        spoke_edges.push(spur.edge);

        // Peel sector k−1 with mef(spoke half-edge pole → w_{k−1}, rim[k]).
        // Each mef mints the spoke edge pole ↔ w_k and the sector face
        // whose loop is {rim[k−1], spoke w_k → pole, spoke pole → w_{k−1}}.
        let mut sector_faces = Vec::with_capacity(count);
        let mut pole_side = spur.half_edge_twin;
        for &rim_half_edge in rim.iter().skip(1) {
            let sector = self.mef(pole_side, rim_half_edge)?;
            sector_faces.push(sector.face);
            spoke_edges.push(sector.edge);
            pole_side = sector.half_edge;
        }
        // The residual cap loop {he_out w_0 → pole, pole → w_{N−1},
        // rim[N−1]} is the last sector, still on the original cap face.
        sector_faces.push(cap.face);

        Ok(BoundaryRingFanResult {
            pole_vertex,
            sector_faces,
            spoke_edges,
            rim_half_edges: rim,
        })
    }

    fn verify_boundary_ring_fan_postconditions(
        &self,
        before: EulerTopologyCounts,
        result: &BoundaryRingFanResult,
    ) -> Result<(), CoreError> {
        let count = result.rim_half_edges.len() as isize;
        EulerDelta::verify(
            "cap_boundary_ring_fan",
            before,
            EulerTopologyCounts::from_brep(self),
            EulerDelta {
                vertices: 1,
                edges: count,
                half_edges: 3 * count,
                loops: count,
                faces: count,
                shells: 0,
                solids: 0,
                boundary_edges: -count,
            },
        )?;
        for (index, (&face, &rim_half_edge)) in result
            .sector_faces
            .iter()
            .zip(&result.rim_half_edges)
            .enumerate()
        {
            let loop_ = self
                .faces
                .get(face)
                .map(|face_data| face_data.outer_loop)
                .ok_or_else(|| {
                    CoreError::PreconditionFailed(format!(
                        "cap_boundary_ring_fan sector {index} face is not live"
                    ))
                })?;
            verify_loop_contains_half_edges(self, loop_, &[rim_half_edge])?;
            let mut cycle_len = 0usize;
            let start = self
                .loops
                .get(loop_)
                .map(|loop_data| loop_data.half_edge)
                .ok_or_else(|| {
                    CoreError::PreconditionFailed(format!(
                        "cap_boundary_ring_fan sector {index} loop is not live"
                    ))
                })?;
            let mut cursor = start;
            loop {
                cycle_len += 1;
                cursor = self
                    .half_edges
                    .get(cursor)
                    .map(|half_edge| half_edge.next)
                    .ok_or_else(|| {
                        CoreError::PreconditionFailed(format!(
                            "cap_boundary_ring_fan sector {index} loop walk hit a dead half-edge"
                        ))
                    })?;
                if cursor == start {
                    break;
                }
                if cycle_len > 3 {
                    break;
                }
            }
            if cycle_len != 3 {
                return Err(CoreError::PreconditionFailed(format!(
                    "cap_boundary_ring_fan sector {index} loop has {cycle_len} half-edges, \
                     expected a triangle"
                )));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "cap_boundary_ring_fan_tests.rs"]
mod cap_boundary_ring_fan_tests;
