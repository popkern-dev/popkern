//! Boundary-ring capping — close a hole (a closed ring of free
//! boundary edges) with one new face.
//!
//! A shell hole in popkern is a set of edges in `BRep::boundary_edges`
//! whose surviving half-edge has a null-sentinel twin. This operator
//! manufactures the missing hole-side twin on every ring edge, wires
//! them into one new loop owned by one new face, and removes the edges
//! from `boundary_edges`. Once the set empties, the strict manifold /
//! Euler-Poincaré checks re-arm — a capped body is fully validated,
//! not grandfathered.
//!
//! Delta: V 0, E 0, HE +N, L +1, F +1, S 0, Sol 0, boundary_edges −N.
//! χ effect: ΔV − ΔE + ΔF = +1, the closure of one boundary component.

use crate::error::CoreError;
use crate::euler::contract::{EulerDelta, EulerTopologyCounts, verify_loop_contains_half_edges};
use crate::topology::{
    BRep, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, ShellHandle,
    VertexHandle,
};

/// Result of capping a closed boundary ring with a single face.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundaryRingCapResult {
    /// The new cap face (surfaceless; geometry is the caller's next step).
    pub face: FaceHandle,
    /// The new cap loop.
    pub loop_: LoopHandle,
    /// Hole-side half-edges, one per ring edge, in committed cycle
    /// order: `half_edges[k].next == half_edges[k + 1]`.
    pub half_edges: Vec<HalfEdgeHandle>,
}

/// One validated ring edge in cap-cycle order.
#[derive(Debug, Clone, Copy)]
pub(crate) struct BoundaryRingStep {
    /// The boundary edge being consumed.
    pub(crate) edge: EdgeHandle,
    /// The surviving (face-side) half-edge whose twin is the null sentinel.
    pub(crate) surviving: HalfEdgeHandle,
    /// Origin of the hole-side half-edge (end vertex of the surviving side).
    pub(crate) hole_origin: VertexHandle,
}

/// Validated capping plan: cycle-ordered steps plus the owning shell.
#[derive(Debug, Clone)]
pub(crate) struct BoundaryRingPlan {
    pub(crate) steps: Vec<BoundaryRingStep>,
    pub(crate) shell: ShellHandle,
}

impl BRep {
    /// Cap a closed ring of free boundary edges with one new face.
    ///
    /// `ring` is treated as a set — the operator derives the cycle
    /// order itself from vertex identity, starting at `ring[0]`. The
    /// cap loop traverses the hole opposite to the surrounding face
    /// loops, so winding stays outward-consistent by construction.
    ///
    /// The new face is surfaceless with no pcurves; attaching geometry
    /// and re-establishing coherence is the caller's next step.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if a ring edge or its surviving
    ///   half-edge is not live.
    /// - [`CoreError::PreconditionFailed`] if the ring has fewer than
    ///   3 edges, repeats an edge, contains an edge outside
    ///   `boundary_edges`, contains a registered boundary edge whose
    ///   twin is unexpectedly live (corruption), does not close into a
    ///   single cycle by vertex identity, shares a vertex between two
    ///   ring positions (figure-eight), or spans more than one shell.
    pub fn cap_boundary_ring(
        &mut self,
        ring: &[EdgeHandle],
    ) -> Result<BoundaryRingCapResult, CoreError> {
        let plan = self.plan_boundary_ring_cap(ring)?;
        let before = EulerTopologyCounts::from_brep(self);
        let mut candidate = self.clone();
        let result = candidate.apply_boundary_ring_cap(&plan);
        candidate.verify_boundary_ring_cap_postconditions(before, &plan, &result)?;
        *self = candidate;
        Ok(result)
    }

    /// Validate a boundary ring and derive its cap-cycle order.
    pub(crate) fn plan_boundary_ring_cap(
        &self,
        ring: &[EdgeHandle],
    ) -> Result<BoundaryRingPlan, CoreError> {
        if ring.len() < 3 {
            return Err(CoreError::PreconditionFailed(format!(
                "boundary ring cap needs at least 3 edges, got {}",
                ring.len()
            )));
        }
        let mut seen = std::collections::BTreeSet::new();
        for &edge in ring {
            if !seen.insert(edge) {
                return Err(CoreError::PreconditionFailed(format!(
                    "boundary ring cap edge {edge:?} is listed twice"
                )));
            }
        }

        // Resolve each edge to its surviving half-edge and directed span.
        let mut spans = Vec::with_capacity(ring.len());
        let mut shell: Option<ShellHandle> = None;
        for &edge in ring {
            if !self.boundary_edges.contains(&edge) {
                return Err(CoreError::PreconditionFailed(format!(
                    "boundary ring cap edge {edge:?} is not a registered boundary edge"
                )));
            }
            let edge_data = self.edges.get(edge).ok_or_else(|| {
                CoreError::InvalidHandle(format!("boundary ring cap edge {edge:?} is not live"))
            })?;
            let surviving = edge_data.half_edge;
            let surviving_data = self.half_edges.get(surviving).ok_or_else(|| {
                CoreError::InvalidHandle(format!(
                    "boundary ring cap edge {edge:?} has no live surviving half-edge"
                ))
            })?;
            if self.half_edges.contains_key(surviving_data.twin) {
                return Err(CoreError::PreconditionFailed(format!(
                    "boundary ring cap edge {edge:?} is registered as boundary but its twin is live"
                )));
            }
            let end = self
                .half_edges
                .get(surviving_data.next)
                .map(|next| next.origin)
                .ok_or_else(|| {
                    CoreError::InvalidHandle(format!(
                        "boundary ring cap edge {edge:?} surviving half-edge has no live next"
                    ))
                })?;
            let owner_shell = self
                .loops
                .get(surviving_data.loop_)
                .and_then(|loop_data| self.faces.get(loop_data.face))
                .map(|face| face.shell)
                .ok_or_else(|| {
                    CoreError::InvalidHandle(format!(
                        "boundary ring cap edge {edge:?} surviving half-edge has no owning face"
                    ))
                })?;
            match shell {
                None => shell = Some(owner_shell),
                Some(existing) if existing != owner_shell => {
                    return Err(CoreError::PreconditionFailed(format!(
                        "boundary ring cap spans shells {existing:?} and {owner_shell:?}; \
                         capping is a single-shell operation"
                    )));
                }
                Some(_) => {}
            }
            // The hole-side half-edge runs end -> start.
            spans.push(BoundaryRingStep {
                edge,
                surviving,
                hole_origin: end,
            });
        }
        let shell = shell.ok_or_else(|| {
            CoreError::PreconditionFailed("boundary ring cap resolved no owning shell".into())
        })?;

        // Chain hole-side spans tip-to-tail by vertex identity. The
        // hole-side half-edge of span k ends at its surviving origin;
        // the next span must originate there.
        let mut by_hole_origin = std::collections::BTreeMap::new();
        for (index, span) in spans.iter().enumerate() {
            if by_hole_origin.insert(span.hole_origin, index).is_some() {
                return Err(CoreError::PreconditionFailed(format!(
                    "boundary ring cap vertex {:?} starts two hole-side edges (figure-eight ring)",
                    span.hole_origin
                )));
            }
        }
        let mut steps = Vec::with_capacity(spans.len());
        let mut visited = std::collections::BTreeSet::new();
        let mut cursor = 0usize;
        for _ in 0..spans.len() {
            if !visited.insert(cursor) {
                return Err(CoreError::PreconditionFailed(format!(
                    "boundary ring closes after {} edges without consuming all {} listed edges \
                     (disconnected ring components)",
                    steps.len(),
                    spans.len()
                )));
            }
            let span = spans[cursor];
            steps.push(span);
            let next_origin = self
                .half_edges
                .get(span.surviving)
                .map(|surviving| surviving.origin)
                .ok_or_else(|| {
                    CoreError::InvalidHandle(format!(
                        "boundary ring cap edge {:?} surviving half-edge vanished mid-plan",
                        span.edge
                    ))
                })?;
            cursor = *by_hole_origin.get(&next_origin).ok_or_else(|| {
                CoreError::PreconditionFailed(format!(
                    "boundary ring is not closed: no ring edge continues from vertex \
                     {next_origin:?}"
                ))
            })?;
        }
        if cursor != 0 {
            return Err(CoreError::PreconditionFailed(
                "boundary ring does not close into a single cycle over all listed edges".into(),
            ));
        }
        Ok(BoundaryRingPlan { steps, shell })
    }

    fn apply_boundary_ring_cap(&mut self, plan: &BoundaryRingPlan) -> BoundaryRingCapResult {
        let face = self.faces.insert(Face {
            outer_loop: Default::default(),
            inner_loops: vec![],
            shell: plan.shell,
            surface: None,
            name: None,
            tolerance: None,
            role: Default::default(),
        });
        let loop_ = self.loops.insert(Loop {
            half_edge: Default::default(),
            face,
        });
        let mut hole_sides = Vec::with_capacity(plan.steps.len());
        for step in &plan.steps {
            let hole_side = self.half_edges.insert(HalfEdge {
                origin: step.hole_origin,
                twin: step.surviving,
                next: Default::default(),
                prev: Default::default(),
                edge: step.edge,
                loop_,
                pcurve: None,
            });
            self.half_edges[step.surviving].twin = hole_side;
            self.boundary_edges.remove(&step.edge);
            hole_sides.push(hole_side);
        }
        let count = hole_sides.len();
        for (index, &hole_side) in hole_sides.iter().enumerate() {
            self.half_edges[hole_side].next = hole_sides[(index + 1) % count];
            self.half_edges[hole_side].prev = hole_sides[(index + count - 1) % count];
        }
        self.loops[loop_].half_edge = hole_sides[0];
        self.faces[face].outer_loop = loop_;
        BoundaryRingCapResult {
            face,
            loop_,
            half_edges: hole_sides,
        }
    }

    fn verify_boundary_ring_cap_postconditions(
        &self,
        before: EulerTopologyCounts,
        plan: &BoundaryRingPlan,
        result: &BoundaryRingCapResult,
    ) -> Result<(), CoreError> {
        let count = plan.steps.len() as isize;
        EulerDelta::verify(
            "cap_boundary_ring",
            before,
            EulerTopologyCounts::from_brep(self),
            EulerDelta {
                vertices: 0,
                edges: 0,
                half_edges: count,
                loops: 1,
                faces: 1,
                shells: 0,
                solids: 0,
                boundary_edges: -count,
            },
        )?;
        verify_loop_contains_half_edges(self, result.loop_, &result.half_edges)?;
        for (step, &hole_side) in plan.steps.iter().zip(&result.half_edges) {
            if self.boundary_edges.contains(&step.edge) {
                return Err(CoreError::PreconditionFailed(format!(
                    "cap_boundary_ring left edge {:?} in boundary_edges",
                    step.edge
                )));
            }
            let surviving_twin = self
                .half_edges
                .get(step.surviving)
                .map(|surviving| surviving.twin);
            let hole_twin = self.half_edges.get(hole_side).map(|hole| hole.twin);
            if surviving_twin != Some(hole_side) || hole_twin != Some(step.surviving) {
                return Err(CoreError::PreconditionFailed(format!(
                    "cap_boundary_ring twin pairing broken on edge {:?}",
                    step.edge
                )));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "cap_boundary_ring_tests.rs"]
mod cap_boundary_ring_tests;
