//! BRep — the top-level boundary representation container.

use crate::geometry::edge_triple::BoundarySide;
use crate::geometry::pcurve::Pcurve;
use crate::naming::{NamePool, OperationLabels};
use crate::precision::{
    DegradationTracker, EntityHandle, OperationPrecisionReport, PrecisionHealth, PrecisionSummary,
};
use crate::tolerance::Tolerance;
use crate::topology::brep_face_use_trim_storage::FaceUseTrimKey;
use crate::topology::*;
use serde::{Deserialize, Serialize};
use slotmap::SlotMap;

/// A boundary representation model.
///
/// Owns all topological entities in typed arenas (SlotMaps).
/// All topology modifications go through Euler operators.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRep {
    /// All vertices.
    pub vertices: SlotMap<VertexHandle, Vertex>,
    /// All half-edges.
    pub half_edges: SlotMap<HalfEdgeHandle, HalfEdge>,
    /// All edges.
    pub edges: SlotMap<EdgeHandle, Edge>,
    /// All loops.
    pub loops: SlotMap<LoopHandle, Loop>,
    /// All faces.
    pub faces: SlotMap<FaceHandle, Face>,
    /// All shells.
    pub shells: SlotMap<ShellHandle, Shell>,
    /// All solids.
    pub solids: SlotMap<SolidHandle, Solid>,
    /// Geometric tolerance for this model.
    pub tolerance: Tolerance,
    /// Pool of stable topological names for entities in this model.
    #[serde(default)]
    pub name_pool: NamePool,
    /// Human-readable labels for modelling operations.
    #[serde(default)]
    pub operation_labels: OperationLabels,
    /// Maximum imprecision allowed before the model is considered invalid.
    #[serde(default = "BRep::default_max_imprecision")]
    pub max_imprecision: f64,
    /// Edges with a null-sentinel twin (one face only). When
    /// non-empty, manifold / Euler validation is relaxed.
    #[serde(default)]
    pub boundary_edges: std::collections::BTreeSet<EdgeHandle>,
    /// Radial cycle for N>2 non-manifold edges. Empty for manifold
    /// BReps; [`BRep::radial_next`] falls back to `he.twin`.
    /// Serialized as key-sorted pairs: handle-keyed maps break JSON
    /// (non-string keys) and unsorted dumps break fingerprint
    /// determinism.
    #[serde(default, with = "super::serde_handle_map")]
    pub nonmanifold_radial: std::collections::HashMap<HalfEdgeHandle, HalfEdgeHandle>,
    /// Caller-authored evidence that a half-edge is a canonical boundary of
    /// its owning face surface. This is semantic construction data, not an
    /// inferred repair result: downstream repair/export stages can consume it
    /// before falling back to lossy pcurve or endpoint inference.
    #[serde(default, with = "super::serde_handle_map")]
    pub surface_boundary_hints: std::collections::HashMap<HalfEdgeHandle, BoundarySide>,
    /// Construction-time reasons a live half-edge was inspected for canonical
    /// boundary ownership but could not be accepted as one.
    #[serde(default, with = "super::serde_handle_map")]
    pub surface_boundary_declines:
        std::collections::HashMap<HalfEdgeHandle, SurfaceBoundaryHintDecline>,
    /// Durable audit trail for surface-boundary classification and shared
    /// trim-edge validation decisions. Repair knowledge collectors promote
    /// these records into structured facts so construction diagnostics are not
    /// hidden behind local trace output.
    #[serde(default)]
    pub surface_boundary_audit_events: Vec<SurfaceBoundaryAuditEvent>,
    /// Face-use-local pcurve authority for loop uses.
    #[serde(default)]
    pub face_use_pcurves: std::collections::BTreeMap<FaceUseTrimKey, Pcurve>,
    /// Tracks precision degradation across operations.
    #[serde(skip)]
    degradation_tracker: Option<DegradationTracker>,
    /// Transient L1/L3 obligations recorded via
    /// [`BRep::record_pending_authoring`]. MUST be empty when the
    /// BRep returns to a steady-state boundary (operator return,
    /// validation pass entry, BlendResult-shaped struct).
    #[serde(skip)]
    pub pending_coherence: Vec<crate::topology::coherence::PendingCoherenceMarker>,
    /// Source location of each edge's creator (file, line). Populated
    /// by [`BRep::insert_edge_tracked`] via `#[track_caller]`. Read by
    /// the Stage C0 coherence-attribution probe to map L1 violations
    /// to the operator that authored the edge. Edges inserted via the
    /// raw `edges.insert()` path will be absent from this map.
    #[serde(skip)]
    pub edge_creator_locations: std::collections::HashMap<EdgeHandle, (&'static str, u32)>,
}

impl BRep {
    fn default_max_imprecision() -> f64 {
        Tolerance::default().linear() * 100.0
    }

    /// Create a new empty BRep with the given tolerance.
    pub fn new(tolerance: Tolerance) -> Self {
        let max_imprecision = tolerance.linear() * 100.0;
        Self {
            vertices: SlotMap::with_key(),
            half_edges: SlotMap::with_key(),
            edges: SlotMap::with_key(),
            loops: SlotMap::with_key(),
            faces: SlotMap::with_key(),
            shells: SlotMap::with_key(),
            solids: SlotMap::with_key(),
            tolerance,
            name_pool: NamePool::new(),
            operation_labels: OperationLabels::default(),
            max_imprecision,
            boundary_edges: std::collections::BTreeSet::new(),
            nonmanifold_radial: std::collections::HashMap::new(),
            surface_boundary_hints: std::collections::HashMap::new(),
            surface_boundary_declines: std::collections::HashMap::new(),
            surface_boundary_audit_events: Vec::new(),
            face_use_pcurves: std::collections::BTreeMap::new(),
            degradation_tracker: None,
            pending_coherence: Vec::new(),
            edge_creator_locations: std::collections::HashMap::new(),
        }
    }

    /// Insert an edge into the slot map, capturing the caller's
    /// source location for Stage C0 coherence-attribution probes.
    ///
    /// The `(file, line)` of the immediate non-`#[track_caller]`
    /// caller is recorded in [`Self::edge_creator_locations`]. Use
    /// this in place of `self.edges.insert(...)` at every edge
    /// creation site that should be visible in attribution
    /// histograms; the raw insert path stays usable for tests and
    /// out-of-band repair surgery that should not be attributed.
    #[track_caller]
    pub fn insert_edge_tracked(&mut self, edge: Edge) -> EdgeHandle {
        let handle = self.edges.insert(edge);
        let loc = core::panic::Location::caller();
        self.edge_creator_locations
            .insert(handle, (loc.file(), loc.line()));
        handle
    }
}

/// Machine-readable reason a construction stage declined to record canonical
/// surface-boundary ownership for a half-edge.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum SurfaceBoundaryHintDeclineReason {
    /// The half-edge had no pcurve to classify on its owning surface.
    MissingPcurve,
    /// The committed pcurve does not lie on a rectangular canonical boundary.
    NotOnCanonicalBoundary,
    /// Authored boundary evidence disagrees with the half-edge pcurve.
    AuthoredHintBoundaryMismatch,
    /// Authored boundary evidence disagrees with the topological endpoints.
    AuthoredHintEndpointMismatch,
}

/// Structured decline emitted when construction inspected a half-edge for
/// canonical surface-boundary ownership and declined it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SurfaceBoundaryHintDecline {
    /// Machine-readable reason.
    pub reason: SurfaceBoundaryHintDeclineReason,
    /// Best sampled deviation from any canonical boundary in owning UV space.
    pub best_boundary_deviation: f64,
    /// Best sampled canonical boundary, when one was measurable.
    pub best_boundary: Option<BoundarySide>,
    /// Maximum gap between topology endpoints and the authored boundary trace.
    pub endpoint_gap: f64,
}

impl BRep {
    /// Record that `half_edge` follows a canonical boundary of its owning face
    /// surface. The hint is accepted only for live half-edges so stale handles
    /// cannot be introduced by construction code.
    pub fn set_surface_boundary_hint(
        &mut self,
        half_edge: HalfEdgeHandle,
        boundary: BoundarySide,
    ) -> bool {
        if !self.half_edges.contains_key(half_edge) {
            return false;
        }
        self.surface_boundary_declines.remove(&half_edge);
        self.surface_boundary_hints.insert(half_edge, boundary);
        true
    }

    /// Remove authored surface-boundary evidence for a half-edge.
    pub fn clear_surface_boundary_hint(&mut self, half_edge: HalfEdgeHandle) {
        self.surface_boundary_hints.remove(&half_edge);
    }

    /// Record why canonical boundary ownership was not accepted for a live
    /// half-edge inspected by construction or repair.
    pub fn set_surface_boundary_decline(
        &mut self,
        half_edge: HalfEdgeHandle,
        decline: SurfaceBoundaryHintDecline,
    ) -> bool {
        if !self.half_edges.contains_key(half_edge) {
            return false;
        }
        self.surface_boundary_hints.remove(&half_edge);
        self.surface_boundary_declines.insert(half_edge, decline);
        true
    }

    /// Return caller-authored canonical boundary evidence for a live half-edge.
    pub fn surface_boundary_hint(&self, half_edge: HalfEdgeHandle) -> Option<BoundarySide> {
        self.half_edges
            .contains_key(half_edge)
            .then(|| self.surface_boundary_hints.get(&half_edge).copied())
            .flatten()
    }

    /// Return the structured construction-time decline for a live half-edge.
    pub fn surface_boundary_decline(
        &self,
        half_edge: HalfEdgeHandle,
    ) -> Option<&SurfaceBoundaryHintDecline> {
        self.half_edges
            .contains_key(half_edge)
            .then(|| self.surface_boundary_declines.get(&half_edge))
            .flatten()
    }

    /// Record a durable surface-boundary audit event for repair knowledge and
    /// diagnostic consumers. Stale half-edge handles are rejected.
    pub fn record_surface_boundary_audit_event(
        &mut self,
        event: SurfaceBoundaryAuditEvent,
    ) -> bool {
        if !self.half_edges.contains_key(event.half_edge) {
            return false;
        }
        self.surface_boundary_audit_events.push(event);
        true
    }

    /// True when boundary edges exist (open shell / partial solid).
    pub fn has_boundary_edges(&self) -> bool {
        !self.boundary_edges.is_empty()
    }

    /// Collect vertex positions from the outer loop of a face by walking half-edges.
    ///
    /// Returns an empty `Vec` if the face handle is invalid, the loop is missing,
    /// or the loop has no half-edges.  This is the canonical implementation —
    /// downstream crates should call this instead of re-implementing the walk.
    pub fn face_outer_loop_positions(
        &self,
        face_handle: FaceHandle,
    ) -> Vec<crate::geometry::Point3> {
        let Some(face) = self.faces.get(face_handle) else {
            return Vec::new();
        };
        let Some(loop_) = self.loops.get(face.outer_loop) else {
            return Vec::new();
        };
        let start_he = loop_.half_edge;
        if !self.half_edges.contains_key(start_he) {
            return Vec::new();
        }
        let mut positions = Vec::new();
        let mut current = start_he;
        loop {
            let he = &self.half_edges[current];
            let vertex = &self.vertices[he.origin];
            positions.push(vertex.position);
            current = he.next;
            if current == start_he {
                break;
            }
        }
        positions
    }

    /// Collect sample points along the outer loop's edge curves.
    ///
    /// Walks the outer loop's half-edges. For each edge that has curve
    /// geometry, samples the curve at `n` equally-spaced parameter values
    /// and collects the 3D points. Includes vertex positions as well.
    ///
    /// This is more representative than [`face_outer_loop_positions`] for
    /// faces with self-loop edges (e.g., circular caps) where the vertex
    /// count is 1 but the curve traverses a full circle.
    pub fn face_boundary_sample_points(
        &self,
        face_handle: FaceHandle,
        samples_per_edge: usize,
    ) -> Vec<crate::geometry::Point3> {
        use crate::geometry::curve::CurveEval;

        let Some(face) = self.faces.get(face_handle) else {
            return Vec::new();
        };
        let Some(loop_) = self.loops.get(face.outer_loop) else {
            return Vec::new();
        };
        let start_he = loop_.half_edge;
        if !self.half_edges.contains_key(start_he) {
            return Vec::new();
        }

        let mut points = Vec::new();
        let mut current = start_he;
        loop {
            let he = &self.half_edges[current];
            // Always include the vertex position.
            points.push(self.vertices[he.origin].position);

            // If the edge has curve geometry, sample along it.
            let edge = &self.edges[he.edge];
            if let Some(ref curve) = edge.curve {
                let n = samples_per_edge.max(2);
                for i in 1..n {
                    let t = curve.t_start + (curve.t_end - curve.t_start) * (i as f64 / n as f64);
                    points.push(curve.geometry.point_at(t));
                }
            }

            current = he.next;
            if current == start_he {
                break;
            }
        }
        points
    }

    /// Get the effective tolerance for a vertex (entity-level or BRep default).
    pub fn vertex_tolerance(&self, vh: VertexHandle) -> f64 {
        self.vertices[vh]
            .tolerance
            .unwrap_or(self.tolerance.linear())
    }

    /// Get the effective tolerance for an edge (entity-level or BRep default).
    pub fn edge_tolerance(&self, eh: EdgeHandle) -> f64 {
        self.edges[eh].tolerance.unwrap_or(self.tolerance.linear())
    }

    /// Get the effective tolerance for a face (entity-level or BRep default).
    pub fn face_tolerance(&self, fh: FaceHandle) -> f64 {
        self.faces[fh].tolerance.unwrap_or(self.tolerance.linear())
    }

    /// Set entity tolerance — debug-asserts value >= current (monotonicity).
    pub fn set_entity_tolerance(&mut self, handle: EntityHandle, value: f64) {
        match handle {
            EntityHandle::Vertex(vh) => {
                debug_assert!(value >= self.vertex_tolerance(vh), "monotonicity");
                self.vertices[vh].tolerance = Some(value);
            }
            EntityHandle::Edge(eh) => {
                debug_assert!(value >= self.edge_tolerance(eh), "monotonicity");
                self.edges[eh].tolerance = Some(value);
            }
            EntityHandle::Face(fh) => {
                debug_assert!(value >= self.face_tolerance(fh), "monotonicity");
                self.faces[fh].tolerance = Some(value);
            }
        }
    }

    /// Intentionally tighten tolerance after recovery (only code path that reduces).
    pub fn recover_entity_tolerance(&mut self, handle: EntityHandle, value: f64) {
        match handle {
            EntityHandle::Vertex(vh) => self.vertices[vh].tolerance = Some(value),
            EntityHandle::Edge(eh) => self.edges[eh].tolerance = Some(value),
            EntityHandle::Face(fh) => self.faces[fh].tolerance = Some(value),
        }
    }

    /// Enable degradation tracking.
    pub fn enable_degradation_tracking(&mut self, window_size: usize) {
        self.degradation_tracker = Some(DegradationTracker::new(window_size));
    }

    /// Record an operation report.
    pub fn record_operation(&mut self, op_name: &'static str, report: &OperationPrecisionReport) {
        if let Some(ref mut tracker) = self.degradation_tracker {
            tracker.record(op_name, report);
        }
    }

    /// Compute aggregate precision health.
    pub fn precision_summary(&self) -> PrecisionSummary {
        let mut worst: Option<(EntityHandle, f64)> = None;
        let mut degraded_count = 0usize;
        let mut critical_count = 0usize;
        let degraded_thresh = self.max_imprecision * 0.5;
        let critical_thresh = self.max_imprecision * 0.9;

        // Scan vertices.
        for (vh, v) in &self.vertices {
            if let Some(tol) = v.tolerance {
                if tol >= critical_thresh {
                    critical_count += 1;
                } else if tol >= degraded_thresh {
                    degraded_count += 1;
                }
                let entity = EntityHandle::Vertex(vh);
                if worst.is_none_or(|(_, w)| tol > w) {
                    worst = Some((entity, tol));
                }
            }
        }
        // Scan edges.
        for (eh, e) in &self.edges {
            if let Some(tol) = e.tolerance {
                if tol >= critical_thresh {
                    critical_count += 1;
                } else if tol >= degraded_thresh {
                    degraded_count += 1;
                }
                let entity = EntityHandle::Edge(eh);
                if worst.is_none_or(|(_, w)| tol > w) {
                    worst = Some((entity, tol));
                }
            }
        }
        // Scan faces.
        for (fh, face) in &self.faces {
            if let Some(tol) = face.tolerance {
                if tol >= critical_thresh {
                    critical_count += 1;
                } else if tol >= degraded_thresh {
                    degraded_count += 1;
                }
                let entity = EntityHandle::Face(fh);
                if worst.is_none_or(|(_, w)| tol > w) {
                    worst = Some((entity, tol));
                }
            }
        }

        let health = if critical_count > 0 {
            PrecisionHealth::Critical
        } else if degraded_count > 0 {
            PrecisionHealth::Degraded
        } else {
            PrecisionHealth::Good
        };

        let health_score = match worst {
            Some((_, w)) => (1.0 - w / self.max_imprecision).max(0.0),
            None => 1.0,
        };

        let estimated_remaining = self.degradation_tracker.as_ref().and_then(|t| {
            let current = worst.map(|(_, w)| w).unwrap_or(0.0);
            t.estimated_operations_remaining(current, self.max_imprecision)
        });

        PrecisionSummary {
            health_score,
            health,
            worst_entity: worst,
            degraded_count,
            critical_count,
            estimated_remaining,
        }
    }

    /// List entities sorted by imprecision (worst first).
    ///
    /// Only entities with explicit (non-`None`) tolerance are included.
    pub fn degraded_entities(&self) -> Vec<(EntityHandle, f64)> {
        let mut result = Vec::new();
        for (vh, v) in &self.vertices {
            if let Some(tol) = v.tolerance {
                result.push((EntityHandle::Vertex(vh), tol));
            }
        }
        for (eh, e) in &self.edges {
            if let Some(tol) = e.tolerance {
                result.push((EntityHandle::Edge(eh), tol));
            }
        }
        for (fh, face) in &self.faces {
            if let Some(tol) = face.tolerance {
                result.push((EntityHandle::Face(fh), tol));
            }
        }
        result.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        result
    }
}

#[cfg(test)]
#[path = "brep_tests.rs"]
mod brep_tests;
#[path = "brep_transforms.rs"]
mod brep_transforms;
