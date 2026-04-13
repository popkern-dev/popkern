//! BRep — the top-level boundary representation container.

use crate::naming::{NamePool, OperationLabels};
use crate::precision::{
    DegradationTracker, EntityHandle, OperationPrecisionReport, PrecisionHealth, PrecisionSummary,
};
use crate::tolerance::Tolerance;
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
    /// Non-manifold radial cycle for half-edges whose parent edge is
    /// shared by more than two half-edges (N > 2 case, e.g. the
    /// crossing line of a transversal sheet-sheet union).
    ///
    /// For 2-manifold edges the cycle is the pair `(he, he.twin)` and
    /// this map is empty — [`BRep::radial_next`] falls back to
    /// `he.twin` when a half-edge is absent from the map. For
    /// non-manifold edges with N > 2 half-edges, the map contains
    /// every half-edge of the group and forms a full cyclic chain
    /// `h1 → h2 → h3 → ... → hN → h1`.
    #[serde(default)]
    pub nonmanifold_radial: std::collections::HashMap<HalfEdgeHandle, HalfEdgeHandle>,
    /// Tracks precision degradation across operations.
    #[serde(skip)]
    degradation_tracker: Option<DegradationTracker>,
}

impl BRep {
    /// Default max imprecision sentinel used during deserialization.
    fn default_max_imprecision() -> f64 {
        // Deserialization path — tolerance is not yet known, use a safe
        // fallback.  The caller should set this after deserialization if
        // the BRep's tolerance differs from the default.
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
            nonmanifold_radial: std::collections::HashMap::new(),
            degradation_tracker: None,
        }
    }
}

use crate::geometry::transform::Transform3;

impl BRep {
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

    /// Apply a rigid transformation to the entire BRep, returning a new copy.
    ///
    /// Transforms vertex positions, surface geometry on every face, edge
    /// curve geometry, and NURBS control points on pcurves (for non-planar
    /// faces; pcurves on planar faces live in UV space and are
    /// parameterization-invariant for rigid transforms that keep the plane
    /// canonical form, so they are regenerated via `attach_planar_pcurves`).
    pub fn transform(&self, t: &Transform3) -> Result<BRep, crate::error::CoreError> {
        let mut result = self.clone();

        // 1. Transform vertex positions.
        for (_, vertex) in &mut result.vertices {
            vertex.position = t.apply_point(&vertex.position);
        }

        // 2. Transform face surface geometry.
        for (_, face) in &mut result.faces {
            if let Some(ref mut sof) = face.surface {
                sof.geometry = transform_surface(&sof.geometry, t);
            }
        }

        // 3. Transform edge curve geometry.
        for (_, edge) in &mut result.edges {
            if let Some(ref mut curve_on_edge) = edge.curve {
                curve_on_edge.geometry = transform_curve(&curve_on_edge.geometry, t);
            }
        }

        // 4. Regenerate pcurves for planar faces (UV coordinates change
        //    when the plane origin/axes change).
        result.attach_planar_pcurves()?;

        Ok(result)
    }
}

/// Transform a `SurfaceGeometry` by a rigid transform.
fn transform_surface(
    surface: &crate::geometry::surface::SurfaceGeometry,
    t: &Transform3,
) -> crate::geometry::surface::SurfaceGeometry {
    use crate::geometry::surface::SurfaceGeometry;

    match surface {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        } => SurfaceGeometry::Plane {
            origin: t.apply_point(origin),
            normal: t.apply_direction(normal),
            u_axis: t.apply_direction(u_axis),
            v_axis: t.apply_direction(v_axis),
        },
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => SurfaceGeometry::Cylinder {
            origin: t.apply_point(origin),
            axis: t.apply_direction(axis),
            radius: *radius,
        },
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => SurfaceGeometry::Cone {
            apex: t.apply_point(apex),
            axis: t.apply_direction(axis),
            half_angle: *half_angle,
        },
        SurfaceGeometry::Sphere { center, radius } => SurfaceGeometry::Sphere {
            center: t.apply_point(center),
            radius: *radius,
        },
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => SurfaceGeometry::Torus {
            center: t.apply_point(center),
            axis: t.apply_direction(axis),
            major_radius: *major_radius,
            minor_radius: *minor_radius,
        },
        SurfaceGeometry::Nurbs(nurbs) => {
            let control_points = nurbs
                .control_points
                .iter()
                .map(|p| t.apply_point(p))
                .collect();
            SurfaceGeometry::Nurbs(crate::geometry::NurbsSurface {
                degree_u: nurbs.degree_u,
                degree_v: nurbs.degree_v,
                knots_u: nurbs.knots_u.clone(),
                knots_v: nurbs.knots_v.clone(),
                control_points,
                num_u: nurbs.num_u,
                num_v: nurbs.num_v,
                weights: nurbs.weights.clone(),
            })
        }
    }
}

/// Transform a `CurveGeometry` by a rigid transform.
fn transform_curve(
    curve: &crate::geometry::curve::CurveGeometry,
    t: &Transform3,
) -> crate::geometry::curve::CurveGeometry {
    use crate::geometry::curve::CurveGeometry;

    match curve {
        CurveGeometry::Line { origin, direction } => CurveGeometry::Line {
            origin: t.apply_point(origin),
            direction: t.apply_direction(direction),
        },
        CurveGeometry::Circle {
            center,
            normal,
            radius,
        } => CurveGeometry::Circle {
            center: t.apply_point(center),
            normal: t.apply_direction(normal),
            radius: *radius,
        },
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle,
            end_angle,
        } => CurveGeometry::Arc {
            center: t.apply_point(center),
            normal: t.apply_direction(normal),
            radius: *radius,
            start_angle: *start_angle,
            end_angle: *end_angle,
        },
        CurveGeometry::Nurbs(nurbs) => {
            let control_points = nurbs
                .control_points
                .iter()
                .map(|p| t.apply_point(p))
                .collect();
            CurveGeometry::Nurbs(crate::geometry::NurbsCurve {
                degree: nurbs.degree,
                knots: nurbs.knots.clone(),
                control_points,
                weights: nurbs.weights.clone(),
            })
        }
    }
}

#[cfg(test)]
#[path = "brep_tests.rs"]
mod brep_tests;
