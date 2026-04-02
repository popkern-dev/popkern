#![warn(missing_docs)]

//! popkern-sketch: Sketch-to-BRep bridge for the popkern CAD kernel.
//!
//! Converts solved 2D parametric sketches (from `popkern-constraint`) on 3D
//! planes into BRep faces (in `popkern-core`). The pipeline is:
//!
//! 1. Define a [`SketchPlane`] in 3D space
//! 2. Build and solve a 2D sketch using `popkern-constraint`
//! 3. Convert the sketch into a BRep face via [`SketchOnPlane::to_brep`]

pub(crate) mod brep_builder;
pub(crate) mod containment;
pub(crate) mod curve_convert;
pub(crate) mod edge_graph;
pub mod error;
pub(crate) mod loop_tracer;
pub(crate) mod nesting;
pub mod plane;
pub(crate) mod reference_import;
pub(crate) mod self_intersection;

pub use error::{SketchError, SketchWarning};
pub use plane::SketchPlane;

// Re-export frequently used types from dependencies for convenience.
pub use popkern_constraint::{Constraint, Sketch};
pub use popkern_core::tolerance::Tolerance;
pub use popkern_core::topology::{BRep, FaceHandle, LoopHandle};

use std::collections::HashMap;

use popkern_constraint::ElementHandle;
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};

/// A 2D parametric sketch placed on a 3D plane.
///
/// Wraps a [`Sketch`] and a [`SketchPlane`] together so that the
/// sketch-to-BRep conversion has all the information it needs.
pub struct SketchOnPlane {
    /// The 2D parametric sketch.
    pub sketch: Sketch,
    /// The 3D plane the sketch lives on.
    pub plane: SketchPlane,
    /// Names assigned to non-construction sketch elements.
    ///
    /// Populated by [`assign_element_names`](Self::assign_element_names)
    /// before calling [`to_brep`](Self::to_brep) with tracking enabled.
    pub element_names: HashMap<ElementHandle, EntityName>,
}

/// A face region in the output BRep, consisting of a face and its holes.
#[derive(Debug, Clone)]
pub struct FaceRegion {
    /// Handle to the front face in the output BRep.
    pub face: FaceHandle,
    /// Handles to inner loops (holes) on this face.
    pub holes: Vec<LoopHandle>,
}

/// Result of converting a sketch to a BRep.
pub struct SketchBrepResult {
    /// The constructed BRep sheet body.
    pub brep: BRep,
    /// One entry per face region (outer boundary + holes).
    pub regions: Vec<FaceRegion>,
    /// Non-fatal warnings produced during conversion.
    pub warnings: Vec<SketchWarning>,
}

impl SketchOnPlane {
    /// Create a new sketch on the given plane.
    ///
    /// The sketch is initially empty; add elements and constraints via
    /// the [`sketch`](Self::sketch) field.
    pub fn new(plane: SketchPlane) -> Self {
        Self {
            sketch: Sketch::new(),
            plane,
            element_names: HashMap::new(),
        }
    }

    /// Assign stable names to all non-construction sketch elements.
    ///
    /// Call this before [`to_brep`](Self::to_brep) when tracking is enabled.
    /// Each non-construction curve element receives a `Generated` name
    /// with the given `operation` id, indexed in iteration order.
    pub fn assign_element_names(&mut self, operation: OperationId) {
        self.element_names.clear();
        let curves = self.sketch.non_construction_curves();
        for (i, elem) in curves.into_iter().enumerate() {
            let name = EntityName {
                operation,
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: i as u32,
                },
                ancestors: Vec::new(),
            };
            self.element_names.insert(elem, name);
        }
    }

    /// Convert the solved sketch into a BRep sheet body.
    ///
    /// Runs the full pipeline:
    /// 1. Build edge graph (merge endpoints, filter degenerate edges)
    /// 2. Trace closed loops via smallest-angle traversal
    /// 3. Check for self-intersecting profiles
    /// 4. Classify nesting (faces vs. holes)
    /// 5. Build BRep topology with curves and surfaces
    ///
    /// When `tracking` is `Some`, assigns stable names to all BRep entities
    /// and records them in the tracking context. Call
    /// [`assign_element_names`](Self::assign_element_names) beforehand to
    /// propagate sketch element names as ancestors.
    ///
    /// # Errors
    ///
    /// - [`SketchError::NoClosedLoops`] if the sketch has no closed loops
    /// - [`SketchError::SelfIntersectingProfile`] if any loop crosses itself
    /// - [`SketchError::InvalidTopology`] if BRep construction fails
    pub fn to_brep(
        &self,
        tolerance: Tolerance,
        tracking: Option<&mut TrackingContext>,
    ) -> Result<SketchBrepResult, SketchError> {
        // 1. Build edge graph.
        let (graph, warnings) = edge_graph::build_edge_graph(&self.sketch, &tolerance);
        if graph.edges.is_empty() {
            return Err(SketchError::NoClosedLoops);
        }

        // 2. Trace loops.
        let loops = loop_tracer::trace_loops(&graph, &self.sketch);
        if loops.is_empty() {
            return Err(SketchError::NoClosedLoops);
        }

        // 3. Self-intersection check.
        self_intersection::check_self_intersection(&loops, &self.sketch, &tolerance)?;

        // 4. Nesting classification.
        let roots = nesting::classify_nesting(&loops, &self.sketch);
        if roots.is_empty() {
            return Err(SketchError::NoClosedLoops);
        }

        // 5. Build BRep.
        let (mut brep, built_regions) =
            brep_builder::build_brep(&roots, &loops, &self.sketch, &self.plane, &tolerance)?;

        // 6. Assign stable names when tracking is active.
        if let Some(ctx) = tracking {
            assign_brep_names(&mut brep, ctx);
        }

        let regions = built_regions
            .into_iter()
            .map(|r| FaceRegion {
                face: r.face,
                holes: r.holes,
            })
            .collect();

        Ok(SketchBrepResult {
            brep,
            regions,
            warnings,
        })
    }
}

#[cfg(test)]
#[path = "lib_tests.rs"]
mod lib_tests;

/// Assign stable names to all entities in a freshly-built sketch BRep.
///
/// Faces get `Generated` names, edges get `Generated` names, and
/// vertices get `Junction` names. All are recorded in the tracking context.
fn assign_brep_names(brep: &mut BRep, ctx: &mut TrackingContext) {
    use popkern_core::naming_criteria::{
        capture_edge, capture_face, capture_vertex, GeometricCriteria,
    };

    // Collect handles first so we can mutate the BRep in the loop.
    let face_handles: Vec<_> = brep.faces.keys().collect();
    for (i, fh) in face_handles.iter().enumerate() {
        let name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: Vec::new(),
        };
        let criteria = capture_face(brep, *fh).map(GeometricCriteria::Face);
        let name_id = brep.name_pool.intern(name.clone(), criteria);
        brep.faces[*fh].name = Some(name_id);
        ctx.record_generated(name);
    }

    // Edges: Generated
    let edge_handles: Vec<_> = brep.edges.keys().collect();
    for (i, eh) in edge_handles.iter().enumerate() {
        let name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Lateral,
                index: i as u32,
            },
            ancestors: Vec::new(),
        };
        let criteria = capture_edge(brep, *eh).map(GeometricCriteria::Edge);
        let name_id = brep.name_pool.intern(name.clone(), criteria);
        brep.edges[*eh].name = Some(name_id);
        ctx.record_generated(name);
    }

    // Vertices: Junction
    let vertex_handles: Vec<_> = brep.vertices.keys().collect();
    for (i, vh) in vertex_handles.iter().enumerate() {
        let name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Junction,
                index: i as u32,
            },
            ancestors: Vec::new(),
        };
        let criteria = capture_vertex(brep, *vh).map(GeometricCriteria::Vertex);
        let name_id = brep.name_pool.intern(name.clone(), criteria);
        brep.vertices[*vh].name = Some(name_id);
        ctx.record_generated(name);
    }
}
