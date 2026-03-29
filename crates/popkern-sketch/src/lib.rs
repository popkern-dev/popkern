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
pub(crate) mod self_intersection;

pub use error::{SketchError, SketchWarning};
pub use plane::SketchPlane;

// Re-export frequently used types from dependencies for convenience.
pub use popkern_constraint::{Constraint, Sketch};
pub use popkern_core::tolerance::Tolerance;
pub use popkern_core::topology::{BRep, FaceHandle, LoopHandle};

/// A 2D parametric sketch placed on a 3D plane.
///
/// Wraps a [`Sketch`] and a [`SketchPlane`] together so that the
/// sketch-to-BRep conversion has all the information it needs.
pub struct SketchOnPlane {
    /// The 2D parametric sketch.
    pub sketch: Sketch,
    /// The 3D plane the sketch lives on.
    pub plane: SketchPlane,
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
    /// # Errors
    ///
    /// - [`SketchError::NoClosedLoops`] if the sketch has no closed loops
    /// - [`SketchError::SelfIntersectingProfile`] if any loop crosses itself
    /// - [`SketchError::InvalidTopology`] if BRep construction fails
    pub fn to_brep(&self, tolerance: Tolerance) -> Result<SketchBrepResult, SketchError> {
        // 1. Build edge graph.
        let (graph, warnings) = edge_graph::build_edge_graph(&self.sketch, &tolerance);
        if graph.edges.is_empty() {
            return Err(SketchError::NoClosedLoops);
        }

        // 2. Trace loops.
        let loops = loop_tracer::trace_loops(&graph);
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
        let (brep, built_regions) =
            brep_builder::build_brep(&roots, &loops, &self.sketch, &self.plane, &tolerance)?;

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
