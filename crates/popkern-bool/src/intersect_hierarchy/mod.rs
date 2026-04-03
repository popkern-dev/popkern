//! Sub-face intersection hierarchy: VV → VE → EE → VF → EF.
//!
//! Computes all pairwise intersections between two BRep bodies below the
//! face-face level, proceeding from lowest to highest dimension. Results
//! are collected into a [`HierarchyResult`] consumed by the FF intersection
//! phase and by imprinting.

#[allow(dead_code)]
pub mod result;

mod edge_edge;
mod edge_face;
mod vertex_edge;
mod vertex_face;
mod vertex_vertex;

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::topology::{BRep, FaceHandle};

use crate::tolerance::BoolTolerance;

pub use result::HierarchyResult;

/// Compute the full sub-FF intersection hierarchy between two bodies.
///
/// Runs VV → VE → EE → VF → EF in order. Each level's results are
/// accumulated into the returned [`HierarchyResult`].
///
/// `aabbs_a` and `aabbs_b` are precomputed face AABBs used for EF filtering.
pub(crate) fn compute_hierarchy(
    a: &BRep,
    b: &BRep,
    aabbs_a: &[(FaceHandle, Aabb)],
    aabbs_b: &[(FaceHandle, Aabb)],
    tol: &BoolTolerance,
) -> HierarchyResult {
    let snap = tol.snap;

    // Level 1: Vertex-Vertex
    let vv_hits = vertex_vertex::detect_vertex_vertex(a, b, snap);

    // Level 2: Vertex-Edge
    let ve_hits = vertex_edge::detect_vertex_edge(a, b, snap);

    // Level 3: Edge-Edge
    let ee_hits = edge_edge::detect_edge_edge(a, b, &tol.core);

    // Level 4: Vertex-Face
    let vf_hits = vertex_face::detect_vertex_face(a, b, aabbs_a, aabbs_b, snap);

    // Level 5: Edge-Face
    let ef_hits = edge_face::detect_edge_face(a, b, aabbs_a, aabbs_b, &tol.core);

    HierarchyResult {
        vv_hits,
        ve_hits,
        ee_hits,
        vf_hits,
        ef_hits,
    }
}

#[cfg(test)]
#[path = "hierarchy_tests.rs"]
mod tests;
