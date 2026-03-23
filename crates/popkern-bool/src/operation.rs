//! Pipeline orchestrator: connects all 5 stages of the boolean pipeline.
//!
//! The [`execute_boolean`] function drives intersection, split, classify,
//! select, and merge in sequence, handling the no-intersection shortcut
//! and accumulating diagnostics on failure.

use popkern_core::{tolerance::Tolerance, topology::BRep};

use crate::classify::{classify_fragments, point_in_solid};
use crate::error::{BoolDiagnostics, BoolError, BoolOperation, PointClassification};
use crate::intersect::find_intersections;
use crate::merge::merge_faces;
use crate::select::{Operand, select_faces};
use crate::split::split_along_curves;

/// Execute a boolean operation on two BRep solids.
///
/// Runs the full five-stage pipeline: intersect, split, classify, select,
/// and merge. When no intersection curves are found, takes a shortcut based
/// on the spatial relationship between the two solids.
///
/// # Arguments
///
/// * `op` — the boolean set operation to perform.
/// * `a` — the first (left-hand) operand solid.
/// * `b` — the second (right-hand) operand solid.
/// * `tolerance` — geometric tolerance for all distance and coincidence tests.
///
/// # Errors
///
/// Returns [`BoolError`] describing the failure, wrapping partial diagnostics
/// from whichever stage failed.
pub(crate) fn execute_boolean(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    // Stage 1: Intersect — compute face-face intersection curves.
    let intersection_result =
        find_intersections(&a, &b, tolerance).map_err(|diag| diag.into_error())?;

    // No-intersection shortcut.
    if intersection_result.records.is_empty() {
        return handle_no_intersection(op, a, b, &intersection_result.aabbs_b, tolerance);
    }

    // Stage 2: Split — trim faces along intersection curves.
    let split_result = split_along_curves(a, b, &intersection_result.records, tolerance)
        .map_err(|diag| diag.into_error())?;

    // Stage 3: Classify — label each fragment as inside/outside/on-boundary.
    let classifications_a = classify_fragments(
        &split_result.a,
        &split_result.b,
        &intersection_result.aabbs_b,
        &split_result.fragments_a.intersection_edges,
        tolerance,
    )?;

    let classifications_b = classify_fragments(
        &split_result.b,
        &split_result.a,
        &intersection_result.aabbs_a,
        &split_result.fragments_b.intersection_edges,
        tolerance,
    )?;

    // Build the combined classification list for Stage 4.
    let mut all_classifications = Vec::new();
    for (face, class) in &classifications_a {
        all_classifications.push((*face, Operand::A, *class));
    }
    for (face, class) in &classifications_b {
        all_classifications.push((*face, Operand::B, *class));
    }

    // Stage 4: Select — choose which fragments belong to the result.
    let selection = select_faces(op, &all_classifications);

    // Stage 5: Merge — stitch selected fragments into a valid BRep.
    let result =
        merge_faces(&split_result, &selection, tolerance).map_err(|diag| diag.into_error())?;

    Ok(result)
}

/// Handle the case where no intersection curves were found between the solids.
///
/// Determines the spatial relationship (disjoint, containment, or coincident)
/// by classifying a sample point from A against B, then returns the
/// appropriate result for the given operation.
fn handle_no_intersection(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    aabbs_b: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    // Pick a vertex from A to classify against B.
    let sample_point = a
        .vertices
        .values()
        .next()
        .map(|v| v.position)
        .ok_or_else(|| BoolError::DegenerateInput("solid A has no vertices".to_string()))?;

    let classification = point_in_solid(&sample_point, &b, aabbs_b, tolerance)?;

    match classification {
        PointClassification::Outside => {
            // Disjoint: A is entirely outside B.
            handle_disjoint(op, a, b, tolerance)
        }
        PointClassification::Inside => {
            // Containment: A is entirely inside B.
            handle_containment(op, a, b, tolerance)
        }
        PointClassification::OnBoundary => {
            // Coincident boundary — ambiguous case.
            Err(BoolError::CoincidentFaces {
                face_a: a.faces.keys().next().unwrap_or_default(),
                face_b: b.faces.keys().next().unwrap_or_default(),
            })
        }
    }
}

/// Handle disjoint solids (A is entirely outside B) with no intersections.
fn handle_disjoint(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    match op {
        BoolOperation::Union => {
            // Combine both BReps into one.
            combine_breps(a, b)
        }
        BoolOperation::Intersection => {
            // Empty result.
            Ok(BRep::new(*tolerance))
        }
        BoolOperation::Subtraction => {
            // A minus nothing = A.
            Ok(a)
        }
    }
}

/// Handle containment (A is entirely inside B) with no intersections.
fn handle_containment(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    _tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    match op {
        BoolOperation::Union => {
            // B fully contains A, so union is just B.
            Ok(b)
        }
        BoolOperation::Intersection => {
            // A is inside B, so intersection is A.
            Ok(a)
        }
        BoolOperation::Subtraction => {
            // A − B when A is inside B would produce a hollow solid.
            // This requires KFMRH which is not yet implemented.
            Err(BoolError::MergeFailed {
                stage: "no_intersection".to_string(),
                reason: "subtraction with contained solid requires KFMRH (not yet implemented)"
                    .to_string(),
            })
        }
    }
}

/// Combine two disjoint BReps into a single BRep by copying all topology
/// from `b` into `a`.
fn combine_breps(mut a: BRep, b: BRep) -> Result<BRep, BoolError> {
    use popkern_core::topology::{Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};
    use std::collections::HashMap;

    // Build handle maps for B → new handles in A.
    let mut v_map: HashMap<
        popkern_core::topology::VertexHandle,
        popkern_core::topology::VertexHandle,
    > = HashMap::new();
    let mut he_map: HashMap<
        popkern_core::topology::HalfEdgeHandle,
        popkern_core::topology::HalfEdgeHandle,
    > = HashMap::new();
    let mut e_map: HashMap<popkern_core::topology::EdgeHandle, popkern_core::topology::EdgeHandle> =
        HashMap::new();
    let mut l_map: HashMap<popkern_core::topology::LoopHandle, popkern_core::topology::LoopHandle> =
        HashMap::new();
    let mut f_map: HashMap<popkern_core::topology::FaceHandle, popkern_core::topology::FaceHandle> =
        HashMap::new();
    let mut sh_map: HashMap<
        popkern_core::topology::ShellHandle,
        popkern_core::topology::ShellHandle,
    > = HashMap::new();
    let mut so_map: HashMap<
        popkern_core::topology::SolidHandle,
        popkern_core::topology::SolidHandle,
    > = HashMap::new();

    // Copy vertices.
    for (old_h, v) in &b.vertices {
        let new_h = a.vertices.insert(Vertex {
            position: v.position,
            outgoing: popkern_core::topology::HalfEdgeHandle::default(),
        });
        v_map.insert(old_h, new_h);
    }

    // Copy edges (with placeholder half_edge).
    for (old_h, e) in &b.edges {
        let new_h = a.edges.insert(Edge {
            half_edge: popkern_core::topology::HalfEdgeHandle::default(),
            curve: e.curve.clone(),
        });
        e_map.insert(old_h, new_h);
    }

    // Copy half-edges (with placeholders).
    for (old_h, _he) in &b.half_edges {
        let new_h = a.half_edges.insert(HalfEdge {
            origin: popkern_core::topology::VertexHandle::default(),
            twin: popkern_core::topology::HalfEdgeHandle::default(),
            next: popkern_core::topology::HalfEdgeHandle::default(),
            prev: popkern_core::topology::HalfEdgeHandle::default(),
            edge: popkern_core::topology::EdgeHandle::default(),
            loop_: popkern_core::topology::LoopHandle::default(),
        });
        he_map.insert(old_h, new_h);
    }

    // Copy loops (with placeholders).
    for (old_h, _l) in &b.loops {
        let new_h = a.loops.insert(Loop {
            half_edge: popkern_core::topology::HalfEdgeHandle::default(),
            face: popkern_core::topology::FaceHandle::default(),
        });
        l_map.insert(old_h, new_h);
    }

    // Copy solids (with placeholders).
    for (old_h, _s) in &b.solids {
        let new_h = a.solids.insert(Solid {
            shell: popkern_core::topology::ShellHandle::default(),
        });
        so_map.insert(old_h, new_h);
    }

    // Copy shells (with placeholders).
    for (old_h, _s) in &b.shells {
        let new_h = a.shells.insert(Shell {
            face: popkern_core::topology::FaceHandle::default(),
            solid: popkern_core::topology::SolidHandle::default(),
        });
        sh_map.insert(old_h, new_h);
    }

    // Copy faces (with placeholders).
    for (old_h, f) in &b.faces {
        let new_h = a.faces.insert(Face {
            outer_loop: popkern_core::topology::LoopHandle::default(),
            inner_loops: Vec::new(),
            shell: popkern_core::topology::ShellHandle::default(),
            surface: f.surface.clone(),
        });
        f_map.insert(old_h, new_h);
    }

    // Now fix up all the cross-references.

    // Fix half-edges.
    for (old_h, he) in &b.half_edges {
        let new_h = he_map[&old_h];
        a.half_edges[new_h].origin = v_map[&he.origin];
        a.half_edges[new_h].twin = he_map[&he.twin];
        a.half_edges[new_h].next = he_map[&he.next];
        a.half_edges[new_h].prev = he_map[&he.prev];
        a.half_edges[new_h].edge = e_map[&he.edge];
        a.half_edges[new_h].loop_ = l_map[&he.loop_];
    }

    // Fix vertices (outgoing).
    for (old_h, v) in &b.vertices {
        let new_h = v_map[&old_h];
        if b.half_edges.contains_key(v.outgoing) {
            a.vertices[new_h].outgoing = he_map[&v.outgoing];
        }
    }

    // Fix edges (half_edge).
    for (old_h, e) in &b.edges {
        let new_h = e_map[&old_h];
        a.edges[new_h].half_edge = he_map[&e.half_edge];
    }

    // Fix loops.
    for (old_h, l) in &b.loops {
        let new_h = l_map[&old_h];
        a.loops[new_h].half_edge = he_map[&l.half_edge];
        a.loops[new_h].face = f_map[&l.face];
    }

    // Fix faces.
    for (old_h, f) in &b.faces {
        let new_h = f_map[&old_h];
        a.faces[new_h].outer_loop = l_map[&f.outer_loop];
        a.faces[new_h].shell = sh_map[&f.shell];
        a.faces[new_h].inner_loops = f.inner_loops.iter().map(|il| l_map[il]).collect();
    }

    // Fix shells.
    for (old_h, s) in &b.shells {
        let new_h = sh_map[&old_h];
        a.shells[new_h].face = f_map[&s.face];
        a.shells[new_h].solid = so_map[&s.solid];
    }

    // Fix solids.
    for (old_h, s) in &b.solids {
        let new_h = so_map[&old_h];
        a.solids[new_h].shell = sh_map[&s.shell];
    }

    Ok(a)
}

impl BoolDiagnostics {
    /// Extract the error from diagnostics, or produce a generic error if none.
    fn into_error(self) -> BoolError {
        self.error.unwrap_or_else(|| BoolError::MergeFailed {
            stage: format!("{:?}", self.failed_stage),
            reason: "unknown error in boolean pipeline".to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    /// Translate all vertices of a BRep by a vector offset.
    fn translate(mut brep: BRep, offset: Vec3) -> BRep {
        for (_, vertex) in &mut brep.vertices {
            vertex.position = vertex.position + offset;
        }
        brep
    }

    #[test]
    fn disjoint_cubes_union_produces_valid_result() {
        let a = BRep::unit_cube(tol()).expect("unit cube a");
        let b = translate(
            BRep::unit_cube(tol()).expect("unit cube b"),
            Vec3::new(5.0, 0.0, 0.0),
        );
        let result = execute_boolean(BoolOperation::Union, a, b, &tol())
            .expect("disjoint union should succeed");
        // Should have 12 faces (6 + 6).
        assert_eq!(result.faces.len(), 12);
    }

    #[test]
    fn disjoint_cubes_intersection_is_empty() {
        let a = BRep::unit_cube(tol()).expect("unit cube a");
        let b = translate(
            BRep::unit_cube(tol()).expect("unit cube b"),
            Vec3::new(5.0, 0.0, 0.0),
        );
        let result = execute_boolean(BoolOperation::Intersection, a, b, &tol())
            .expect("disjoint intersection should succeed");
        assert_eq!(result.faces.len(), 0);
    }

    #[test]
    fn disjoint_cubes_subtraction_returns_a() {
        let a = BRep::unit_cube(tol()).expect("unit cube a");
        let b = translate(
            BRep::unit_cube(tol()).expect("unit cube b"),
            Vec3::new(5.0, 0.0, 0.0),
        );
        let result = execute_boolean(BoolOperation::Subtraction, a, b, &tol())
            .expect("disjoint subtraction should succeed");
        // Should return A unchanged — 6 faces.
        assert_eq!(result.faces.len(), 6);
    }

    #[test]
    fn contained_cube_union_returns_outer() {
        let a = BRep::unit_cube(tol()).expect("unit cube a");
        // Scale B to be much larger, containing A. We approximate this by
        // translating A so its first vertex is inside B.
        let mut b = BRep::unit_cube(tol()).expect("unit cube b");
        // Make B span [-1, 2]^3 by offsetting vertices.
        for (_, v) in &mut b.vertices {
            v.position = popkern_core::geometry::Point3::new(
                v.position.x * 3.0 - 1.0,
                v.position.y * 3.0 - 1.0,
                v.position.z * 3.0 - 1.0,
            );
        }
        let result = execute_boolean(BoolOperation::Union, a, b, &tol())
            .expect("contained union should succeed");
        // Union of contained A inside B should return B (6 faces).
        assert_eq!(result.faces.len(), 6);
    }

    #[test]
    fn contained_cube_intersection_returns_inner() {
        let a = BRep::unit_cube(tol()).expect("unit cube a");
        let mut b = BRep::unit_cube(tol()).expect("unit cube b");
        for (_, v) in &mut b.vertices {
            v.position = popkern_core::geometry::Point3::new(
                v.position.x * 3.0 - 1.0,
                v.position.y * 3.0 - 1.0,
                v.position.z * 3.0 - 1.0,
            );
        }
        let result = execute_boolean(BoolOperation::Intersection, a, b, &tol())
            .expect("contained intersection should succeed");
        // Intersection when A is inside B should return A (6 faces).
        assert_eq!(result.faces.len(), 6);
    }

    #[test]
    fn contained_cube_subtraction_returns_error() {
        let a = BRep::unit_cube(tol()).expect("unit cube a");
        let mut b = BRep::unit_cube(tol()).expect("unit cube b");
        for (_, v) in &mut b.vertices {
            v.position = popkern_core::geometry::Point3::new(
                v.position.x * 3.0 - 1.0,
                v.position.y * 3.0 - 1.0,
                v.position.z * 3.0 - 1.0,
            );
        }
        let result = execute_boolean(BoolOperation::Subtraction, a, b, &tol());
        // Should fail because KFMRH is not yet implemented.
        assert!(result.is_err());
    }
}
