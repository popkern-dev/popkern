//! Helper functions for the boolean pipeline orchestrator.
//!
//! Contains the no-intersection shortcut handlers (disjoint, containment,
//! reverse containment) and the BRep copy/combine utilities.

use std::collections::HashMap;

use popkern_core::{
    tolerance::Tolerance,
    topology::{BRep, FaceHandle},
};

use crate::classify::point_in_solid;
use crate::error::{BoolError, BoolOperation, Classification, PointClassification};
use crate::split::SplitResult;

/// Handle the case where no intersection curves were found between the solids.
///
/// Determines the spatial relationship (disjoint, containment, or coincident)
/// by classifying a sample point from A against B, then returns the
/// appropriate result for the given operation.
pub(super) fn handle_no_intersection(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    aabbs_a: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
    aabbs_b: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    // Pick a vertex from A to classify against B.
    let sample_point_a = a
        .vertices
        .values()
        .next()
        .map(|v| v.position)
        .ok_or_else(|| BoolError::DegenerateInput("solid A has no vertices".to_string()))?;

    let classification = point_in_solid(&sample_point_a, &b, aabbs_b, tolerance)?;

    let result = match classification {
        PointClassification::Outside => {
            // A's sample point is outside B. Check reverse: is B inside A?
            let sample_point_b = b
                .vertices
                .values()
                .next()
                .map(|v| v.position)
                .ok_or_else(|| BoolError::DegenerateInput("solid B has no vertices".to_string()))?;
            let reverse = point_in_solid(&sample_point_b, &a, aabbs_a, tolerance)?;
            match reverse {
                PointClassification::Inside => {
                    // B is inside A — reverse containment.
                    handle_reverse_containment(op, a, b, tolerance)
                }
                _ => {
                    // Truly disjoint.
                    handle_disjoint(op, a, b, tolerance)
                }
            }
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
    };

    // Log validation issues on the no-intersection result.
    if let Ok(ref brep) = result {
        super::log_validation_issues(brep, op);
    }

    result
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
    tolerance: &Tolerance,
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
            // A − B when A is inside B: A is entirely consumed by B, result is empty.
            Ok(BRep::new(*tolerance))
        }
    }
}

/// Handle reverse containment (B is entirely inside A) with no intersections.
fn handle_reverse_containment(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    match op {
        BoolOperation::Union => {
            // A fully contains B, so union is just A.
            Ok(a)
        }
        BoolOperation::Intersection => {
            // B is inside A, so intersection is B.
            Ok(b)
        }
        BoolOperation::Subtraction => {
            // A − B when B is inside A: hollow solid (A with B-shaped void).
            build_hollow_solid(a, b, tolerance)
        }
    }
}

/// Build a hollow solid: A with a B-shaped void inside.
///
/// Copies all of B's topology into A with flipped surface normals, creating
/// an inner void. The result has A's original outer faces plus B's faces
/// with inverted orientation.
fn build_hollow_solid(a: BRep, b: BRep, _tolerance: &Tolerance) -> Result<BRep, BoolError> {
    copy_brep_into(a, b, true)
}

/// Combine two disjoint BReps into a single BRep by copying all topology
/// from `b` into `a`.
///
/// The two inputs should not share geometry — this is a structural
/// concatenation that copies handles from `b` into fresh slotmap slots
/// in `a` without welding at intersection points. Used by the sheet
/// dispatcher to implement non-overlapping `SheetSheetUnion`.
pub(crate) fn combine_breps(a: BRep, b: BRep) -> Result<BRep, BoolError> {
    copy_brep_into(a, b, false)
}

/// Handle the case where the select stage produces 0 kept faces.
///
/// This happens when all fragments from both operands classify the same way
/// (e.g. all Inside). Infers containment from the classification maps and
/// returns the appropriate operand directly. When containment cannot be
/// determined, returns an empty BRep.
pub(super) fn handle_empty_selection(
    op: BoolOperation,
    split_result: SplitResult,
    classifications_a: &HashMap<FaceHandle, Classification>,
    classifications_b: &HashMap<FaceHandle, Classification>,
) -> BRep {
    let all_a_inside = !classifications_a.is_empty()
        && classifications_a
            .values()
            .all(|c| *c == Classification::Inside);
    let all_b_inside = !classifications_b.is_empty()
        && classifications_b
            .values()
            .all(|c| *c == Classification::Inside);

    if all_a_inside {
        // A is entirely inside B.
        return match op {
            BoolOperation::Union => split_result.b,
            BoolOperation::Intersection => split_result.a,
            BoolOperation::Subtraction => BRep::new(split_result.a.tolerance),
        };
    }

    if all_b_inside {
        // B is entirely inside A.
        return match op {
            BoolOperation::Union => split_result.a,
            BoolOperation::Intersection => split_result.b,
            // A minus B when B is inside A — full hollow solid requires
            // combining A with flipped B. For now, return A (the outer shell).
            BoolOperation::Subtraction => split_result.a,
        };
    }

    // Cannot determine containment — return empty result.
    BRep::new(split_result.a.tolerance)
}

/// Copy all topology from `b` into `a`, optionally flipping B's face normals.
///
/// When `flip_b_normals` is true, the surface orientation on every copied face
/// is inverted, producing an inner void shell for hollow-solid subtraction.
fn copy_brep_into(mut a: BRep, b: BRep, flip_b_normals: bool) -> Result<BRep, BoolError> {
    use popkern_core::topology::{Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};

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
            name: v.name,
            tolerance: None,
        });
        v_map.insert(old_h, new_h);
    }

    // Copy edges (with placeholder half_edge).
    for (old_h, e) in &b.edges {
        let new_h = a.edges.insert(Edge {
            half_edge: popkern_core::topology::HalfEdgeHandle::default(),
            curve: e.curve.clone(),
            name: e.name,
            tolerance: None,
            canonical: e.canonical,
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
            pcurve: None,
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

    // Copy faces, optionally flipping surface orientation.
    for (old_h, f) in &b.faces {
        let surface = if flip_b_normals {
            f.surface
                .as_ref()
                .map(|s| popkern_core::geometry::SurfaceOnFace {
                    geometry: s.geometry.clone(),
                    orientation: !s.orientation,
                })
        } else {
            f.surface.clone()
        };
        let new_h = a.faces.insert(Face {
            outer_loop: popkern_core::topology::LoopHandle::default(),
            inner_loops: Vec::new(),
            shell: popkern_core::topology::ShellHandle::default(),
            surface,
            name: f.name,
            tolerance: None,
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
        a.half_edges[new_h].pcurve = he.pcurve.clone();
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
