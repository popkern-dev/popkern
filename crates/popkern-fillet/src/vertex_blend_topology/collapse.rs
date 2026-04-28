//! Phase 6c-C.2c — lens-collapse Euler-op-like primitive.
//!
//! `collapse_apex_lens(brep, lens_face, gen1_he, gen2_he) -> Result<EdgeHandle, ...>`
//!
//! Atomically collapses a 3-edge "lens" face whose two generator
//! edges share **both** adjacent faces (i.e., gen1.twin and
//! gen2.twin are on the same other face — call it F_other).
//!
//! After the op:
//! * The lens face is gone.
//! * gen1 and gen2 are gone.
//! * The shared footprint vertex (where gen1 and gen2 meet) is
//!   gone.
//! * The lens's third edge survives.
//! * F_other's outer loop is spliced to skip the (now-deleted)
//!   gen1.twin and gen2.twin half-edges.
//! * The third edge's lens-side half-edge is re-parented to
//!   F_other's outer loop, replacing the spliced segment.
//!
//! Net effect: the lens region is absorbed by F_other, with the
//! third edge as F_other's new boundary toward F_third (the
//! face on the third edge's other side).
//!
//! # Snapshot-rollback
//!
//! Follows the same pattern as
//! `construct/topology_repair.rs::try_repair_chord_pair`: clone
//! the BRep before mutating; if `validate()` fails post-mutation,
//! restore the snapshot and return `Err`. This makes the op
//! transactional — no partial-state corruption escapes.
//!
//! # Scope
//!
//! **Same-twin-face case only** (gen1.twin.face == gen2.twin.face).
//! The cube_uniform cap-rebuild lens has gen1 and gen2 on
//! *different* cylinder faces (the cross-cylinder junction); that
//! case isn't handled by this primitive. Per design analysis,
//! the cross-cylinder pattern is an artifact of popkern's
//! current fillet pipeline (which produces redundant generators
//! at apex corners) — the commercial-grade fix is upstream
//! (have the fillet pipeline produce the cap-tangent arcs
//! directly), not a more elaborate downstream collapse op.
//!
//! For partial-cube and other apex-blend patterns where the
//! lens *does* have the simple shape, this op delivers the
//! Phase 6c-C.2c lens-collapse step end-to-end.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

/// Collapse a 3-edge lens face whose two generators share the
/// same opposite face. See module docs for full semantics.
///
/// Returns the surviving "third" edge handle on success.
#[allow(dead_code)]
pub(crate) fn collapse_apex_lens(
    brep: &mut BRep,
    lens_face: FaceHandle,
    gen1_he: HalfEdgeHandle,
    gen2_he: HalfEdgeHandle,
) -> Result<EdgeHandle, CollapseError> {
    // ── Preconditions ──────────────────────────────────────────
    let snapshot = brep.clone();

    let (gen1_lens_he, gen2_lens_he, third_lens_he, footprint_v, third_edge) =
        analyze_preconditions(brep, lens_face, gen1_he, gen2_he)?;

    // gen1.twin and gen2.twin must be in the SAME face's loop.
    let gen1_twin = brep.half_edges[gen1_lens_he].twin;
    let gen2_twin = brep.half_edges[gen2_lens_he].twin;
    let f_other_loop = brep.half_edges[gen1_twin].loop_;
    if brep.half_edges[gen2_twin].loop_ != f_other_loop {
        return Err(CollapseError::CrossFaceJunction);
    }
    let _f_other = brep.loops[f_other_loop].face;

    // The third edge must NOT be on the same face as F_other on
    // its other side (otherwise we'd be merging F_other with
    // itself — out of scope for this primitive).
    let third_twin = brep.half_edges[third_lens_he].twin;
    let f_third_loop = brep.half_edges[third_twin].loop_;
    if f_third_loop == f_other_loop {
        return Err(CollapseError::ThirdEdgeFolds);
    }

    // ── Mutation ──────────────────────────────────────────────
    //
    // Step 1: in F_other's loop, find the segment running
    // through the footprint vertex.  In our setup the loop
    // visits gen2.twin → ... → gen1.twin (the direction depends
    // on orientation; we determine empirically).
    //
    // The footprint vertex is the destination of one twin and
    // origin of the other. After surgery, the loop should
    // detour via the third edge's lens-side half-edge instead
    // of going through the footprint.
    //
    // F_other's loop currently:
    //   ... -> A -> gen?_twin -> footprint -> gen?_twin -> B -> ...
    // After surgery:
    //   ... -> A -> third_lens_he -> B -> ...
    let (entry_twin, exit_twin) = if brep.half_edges[gen1_twin].next == gen2_twin {
        (gen1_twin, gen2_twin)
    } else if brep.half_edges[gen2_twin].next == gen1_twin {
        (gen2_twin, gen1_twin)
    } else {
        return Err(CollapseError::TwinNotConsecutive);
    };

    let entry_prev = brep.half_edges[entry_twin].prev;
    let exit_next = brep.half_edges[exit_twin].next;

    // Re-parent third_lens_he to F_other's loop, between
    // entry_prev and exit_next.
    brep.half_edges[entry_prev].next = third_lens_he;
    brep.half_edges[third_lens_he].prev = entry_prev;
    brep.half_edges[third_lens_he].next = exit_next;
    brep.half_edges[exit_next].prev = third_lens_he;
    brep.half_edges[third_lens_he].loop_ = f_other_loop;

    // The third edge's twin (still in F_third's loop) needs
    // its origin/destination preserved — those don't change.
    // Origins of the third edge's two half-edges are V_a and
    // V_b respectively.  After re-parenting, third_lens_he's
    // origin is unchanged.

    // Update F_other's loop's anchor half-edge in case it
    // pointed to a soon-to-be-deleted twin.
    if brep.loops[f_other_loop].half_edge == entry_twin
        || brep.loops[f_other_loop].half_edge == exit_twin
    {
        brep.loops[f_other_loop].half_edge = third_lens_he;
    }

    // Step 2: delete gen1, gen2, the footprint vertex, the
    // lens face, the lens loop, and all 4 generator half-edges.
    let gen1_edge = brep.half_edges[gen1_lens_he].edge;
    let gen2_edge = brep.half_edges[gen2_lens_he].edge;
    let lens_loop = brep.faces[lens_face].outer_loop;

    brep.half_edges.remove(gen1_lens_he);
    brep.half_edges.remove(entry_twin);
    if entry_twin != gen1_twin {
        brep.half_edges.remove(gen1_twin);
    }
    brep.half_edges.remove(gen2_lens_he);
    if gen2_twin != entry_twin && gen2_twin != exit_twin {
        brep.half_edges.remove(gen2_twin);
    }
    if exit_twin != gen1_twin && exit_twin != gen2_twin {
        brep.half_edges.remove(exit_twin);
    }
    brep.edges.remove(gen1_edge);
    brep.edges.remove(gen2_edge);
    brep.vertices.remove(footprint_v);
    brep.loops.remove(lens_loop);
    brep.faces.remove(lens_face);

    // ── Validate ──────────────────────────────────────────────
    if brep.validate().is_err() {
        *brep = snapshot;
        return Err(CollapseError::ValidationFailed);
    }

    Ok(third_edge)
}

/// Verify the lens has exactly 3 outer-loop edges and identify
/// the third edge, the lens-side half-edges of the 3 edges, and
/// the footprint vertex shared by gen1 and gen2.
fn analyze_preconditions(
    brep: &BRep,
    lens_face: FaceHandle,
    gen1_he: HalfEdgeHandle,
    gen2_he: HalfEdgeHandle,
) -> Result<
    (
        HalfEdgeHandle,                       // gen1's lens-side he
        HalfEdgeHandle,                       // gen2's lens-side he
        HalfEdgeHandle,                       // third's lens-side he
        popkern_core::topology::VertexHandle, // footprint
        EdgeHandle,                           // third edge handle
    ),
    CollapseError,
> {
    let face_data = brep
        .faces
        .get(lens_face)
        .ok_or(CollapseError::InvalidLens)?;
    let lens_loop = face_data.outer_loop;
    let start = brep.loops[lens_loop].half_edge;

    let mut hes: Vec<HalfEdgeHandle> = Vec::new();
    let mut cursor = start;
    for _ in 0..8 {
        hes.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            break;
        }
    }
    if hes.len() != 3 {
        return Err(CollapseError::LensNotTriangle { actual: hes.len() });
    }

    if !hes.contains(&gen1_he) || !hes.contains(&gen2_he) {
        return Err(CollapseError::GenNotInLens);
    }
    if gen1_he == gen2_he {
        return Err(CollapseError::DuplicateGen);
    }

    let third_he = *hes
        .iter()
        .find(|&&h| h != gen1_he && h != gen2_he)
        .ok_or(CollapseError::GenNotInLens)?;

    // Footprint = the vertex shared by gen1's destination and
    // gen2's origin (or vice versa, depending on loop direction).
    let gen1_origin = brep.half_edges[gen1_he].origin;
    let gen1_dest = brep.half_edges[brep.half_edges[gen1_he].next].origin;
    let gen2_origin = brep.half_edges[gen2_he].origin;
    let gen2_dest = brep.half_edges[brep.half_edges[gen2_he].next].origin;
    let footprint = if gen1_dest == gen2_origin {
        gen1_dest
    } else if gen2_dest == gen1_origin {
        gen2_dest
    } else {
        return Err(CollapseError::GensDontMeetAtSharedVertex);
    };

    let third_edge = brep.half_edges[third_he].edge;
    Ok((gen1_he, gen2_he, third_he, footprint, third_edge))
}

/// Errors returned by `collapse_apex_lens`. All are recoverable
/// in the sense that the BRep is restored to its pre-call state
/// (snapshot rollback) — callers can fall back to a different
/// surgery path.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub(crate) enum CollapseError {
    /// `lens_face` is not a valid face handle.
    InvalidLens,
    /// `lens_face` doesn't have exactly 3 outer-loop edges.
    LensNotTriangle { actual: usize },
    /// One of `gen1_he`, `gen2_he` is not in `lens_face`'s loop.
    GenNotInLens,
    /// `gen1_he == gen2_he`.
    DuplicateGen,
    /// `gen1` and `gen2` don't share a vertex (origin/dest).
    GensDontMeetAtSharedVertex,
    /// `gen1.twin.face != gen2.twin.face` — the cross-face
    /// junction case. Out of scope for this primitive (see
    /// module docs).
    CrossFaceJunction,
    /// The third edge's twin is in F_other (the merge target),
    /// meaning the surgery would fold the third edge onto its
    /// own face. Not handled.
    ThirdEdgeFolds,
    /// `gen1.twin` and `gen2.twin` aren't consecutive in
    /// F_other's loop. The lens-collapse pattern requires them
    /// to be adjacent (separated only by the footprint vertex).
    TwinNotConsecutive,
    /// Post-mutation `validate()` failed; snapshot restored.
    ValidationFailed,
}
