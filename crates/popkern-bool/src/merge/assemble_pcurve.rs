//! Pcurve selection and direction-honest installation for the
//! graph-merge assembly (split from `assemble.rs` under the 500-line
//! structural budget).
//!
//! DIRECTION-HONEST INSTALL: graph half-edges store pcurves in the
//! RECORD's direction on BOTH sides (insert wires curve_uv1/curve_uv2
//! verbatim), while a brep half-edge's pcurve contract is its own
//! traversal direction. Convention-based reversal picked wrong
//! whenever the borrowing side and the record direction disagreed,
//! and the pcurve endpoint surgery in
//! `recompute_pcurves_where_missing` then "healed" the reversed
//! pcurve's end control points onto the OPPOSITE vertices — an
//! end-swapped pcurve whose interior runs backward, torn 2.9 at mesh
//! time (PK-485, the chamfered-puck mixed-sense rescue). Orientation
//! is decided by sampled 3D comparison — the same authority the
//! cellular assemble uses — never by ownership convention.

use popkern_core::constructors::BRepBuilder;
use popkern_core::topology::{BRep, HalfEdgeHandle};

use super::graph::{FaceRef, IntersectionGraph, Operand};

/// Install `pc` on `he`, oriented to the half-edge's traversal by
/// sampled 3D comparison against the face's surface.
pub(super) fn install_oriented(
    bld: &mut BRepBuilder,
    he: HalfEdgeHandle,
    pc: popkern_core::geometry::pcurve::Pcurve,
    fr: &FaceRef,
    brep_a: &BRep,
    brep_b: &BRep,
) {
    let brep = match fr.operand {
        Operand::A => brep_a,
        Operand::B => brep_b,
    };
    let surface = brep
        .faces
        .get(fr.face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone());
    let pc = crate::cellular::assemble_edge_geometry::orient_pcurve_to_half_edge(
        bld.brep(),
        he,
        pc,
        surface.as_ref(),
    );
    bld.set_half_edge_pcurve(he, pc);
}

/// The pcurve for `ghe` in the UV space of `face_ref`'s face: the
/// side whose registered face matches carries it; a twin borrow is
/// reversed into `ghe`'s walk direction (final authority stays the
/// sampled orientation in [`install_oriented`]). For shared
/// intersection edges the pcurve on each graph half-edge corresponds
/// to a specific face's UV space (forward -> face A, twin -> face B).
pub(super) fn correct_pcurve_for_face(
    graph: &IntersectionGraph,
    ghe: usize,
    face_ref: &FaceRef,
) -> Option<popkern_core::geometry::pcurve::Pcurve> {
    let ghe_face = &graph.half_edges[ghe].face;
    let twin = graph.half_edges[ghe].twin;
    let twin_face = &graph.half_edges[twin].face;

    // If this graph half-edge's face matches the fragment's face,
    // its pcurve is in the correct UV space.
    if ghe_face.as_ref() == Some(face_ref) {
        return graph.half_edges[ghe].pcurve.clone();
    }
    // If the twin's face matches, the twin's pcurve is in the
    // correct UV space but parameterized for the OPPOSITE walk
    // direction — reverse it. Verbatim borrowing silently flipped
    // closed rims (both rims of a boolean-rebuilt cylinder wall
    // wound the same way; endpoint matching cannot catch direction
    // on a closed curve).
    if twin_face.as_ref() == Some(face_ref) {
        return graph.half_edges[twin].pcurve.as_ref().map(|p| p.reversed());
    }
    // Neither matches (original edges, not shared) — use the graph
    // half-edge's own pcurve as-is.
    graph.half_edges[ghe].pcurve.clone()
}
