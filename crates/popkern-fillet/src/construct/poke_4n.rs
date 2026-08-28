//! 4N subdivision of an N-gonal blend face — KP13 bi-5 macro-patch
//! layout.
//!
//! The KP13 bi-5 cap produces `4N` bi-9 Gregory patches arranged in
//! a 2 × 2 piece layout per polar sector (paper's O, L, R, M pieces,
//! here labelled `C_star, O_L, O_C, O_R`). To land these as BRep
//! sub-faces, we need a poke that splits the N-gon into `4N` quads
//! (not `N`, as [`poke_n_gon_face_into_n_quads`] does).
//!
//! Algorithm (composes on top of the N-way poke):
//!
//! 1. Call [`poke_n_gon_face_into_n_quads`] → `N` sector quads + a
//!    new central vertex `C`.
//! 2. Split each of the `N` radial edges (`C ↔ V_{2k}`, shared
//!    between sectors `k` and `k+1`) at `t = 0.5` → `N` midpoint
//!    vertices `M_rad_k`.
//! 3. For each sector `k`:
//!    a. Split its 2 outer edges at `t = 0.5` → 2 midpoints
//!    `M_out_L, M_out_R` per sector (unique, not shared).
//!    b. Walk the sector's 8-HE loop starting from `C`.
//!    c. Insert a centroid vertex `X_k` inside the sector with
//!    `mev`, then run 3 `mef`s to cut `X_k` to the other 3
//!    midpoints — producing 4 quads.
//!
//! Euler invariant:
//!   ΔV = +4N (N radial mids + 2N outer mids + N centroids)
//!   ΔE = +7N (N radial splits + 2N outer splits + 4N interior)
//!   ΔF = +3N (each sector's 1 face → 4 faces)
//!   ΔV − ΔE + ΔF = 4N − 7N + 3N = 0 ✓
//!
//! # Piece labelling
//!
//! Each sector's 4 sub-quads land in the returned `sub_faces` in
//! canonical `(sector, piece)` row-major order:
//!
//! ```text
//!   sub_faces[4k + 0] = O_L     (left-outer, at corner V_{2k-2})
//!   sub_faces[4k + 1] = O_C     (outer-center, at midpoint V_{2k-1})
//!   sub_faces[4k + 2] = O_R     (right-outer, at corner V_{2k})
//!   sub_faces[4k + 3] = C_star  (star-side, at central vertex C)
//! ```
//!
//! This order matches the `(sector, piece)` layout produced by
//! [`try_kp_bi5_macro_cap_as_gregory`] — see the module header of
//! `kp_bi5_macro_apply.rs`.

#![allow(dead_code)]

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::poke_hex::poke_n_gon_face_into_n_quads;

/// Result of [`poke_n_gon_face_into_4n_quads`].
#[derive(Debug)]
pub(crate) struct HexPoke4nResult {
    /// `4N` sub-face handles in `(sector, piece)` row-major order.
    /// See module header for the piece labelling.
    pub sub_faces: Vec<FaceHandle>,
    /// The central vertex `C` created by the inner N-way poke —
    /// shared by all `N` sectors' `C_star` pieces.
    pub central_vertex: VertexHandle,
}

/// Split an N-gonal blend face into `4N` quadrilateral sub-faces in
/// the KP13 bi-5 macro-patch 2 × 2 layout. See module header for
/// the algorithm and piece labelling.
pub(crate) fn poke_n_gon_face_into_4n_quads(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
) -> Result<HexPoke4nResult, CoreError> {
    // N-way poke.
    let n_result = poke_n_gon_face_into_n_quads(brep, face, central_point)?;
    poke_n_gon_face_into_4n_quads_from_n_poke(brep, n_result)
}

pub(crate) fn poke_n_gon_face_into_expected_4n_quads(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    expected_n: usize,
) -> Result<HexPoke4nResult, CoreError> {
    let n_result = super::poke_hex::poke_n_gon_face_into_expected_n_quads(
        brep,
        face,
        central_point,
        expected_n,
    )?;
    poke_n_gon_face_into_4n_quads_from_n_poke(brep, n_result)
}

fn poke_n_gon_face_into_4n_quads_from_n_poke(
    brep: &mut BRep,
    n_result: super::poke_hex::HexPokeResult,
) -> Result<HexPoke4nResult, CoreError> {
    let n = n_result.sub_faces.len();
    let central_vertex = n_result.central_vertex;

    // Split all N radial edges (shared between sectors).
    //
    // Edges are collected BEFORE any split, since split_edge creates
    // new edges + half-edges that would confuse the vertex-fan walk.
    split_radial_edges_for_4n_poke(brep, central_vertex, n)?;

    // Per-sector: split 2 outer edges, then subdivide into 4.
    let mut sub_faces: Vec<FaceHandle> = Vec::with_capacity(4 * n);
    for sector_face in n_result.sub_faces {
        let four = subdivide_sector_into_four(brep, sector_face, central_vertex)?;
        sub_faces.extend_from_slice(&four);
    }

    Ok(HexPoke4nResult {
        sub_faces,
        central_vertex,
    })
}

fn split_radial_edges_for_4n_poke(
    brep: &mut BRep,
    central_vertex: VertexHandle,
    n: usize,
) -> Result<(), CoreError> {
    let radial_edges = collect_vertex_outgoing_edges(brep, central_vertex)?;
    if radial_edges.len() != n {
        return Err(CoreError::PreconditionFailed(format!(
            "poke_4n: central vertex has {} outgoing edges, expected {n}",
            radial_edges.len()
        )));
    }
    for edge in radial_edges {
        brep.split_edge(edge, 0.5)?;
    }
    Ok(())
}

/// Walk the vertex fan of `vertex` and return the `EdgeHandle`s of
/// all incident edges, in fan order (`he.twin.next` traversal).
fn collect_vertex_outgoing_edges(
    brep: &BRep,
    vertex: VertexHandle,
) -> Result<Vec<EdgeHandle>, CoreError> {
    let start = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| CoreError::InvalidHandle("poke_4n: vertex handle invalid".into()))?
        .outgoing;
    if !brep.half_edges.contains_key(start) {
        return Err(CoreError::PreconditionFailed(
            "poke_4n: central vertex has no outgoing half-edge".into(),
        ));
    }
    let mut edges = Vec::new();
    let bound = brep.half_edges.len();
    let mut cursor = start;
    for _ in 0..bound {
        let he = &brep.half_edges[cursor];
        edges.push(he.edge);
        let twin = he.twin;
        let twin_data = brep
            .half_edges
            .get(twin)
            .ok_or_else(|| CoreError::InvalidHandle("poke_4n: twin half-edge invalid".into()))?;
        cursor = twin_data.next;
        if cursor == start {
            return Ok(edges);
        }
    }
    Err(CoreError::PreconditionFailed(
        "poke_4n: vertex fan did not close within bound".into(),
    ))
}

/// Walk the outer loop of `face` starting from the half-edge whose
/// origin is `start_vertex`. Returns the half-edge handles in CCW
/// order of traversal (length = number of edges in the loop).
fn walk_loop_from_vertex(
    brep: &BRep,
    face: FaceHandle,
    start_vertex: VertexHandle,
) -> Result<Vec<HalfEdgeHandle>, CoreError> {
    let outer_loop = brep
        .faces
        .get(face)
        .ok_or_else(|| CoreError::InvalidHandle("poke_4n: face invalid".into()))?
        .outer_loop;
    let first_he = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| CoreError::InvalidHandle("poke_4n: loop invalid".into()))?
        .half_edge;

    // Search for the HE whose origin is start_vertex.
    let bound = brep.half_edges.len();
    let mut cursor = first_he;
    let mut found: Option<HalfEdgeHandle> = None;
    for _ in 0..bound {
        if brep.half_edges[cursor].origin == start_vertex {
            found = Some(cursor);
            break;
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == first_he {
            break;
        }
    }
    let start = found.ok_or_else(|| {
        CoreError::PreconditionFailed("poke_4n: start vertex not found on sector loop".into())
    })?;

    // Walk one full traversal from the found HE.
    let mut hes = Vec::new();
    let mut cursor = start;
    for _ in 0..bound {
        hes.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return Ok(hes);
        }
    }
    Err(CoreError::PreconditionFailed(
        "poke_4n: sector loop did not close within bound".into(),
    ))
}

/// Subdivide a sector quad into 4 sub-quads.
///
/// Preconditions:
/// * `face` is a sector quad returned by
///   [`poke_n_gon_face_into_n_quads`].
/// * The 2 radial edges of `face` (the ones incident to
///   `central_vertex`) have already been split at `t = 0.5`, so the
///   sector's outer loop now has 6 half-edges.
/// * The 2 outer edges of `face` are still intact at call time.
///
/// Returns the 4 sub-face handles in order `[O_L, O_C, O_R, C_star]`.
fn subdivide_sector_into_four(
    brep: &mut BRep,
    face: FaceHandle,
    central_vertex: VertexHandle,
) -> Result<[FaceHandle; 4], CoreError> {
    // Walk the 6-HE loop (after phase-2 radial splits).
    let hes6 = walk_loop_from_vertex(brep, face, central_vertex)?;
    if hes6.len() != 6 {
        return Err(CoreError::PreconditionFailed(format!(
            "poke_4n: sector expected 6-HE loop after radial splits, got {}",
            hes6.len()
        )));
    }

    // Positions 2 and 3 are the two outer edges (between the two
    // radial midpoints). Capture edges BEFORE split_edge invalidates
    // adjacent HE handles.
    let outer_lo_edge = brep.half_edges[hes6[2]].edge;
    let outer_hi_edge = brep.half_edges[hes6[3]].edge;
    brep.split_edge(outer_lo_edge, 0.5)?;
    brep.split_edge(outer_hi_edge, 0.5)?;

    // Re-walk: now 8 HEs.
    let hes = walk_loop_from_vertex(brep, face, central_vertex)?;
    if hes.len() != 8 {
        return Err(CoreError::PreconditionFailed(format!(
            "poke_4n: sector expected 8-HE loop after outer splits, got {}",
            hes.len()
        )));
    }

    // Corner positions for centroid (even indices hold the 4
    // original corners: C, V_{2k-2}, V_{2k-1}, V_{2k}).
    let p_c = brep.vertices[brep.half_edges[hes[0]].origin].position;
    let p_v0 = brep.vertices[brep.half_edges[hes[2]].origin].position;
    let p_v1 = brep.vertices[brep.half_edges[hes[4]].origin].position;
    let p_v2 = brep.vertices[brep.half_edges[hes[6]].origin].position;
    let centroid = Point3::new(
        (p_c.x + p_v0.x + p_v1.x + p_v2.x) * 0.25,
        (p_c.y + p_v0.y + p_v1.y + p_v2.y) * 0.25,
        (p_c.z + p_v0.z + p_v1.z + p_v2.z) * 0.25,
    );

    // Euler ops per module-header derivation:
    //   mev(hes[1], centroid) → X
    //   mef(mev.half_edge_twin, hes[5])  → F_new1 (= O_C after mef3)
    //   mef(mef1.half_edge,  hes[7])     → F_new2 (= O_R)
    //   mef(mev.half_edge_twin, hes[3])  → F_new3 (= O_L)
    // Remaining face `face` becomes `C_star`.
    let mev_res = brep.mev(hes[1], centroid)?;
    let mef1 = brep.mef(mev_res.half_edge_twin, hes[5])?;
    let mef2 = brep.mef(mef1.half_edge, hes[7])?;
    let mef3 = brep.mef(mev_res.half_edge_twin, hes[3])?;

    Ok([mef3.face, mef1.face, mef2.face, face])
}

#[cfg(test)]
#[path = "poke_4n_tests.rs"]
mod tests;
