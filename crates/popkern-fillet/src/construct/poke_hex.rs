//! §5.3 Round 8 C4 — topology surgery to split a hexagonal blend
//! face into three quadrilateral sub-faces meeting at a new
//! central vertex.
//!
//! The fillet pipeline's `merge_corner_triangles` pass produces a
//! single hexagonal face at every multi-fillet vertex (3 fillet-
//! end arcs + 3 face-edge bridges). For mixed-radius corners we
//! want the three-quad decomposition (`build_three_quad_decomposition`)
//! geometry to land as three BRep faces — one quad per sub-patch,
//! sharing a new central vertex `C` and three spoke edges.
//!
//! ```text
//!           V1 --- arc0 --- V0
//!            \             /
//!      bridge0           bridge2
//!            \           /
//!             V2 - C - V5
//!            /     |     \
//!        arc1    spokes    arc2
//!          /       |         \
//!        V3 - bridge1 - V4
//! ```
//!
//! Algorithm (3 Euler operators):
//! 1. `mev(he[0], C)` — create central vertex with tab edge `V0↔C`.
//! 2. `mef(he_C_to_V0, he[2])` — split off Face B (V2, C, V0, V1).
//! 3. `mef(he_C_to_V2_in_A, he[4])` — split remaining hex into
//!    Face A (V0, C, V4, V5) and Face C (V2, V3, V4, C).
//!
//! Euler invariant: adds 1 vertex, 3 edges, 2 faces →
//! `ΔV − ΔE + ΔF = 0`, preserving the shell's Euler characteristic.
//!
//! Face-handle assignment (returned in the order that matches the
//! `ThreeQuadDecomposition::sub_patches` array):
//! * `face[0]` = Face B — covers arc 0 + bridge 0 (quad_α)
//! * `face[1]` = Face C — covers arc 1 + bridge 1 (quad_β)
//! * `face[2]` = Face A — covers arc 2 + bridge 2 (quad_γ)

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

/// Result of [`poke_hex_face_into_three_quads`]: the 3 sub-face
/// handles in CCW order matching the input hex's arc/bridge pairs,
/// plus the newly-created central vertex.
#[allow(dead_code)] // wired in R8 C5
#[derive(Debug)]
pub(crate) struct HexPokeResult {
    /// The three sub-face handles in the same index order as
    /// `ThreeQuadDecomposition::sub_patches`:
    /// * `[0]` covers arc 0 + bridge 0
    /// * `[1]` covers arc 1 + bridge 1
    /// * `[2]` covers arc 2 + bridge 2
    pub sub_faces: [FaceHandle; 3],
    /// The newly-created central vertex where all three sub-faces
    /// meet.
    pub central_vertex: VertexHandle,
}

/// Split a hexagonal blend face into three quadrilateral sub-faces
/// sharing a new central vertex at `central_point`.
///
/// The face must have an outer loop with **exactly 6 half-edges**
/// (the R8 C1 hex-boundary output). The returned sub-faces
/// inherit the input face's surface geometry — the caller is
/// expected to overwrite each sub-face's `surface` with the
/// corresponding `RectangularGregory` from the decomposition.
///
/// # Errors
///
/// Returns [`CoreError::PreconditionFailed`] if the face does not
/// have a 6-half-edge outer loop, or any underlying Euler
/// operator fails (shouldn't happen if the boundary is a
/// well-formed hexagon).
#[allow(dead_code)] // wired in R8 C5
pub(crate) fn poke_hex_face_into_three_quads(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
) -> Result<HexPokeResult, CoreError> {
    // Walk the outer loop to collect 6 half-edges in CCW order.
    let face_data = brep.faces.get(face).ok_or_else(|| {
        CoreError::InvalidHandle("face handle is invalid in poke_hex_face_into_three_quads".into())
    })?;
    let outer_loop = face_data.outer_loop;
    let loop_data = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| CoreError::InvalidHandle("outer loop handle is invalid".into()))?;
    let start_he = loop_data.half_edge;

    let mut loop_hes = Vec::with_capacity(6);
    let mut cursor = start_he;
    for _ in 0..12 {
        loop_hes.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start_he {
            break;
        }
    }
    if loop_hes.len() != 6 || cursor != start_he {
        return Err(CoreError::PreconditionFailed(format!(
            "expected 6-half-edge hex loop, got {}",
            loop_hes.len(),
        )));
    }
    let he0 = loop_hes[0];
    let he2 = loop_hes[2];
    let he4 = loop_hes[4];

    // Euler-op 1: `mev(he0, central_point)` — adds vertex C with a
    // tab edge V0↔C spliced before `he0` in the loop.
    let mev_result = brep.mev(he0, central_point)?;
    let he_c_to_v0 = mev_result.half_edge_twin; // C → V0
    let central_vertex = mev_result.vertex;

    // Euler-op 2: `mef(he_C_to_V0, he2)` — splits Face B off.
    // Face B: covers arc 0 + bridge 0 (quad_α in the decomposition).
    let mef1_result = brep.mef(he_c_to_v0, he2)?;
    let face_b = mef1_result.face;

    // After mef1:
    // * old loop (was Face A's loop) now contains: new_he (C→V2),
    //   he_V2_to_V3, he_V3_to_V4, he_V4_to_V5, he_V5_to_V0, he_out
    //   (V0→C).
    // * new_he is `mef1_result.half_edge` = the C→V2 half-edge in
    //   Face A's loop.
    // he4 (V4→V5 on the original loop) stays in the old loop — mef
    // only moved the half-edges between `new_he_twin` and
    // `he1=he_C_to_V0` to the new loop.
    let he_c_to_v2 = mef1_result.half_edge;

    // Euler-op 3: `mef(he_C_to_V2, he4)` — splits remaining hex
    // into Face A (covers arc 2 + bridge 2, quad_γ) and Face C
    // (covers arc 1 + bridge 1, quad_β).
    //
    // The `new_he` of this call (C→V4) stays in the **old** loop
    // (which was Face A after the first mef). The new face this
    // call creates is Face C — walking around new_he_twin's loop:
    // V4→C → C→V2 → V2→V3 → V3→V4 — covers arc 1 + bridge 1.
    let mef2_result = brep.mef(he_c_to_v2, he4)?;
    let face_c = mef2_result.face;

    // The original `face` handle is now Face A (covers arc 2 +
    // bridge 2, quad_γ). Its loop after the second mef:
    //   new_he (C→V4) → he4 (V4→V5) → he5 (V5→V0) → he_out (V0→C)
    // Four half-edges, as expected.
    let face_a = face;

    Ok(HexPokeResult {
        sub_faces: [face_b, face_c, face_a],
        central_vertex,
    })
}

#[cfg(test)]
#[path = "poke_hex_tests.rs"]
mod tests;
