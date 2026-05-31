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
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle,
};

/// Result of [`poke_n_gon_face_into_n_quads`]: the `N` sub-face
/// handles in CCW order matching the input polygon's arc/bridge
/// pairs, plus the newly-created central vertex.
#[allow(dead_code)]
#[derive(Debug)]
pub(crate) struct HexPokeResult {
    /// The `N` sub-face handles in arc-index order:
    /// * `sub_faces[k]` covers arc `k` + bridge `k` for `k ∈ 0..N`.
    pub sub_faces: Vec<FaceHandle>,
    /// The newly-created central vertex where all `N` sub-faces meet.
    pub central_vertex: VertexHandle,
}

/// Split an N-gonal blend face into `N` quadrilateral sub-faces
/// sharing a new central vertex at `central_point`.
///
/// Accepts two input shapes on the outer loop:
///
/// 1. **`2N` half-edges** — alternating arc/bridge order, the
///    `R8 C1` N-gon boundary output produced when adjacent-fillet
///    bridges survive the kef/kev Euler ops as topologically-
///    distinct edges. Typical at `N ∈ {3, 4}` on cube / square-
///    pyramid corners.
///
/// 2. **`N` half-edges** — pure fillet-arc loop, no topological
///    bridges. Typical at `N ≥ 5` sharp-apex geometry where the
///    bridge region collapses into the arc edges during kef
///    merging (the `BlendBoundary` still carries 2N geometric
///    segments — arcs + bridge curves on the supporting faces —
///    but the BRep loop itself is N-gonal). Each arc is split at
///    its midpoint parameter (`t=0.5`) to synthesize the missing
///    bridge-endpoint vertices, producing an equivalent `2N`
///    loop before the main algorithm runs.
///
/// The returned sub-faces inherit the input face's surface
/// geometry — the caller is expected to overwrite each sub-face's
/// `surface` with the corresponding `RectangularGregory` from the
/// decomposition.
///
/// Algorithm (`N` Euler operators, after normalising to `2N`):
/// 1. `mev(he[0], C)` — create central vertex with tab edge `V0↔C`.
/// 2. For `i ∈ 1..N`: `mef(he_C_to_V_{2(i-1)}, he[2i])` — split off
///    the `i`-th quad. The `N-1`-th call leaves the final quad
///    (covering arc `N-1` + bridge `N-1`) as the remaining face.
///
/// Euler invariant: adds 1 vertex, `N` edges, `N-1` faces →
/// `ΔV − ΔE + ΔF = 0`, preserving the shell's Euler characteristic.
/// When the N-HE preprocessing runs, it additionally adds
/// `N` vertices and `N` edges via `split_edge` (each split:
/// +1 vertex, +1 edge, 0 faces → also Euler-neutral).
///
/// # Errors
///
/// Returns [`CoreError::PreconditionFailed`] if the face's outer
/// loop doesn't have either `N ∈ [3, 16]` half-edges or `2N` half-
/// edges (even, `N ≥ 3`), or any underlying Euler operator fails.
#[allow(dead_code)]
pub(crate) fn poke_n_gon_face_into_n_quads(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
) -> Result<HexPokeResult, CoreError> {
    poke_n_gon_face_into_n_quads_with_expected(brep, face, central_point, None)
}

pub(crate) fn poke_n_gon_face_into_expected_n_quads(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    expected_n: usize,
) -> Result<HexPokeResult, CoreError> {
    poke_n_gon_face_into_n_quads_with_expected(brep, face, central_point, Some(expected_n))
}

fn poke_n_gon_face_into_n_quads_with_expected(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    expected_n: Option<usize>,
) -> Result<HexPokeResult, CoreError> {
    let outer_loop = brep
        .faces
        .get(face)
        .ok_or_else(|| {
            CoreError::InvalidHandle(
                "face handle is invalid in poke_n_gon_face_into_n_quads".into(),
            )
        })?
        .outer_loop;

    // Initial walk to classify the input shape.
    let mut loop_hes = walk_outer_loop(brep, outer_loop)?;
    let initial_count = loop_hes.len();

    // Classify: 2N (existing path) vs N (normalise via split_edge).
    //
    // Phase 6c-E: N-HE path now accepts `N ≥ 4` (was `≥ 5`).
    // Pyramid_uniform produces 4 outer half-edges (one per
    // slant-face apex-corner) because Phase 2a's
    // `kef_cross_face_v_edges` collapses arc + bridge pairs on
    // each slant face into a single edge.  For N=3 (cube) the
    // bridges still survive so the count is 6; the 2N branch
    // catches that.
    let n = if let Some(expected) = expected_n {
        if initial_count == 2 * expected {
            expected
        } else if initial_count == expected && (4..=16).contains(&expected) {
            normalize_n_he_loop_to_2n(brep, outer_loop, &mut loop_hes, initial_count)?;
            expected
        } else {
            return Err(CoreError::PreconditionFailed(format!(
                "expected {expected} or {} half-edges for requested N, got {initial_count}",
                2 * expected
            )));
        }
    } else if initial_count >= 6 && initial_count.is_multiple_of(2) {
        initial_count / 2
    } else if (4..=16).contains(&initial_count) {
        normalize_n_he_loop_to_2n(brep, outer_loop, &mut loop_hes, initial_count)?;
        initial_count
    } else {
        return Err(CoreError::PreconditionFailed(format!(
            "expected 2N half-edges (even, N ≥ 3) or N half-edges (N ∈ [4, 16]), got {initial_count}"
        )));
    };

    // Euler-op 1: `mev(he[0], central_point)` — adds vertex C with
    // tab edge V_0 ↔ C spliced before `he[0]` in the loop.
    let mev_result = brep.mev(loop_hes[0], central_point)?;
    let central_vertex = mev_result.vertex;
    // The half-edge from C pointing back toward V_0: this is the
    // seed for the first mef.
    let mut prev_he_from_c = mev_result.half_edge_twin;

    // Euler-op 2..N: `mef(he_C_to_V_{2(i-1)}, he[2i])` splits off
    // the i-th sub-face (i = 1..N-1). After each mef, the newly
    // produced `mef.half_edge` is the half-edge `C → V_{2i}` in the
    // remaining loop — seed for the next mef.
    let mut sub_faces: Vec<FaceHandle> = Vec::with_capacity(n);
    for i in 1..n {
        let mef_result = brep.mef(prev_he_from_c, loop_hes[2 * i])?;
        sub_faces.push(mef_result.face);
        prev_he_from_c = mef_result.half_edge;
    }
    // The original `face` handle is now the final sub-face (covers
    // arc `N-1` + bridge `N-1`).
    sub_faces.push(face);

    // Re-order so `sub_faces[k]` covers arc `k` + bridge `k`. The
    // above construction produces:
    //   sub_faces[0] = Face after mef1 (covers arc 0 + bridge 0)
    //   sub_faces[1] = Face after mef2 (covers arc 1 + bridge 1)
    //   ...
    //   sub_faces[N-2] = Face after mef_{N-1} (covers arc N-2 + bridge N-2)
    //   sub_faces[N-1] = original `face` (covers arc N-1 + bridge N-1)
    // — which is already the correct arc-index order.

    Ok(HexPokeResult {
        sub_faces,
        central_vertex,
    })
}

fn normalize_n_he_loop_to_2n(
    brep: &mut BRep,
    outer_loop: LoopHandle,
    loop_hes: &mut Vec<HalfEdgeHandle>,
    initial_count: usize,
) -> Result<(), CoreError> {
    // N-HE input at N ≥ 4. Split each arc at `t=0.5` to reach
    // 2N; the synthetic midpoint vertices become the "corner"
    // vertices where bridge endpoints would have existed in the
    // 2N shape.
    let edges: Vec<EdgeHandle> = loop_hes
        .iter()
        .map(|he| brep.half_edges[*he].edge)
        .collect();
    for edge in edges {
        brep.split_edge(edge, 0.5)?;
    }
    *loop_hes = walk_outer_loop(brep, outer_loop)?;
    if loop_hes.len() != 2 * initial_count {
        return Err(CoreError::PreconditionFailed(format!(
            "N-HE -> 2N split expected {} half-edges, got {}",
            2 * initial_count,
            loop_hes.len()
        )));
    }
    Ok(())
}

/// Back-compat wrapper: N=3 specialization of
/// [`poke_n_gon_face_into_n_quads`].
#[allow(dead_code)]
pub(crate) fn poke_hex_face_into_three_quads(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
) -> Result<HexPokeResult, CoreError> {
    let res = poke_n_gon_face_into_n_quads(brep, face, central_point)?;
    if res.sub_faces.len() != 3 {
        return Err(CoreError::PreconditionFailed(format!(
            "poke_hex: expected 3 sub-faces, got {}",
            res.sub_faces.len()
        )));
    }
    Ok(res)
}

/// Check whether a face's outer loop has a half-edge count
/// compatible with [`poke_n_gon_face_into_n_quads`].
///
/// Accepts `2 * n_sides` (classical arc/bridge alternating loop)
/// at all valences; additionally accepts `n_sides` (pure-arc loop,
/// sharp-apex case — `split_edge` normalises to 2N internally)
/// when `n_sides ≥ 5`, or when `n_sides == 4` and `accept_n4_he`
/// is true (Phase 6c-E pyramid_uniform path under sphere-fit
/// substitution). At `N == 3` only the 2N form is accepted —
/// bridges are expected to survive on cube-like geometries.
#[allow(dead_code)]
pub(crate) fn face_has_2n_outer_half_edges_with(
    brep: &BRep,
    face: FaceHandle,
    n_sides: usize,
    accept_n4_he: bool,
) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return false;
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return false;
    };
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return false;
    }
    let cap = (4 * n_sides + 4).max(12);
    let mut cursor = start;
    for n in 1..=cap {
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            // N-HE acceptance: `n_sides >= 5` always, plus N=4 only
            // when the caller opts in (sphere-fit substitution
            // succeeded — pyramid_uniform's commercial-grade path).
            // Without `accept_n4_he`, N=4 stays on the legacy 2N-only
            // path so pyramid4 validation-sweep fixtures keep their
            // pre-Phase-6c-E behavior.
            let accepts_n_he = (n_sides >= 5 || (n_sides == 4 && accept_n4_he)) && n == n_sides;
            return n == 2 * n_sides || accepts_n_he;
        }
    }
    false
}

/// Backwards-compatible wrapper: legacy callers (validation sweep,
/// non-sphere-fit paths) get the pre-Phase-6c-E gate where N=4 is
/// only accepted in the 2N form.
#[allow(dead_code)]
pub(crate) fn face_has_2n_outer_half_edges(brep: &BRep, face: FaceHandle, n_sides: usize) -> bool {
    face_has_2n_outer_half_edges_with(brep, face, n_sides, false)
}

/// Walk a face's outer loop and collect its half-edges.
/// Capped at 64 entries to guard against malformed loops.
fn walk_outer_loop(
    brep: &BRep,
    outer_loop: popkern_core::topology::LoopHandle,
) -> Result<Vec<HalfEdgeHandle>, CoreError> {
    let loop_data = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| CoreError::InvalidHandle("outer loop handle is invalid".into()))?;
    let start_he = loop_data.half_edge;
    let mut hes = Vec::with_capacity(12);
    let mut cursor = start_he;
    for _ in 0..64 {
        hes.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start_he {
            return Ok(hes);
        }
    }
    Err(CoreError::PreconditionFailed(
        "outer loop did not terminate within 64 steps".into(),
    ))
}

#[cfg(test)]
#[path = "poke_hex_tests.rs"]
mod tests;
