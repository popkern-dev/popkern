//! Promotion: `BlendBoundary` → PnS cap surfaces ready for BRep
//! insertion.
//!
//! End-to-end wrapper over the PnS pipeline:
//!
//! ```text
//! BlendBoundary
//!   → Bi2TensorBorderRing   (pns_adapter)
//!   → pole                  (pns_center eq (5))
//!   → assemble_pns_bi4_cap  (pns_cap_assemble)
//!   → Vec<NurbsSurface>     (pns_bi4_to_nurbs per patch)
//!   → Vec<SurfaceGeometry>  (wrapped)
//! ```
//!
//! The output is `N_fillets · 2 = N_seg` surfaces — one per
//! `BlendBoundary` segment, matching the natural sub-face
//! decomposition of the blend face.
//!
//! Downstream topology insertion (e.g.,
//! `poke_n_gon_face_into_n_quads` + per-sub-face surface
//! assignment) is the caller's responsibility — this module only
//! produces the geometries.

use popkern_core::error::CoreError;
use popkern_core::geometry::gregory::{
    BetaReparam, CapAssembleInput, assemble_pns_bi4_cap, pns_bi4_to_nurbs, pole_from_boundary,
    sigma_center_jet_default,
};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::boundary::BlendBoundary;
use super::pns_adapter::{PnsAdapterError, blend_boundary_to_pns_ring_with_surfaces};
use crate::construct::poke_hex::poke_n_gon_face_into_n_quads;
use crate::fillet_surface::FilletSurfaceResult;

/// Error cases for the PnS promote pipeline.
#[derive(Debug, thiserror::Error)]
pub enum PnsPromoteError {
    /// Adapter step failed (segment count / surface lookup).
    #[error("PnS promote: adapter failed: {0}")]
    Adapter(#[from] PnsAdapterError),
    /// Cap assembly step failed (bad valence, singular F_κ, etc).
    #[error("PnS promote: cap assembly failed: {0}")]
    Assembly(#[from] popkern_core::geometry::gregory::CapAssembleError),
    /// The σ center jet could not be built for this valence.
    #[error("PnS promote: σ center jet failed: {0}")]
    SigmaJet(#[from] popkern_core::geometry::gregory::SigmaJetError),
    /// BRep topology op (poke) failed.
    #[error("PnS promote: BRep operation failed: {0}")]
    Brep(#[from] CoreError),
    /// Mismatched sub-face / surface count after poke.
    #[error(
        "PnS promote: poke produced {n_sub_faces} sub-faces, expected {expected} (adapter surfaces)"
    )]
    SubFaceCountMismatch {
        /// Number of sub-faces after poke.
        n_sub_faces: usize,
        /// Expected count (= number of adapter-produced surfaces).
        expected: usize,
    },
}

/// Result of a full PnS insertion on a blend face.
#[derive(Debug, Clone)]
pub struct PnsInsertResult {
    /// `N` sub-faces after poke, each assigned a PnS NURBS surface.
    pub sub_faces: Vec<FaceHandle>,
    /// Central pole vertex (all N sub-faces meet here).
    pub central_vertex: VertexHandle,
}

/// F_κ order for the cap assembly. 3 is a good default for caps
/// in the supported valence range (n ∈ [5, 10]); could be made
/// configurable if needed.
const DEFAULT_KAPPA: usize = 3;

/// Lowest corner valence the PnS cap is dispatched for.
///
/// Held at 5 for most of this cap's life, for two reasons that have
/// both since been retired:
///
/// - `n = 4` was believed structurally impossible (σ's `1/cos(2π/n)`
///   pole). It is not — the pole is removable and cancels exactly
///   under the paper's Appendix B scalings. See
///   `popkern_core::geometry::gregory::pns_sigma`.
/// - `n = 3` was simply never wired up, even though the adapter has
///   always supported it.
///
/// Both are now measured through the same live-seam gates as the
/// higher valences (`pns_low_valence_gates`), so the floor sits at
/// the construction's true lower bound: a triangle.
///
/// This matters disproportionately because valence tracks corner
/// frequency in reverse — `n = 3` is the cube corner and `n = 4` the
/// pyramid apex, i.e. the two commonest corners in mechanical parts,
/// while `n ≥ 5` is comparatively rare.
pub const PNS_MIN_VALENCE: usize = 3;

/// Produce one `SurfaceGeometry::Nurbs` per `BlendBoundary`
/// segment, ready for assignment to a corresponding sub-face.
///
/// See module docs for the pipeline.
pub fn promote_blend_boundary_to_pns_surfaces(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
) -> Result<Vec<SurfaceGeometry>, PnsPromoteError> {
    let ring = blend_boundary_to_pns_ring_with_surfaces(boundary, fillet_surfaces, brep)?;
    let n = ring.n();
    let pole = pole_from_boundary(&ring);
    let sigma = sigma_center_jet_default(n)?;
    let beta = BetaReparam::bi4(n);
    let input = CapAssembleInput {
        boundary: ring,
        beta,
        sigma,
        pole,
        kappa: DEFAULT_KAPPA,
    };
    let cap = assemble_pns_bi4_cap(&input)?;

    let surfaces: Vec<SurfaceGeometry> = cap
        .patches
        .iter()
        .map(|patch| SurfaceGeometry::Nurbs(pns_bi4_to_nurbs(patch)))
        .collect();
    Ok(surfaces)
}

/// End-to-end PnS insertion: promote a blend boundary to N PnS
/// cap surfaces, poke the `blend_face` into N sub-faces meeting
/// at the central pole, and assign each sub-face its corresponding
/// PnS NURBS surface.
///
/// Returns the N sub-face handles + the central vertex.
///
/// The `sub_faces[k]` entry covers `blend_boundary.segments[2k]`
/// (FilletEnd arc) + `segments[2k+1]` (FaceEdge bridge) — i.e.,
/// the outer boundary of sub-face `k` is the concatenation of arc
/// `k` and bridge `k`.
pub fn insert_pns_cap_into_blend_face(
    brep: &mut BRep,
    blend_face: FaceHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<PnsInsertResult, PnsPromoteError> {
    let surfaces = promote_blend_boundary_to_pns_surfaces(boundary, fillet_surfaces, brep)?;

    // Derive the pole the same way `promote_...` did, so the
    // topology's central vertex sits at the cap's geometric pole.
    let ring = blend_boundary_to_pns_ring_with_surfaces(boundary, fillet_surfaces, brep)?;
    let pole = pole_from_boundary(&ring);

    // The poke op expects either 2N half-edges (even, N ≥ 3) OR N
    // half-edges (N ∈ [5, 16]). Cube corners at N=3, 4 where the
    // bridge edges were consumed by a prior kef end up with N HEs
    // outside that range. Each such HE carries a curve that spans
    // (arc_k + bridge_k); we split it at the arc/bridge junction so
    // the resulting 2N HEs correspond to BlendBoundary segments
    // 1-to-1 (poke then pairs them into sectors aligned with the
    // cap's bi-2 strip midpoint).
    normalize_low_valence_face_for_poke(brep, blend_face, surfaces.len(), boundary)?;

    let poke = poke_n_gon_face_into_n_quads(brep, blend_face, pole)?;
    if poke.sub_faces.len() != surfaces.len() {
        return Err(PnsPromoteError::SubFaceCountMismatch {
            n_sub_faces: poke.sub_faces.len(),
            expected: surfaces.len(),
        });
    }

    for (sub_face, geometry) in poke.sub_faces.iter().zip(surfaces) {
        brep.faces[*sub_face].surface = Some(SurfaceOnFace {
            geometry,
            orientation: true,
        });
    }

    // The sub-faces inherited their boundary pcurves from the pre-poke
    // blend face, whose surface had an entirely different
    // parameterization. Against the bi-4 patches just assigned, those
    // pcurves reach u = 2.70 and v = 1.64 — far outside the patch
    // domain of [0, 1]^2, where a Bezier extrapolates violently
    // (measured 19.287 of pcurve deviation on a unit-scale model).
    // Re-derive them in the new charts.
    for sub_face in &poke.sub_faces {
        let _ = brep.recompute_face_pcurves(*sub_face);
    }

    Ok(PnsInsertResult {
        sub_faces: poke.sub_faces,
        central_vertex: poke.central_vertex,
    })
}

/// If the face has exactly `n` half-edges and `n ∈ {3, 4}` (too
/// few for the poke's normal `N ∈ [5, 16]` split path), split each
/// HE at the arc/bridge junction (from `boundary`) to reach `2n`
/// half-edges. Splitting at the junction — rather than a naive
/// `t = 0.5` — keeps the sub-face loop aligned with the cap's
/// bi-2 strip midpoint (strip[1] = arc/bridge junction point),
/// so the cap's `uv = (0.5, 1)` boundary midpoint matches the
/// sub-face's mid-split vertex position.
fn normalize_low_valence_face_for_poke(
    brep: &mut BRep,
    face: FaceHandle,
    expected_n: usize,
    boundary: &BlendBoundary,
) -> Result<(), CoreError> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::topology::EdgeHandle;
    let outer_loop = brep
        .faces
        .get(face)
        .ok_or_else(|| CoreError::InvalidHandle("face not in BRep".into()))?
        .outer_loop;
    let loop_start = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| CoreError::InvalidHandle("loop not in BRep".into()))?
        .half_edge;
    let mut count = 0usize;
    let mut he_edges: Vec<EdgeHandle> = Vec::new();
    let mut cursor = loop_start;
    loop {
        count += 1;
        he_edges.push(brep.half_edges[cursor].edge);
        cursor = brep.half_edges[cursor].next;
        if cursor == loop_start {
            break;
        }
        if count > 32 {
            return Err(CoreError::PreconditionFailed(
                "face loop exceeds 32 half-edges — likely unterminated".into(),
            ));
        }
    }
    if count != expected_n || !(3..=4).contains(&expected_n) {
        return Ok(());
    }
    let junctions: Vec<popkern_core::geometry::Point3> = (0..expected_n)
        .map(|k| {
            let arc = match &boundary.segments[2 * k] {
                super::boundary::BoundarySegment::FilletEnd { curve, .. } => curve,
                super::boundary::BoundarySegment::FaceEdge { curve, .. } => curve,
            };
            let (_, t_max) = arc.parameter_range();
            arc.point_at(t_max)
        })
        .collect();
    for edge in he_edges {
        let curve_on_edge = brep
            .edges
            .get(edge)
            .and_then(|e| e.curve.as_ref())
            .ok_or_else(|| CoreError::InvalidHandle("edge has no curve".into()))?
            .clone();
        let he_curve = &curve_on_edge.geometry;
        let (t_min, t_max) = he_curve.parameter_range();
        let mut best: Option<(f64, f64)> = None;
        for jp in &junctions {
            if let Some(t) = he_curve.closest_point(jp) {
                if t <= t_min + 1e-9 || t >= t_max - 1e-9 {
                    continue;
                }
                let d = he_curve.point_at(t).distance_to(jp);
                if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                    best = Some((t, d));
                }
            }
        }
        let t_split = best.map(|(t, _)| t).unwrap_or(0.5 * (t_min + t_max));
        brep.split_edge(edge, t_split)?;
    }
    Ok(())
}

/// PnS dispatch is on by default. The `POPKERN_PNS` env var is an
/// *off* switch — set to `0` to disable PnS and force every corner
/// onto the legacy rungs (rational-Gregory ladder, polar-NURBS
/// fallback, etc.). Useful for A/B comparisons or for dodging a
/// suspected PnS regression in the field without a recompile.
pub fn pns_env_enabled() -> bool {
    std::env::var("POPKERN_PNS")
        .map(|v| v != "0")
        .unwrap_or(true)
}

/// Dispatch entry point: try PnS promotion on a blend face, fall
/// through to the caller's next rung on failure. Returns
/// `Some(corner_patches)` on success; `None` if any step rejected
/// the input (invalid valence, segment count, assembly failure,
/// etc.).
///
/// This is gated by `pns_env_enabled()` at the caller's discretion.
/// When invoked and PnS succeeds, it bypasses the three-quad
/// decomposition / K-P promotion ladder entirely.
pub fn try_pns_dispatch(
    brep: &mut BRep,
    blend_face: FaceHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    vertex: VertexHandle,
    classification: crate::vertex_classifier::VertexBlendClassification,
    bridge_fallbacks: Vec<&'static str>,
) -> Result<Vec<crate::types::CornerPatchInfo>, PnsPromoteError> {
    use crate::types::{ContinuityLevel, CornerMethod, CornerPatchInfo};

    // This step used to swallow its failure reason: it returned
    // `Option` and printed the cause behind a `POPKERN_PNS_TRACE`
    // env var, so in any normal run a PnS decline was invisible and
    // the corner silently fell through to the next rung. The reason
    // is now typed and returned, and the caller turns it into a
    // structured `pns.dispatch.declined` diagnostic — which is how
    // the n=3 / n=4 dispatch gaps were finally located.
    let result = insert_pns_cap_into_blend_face(brep, blend_face, boundary, fillet_surfaces)?;

    let n_sides = result.sub_faces.len();
    let patches = result
        .sub_faces
        .into_iter()
        .map(|face| CornerPatchInfo {
            vertex,
            n_sides,
            method: CornerMethod::PolyhedralNetSpline,
            achieved_continuity: ContinuityLevel::G1,
            classification,
            face,
            bridge_fallbacks: bridge_fallbacks.clone(),
            overlap_remediated: false,
            relaxed_continuity_coplanarity: None,
            cap_sphere: None,
            sub_patch_provenance: None,
            setback_cap: false,
            exact_boundary_cap: false,
        })
        .collect();
    Ok(patches)
}

/// Whether N=3 corners are released from direct-apex promotion so the
/// PnS cap can bid for them first.
///
/// Default **off**. This moves a whole valence class onto a different
/// construction route, which is the largest behavioural change in the
/// corner pipeline, so it stays opt-in until the corpus, the polar
/// gates and the direct-apex matrix have all been measured against it.
/// `POPKERN_PNS_FIRST_N3=1` enables it.
pub fn pns_first_at_low_valence() -> bool {
    std::env::var("POPKERN_PNS_FIRST_N3")
        .map(|v| v != "0")
        .unwrap_or(false)
}
