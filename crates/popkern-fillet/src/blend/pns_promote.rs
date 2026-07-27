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
//!   → Vec<SurfaceGeometry::PolyhedralSpline> (net + provenance per patch)
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
    BetaReparam, CapAssembleInput, assemble_pns_bi4_cap, pole_from_boundary,
    sigma_center_jet_default,
};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::boundary::BlendBoundary;
use super::pns_adapter::{PnsAdapterError, blend_boundary_to_pns_ring_with_surfaces};
pub(crate) use super::pns_cap_authoring::reauthor_pns_cap_edge_geometry;
use super::pns_cap_authoring::{author_pns_cap_edge_geometry, measure_sub_face_patch_assignment};
use super::pns_rim_materialize::{
    align_boundary_to_loop, materialize_pns_rim_from_boundary, normalize_low_valence_face_for_poke,
};
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

    // First-class provenance (PK-200): each sector CAN ship as
    // PolyhedralSpline carrying {valence, sector}, so unfillet /
    // modify / STEP can still say WHICH piece of WHICH cap a face is.
    // TEMPORARY CONTAINMENT: default-off behind
    // POPKERN_PNS_PROVENANCE=1 — direct evaluation of the variant
    // perturbs tessellation-derived measures by ~1.1% against the
    // Gregory parity bars (rational_gregory_measure_parity), so the
    // variant's refinement/measure integration must land before this
    // ships by default. The chart-transpose fix and the emit
    // machinery stay exercised by the promote unit tests either way.
    let provenance_enabled = std::env::var("POPKERN_PNS_PROVENANCE")
        .map(|v| v == "1")
        .unwrap_or(false);
    if !provenance_enabled {
        let surfaces: Vec<SurfaceGeometry> = cap
            .patches
            .iter()
            .map(|patch| {
                SurfaceGeometry::Nurbs(popkern_core::geometry::gregory::pns_bi4_to_nurbs(patch))
            })
            .collect();
        return Ok(surfaces);
    }
    let surfaces: Vec<SurfaceGeometry> = cap
        .patches
        .iter()
        .enumerate()
        .map(|(sector, patch)| {
            // The assembler stores `patch[u_c][v_c]`; the variant's
            // net is `net[j][i]` with `i` the u index — transpose or
            // the shipped chart swaps u and v (caught by the promote
            // unit tests' chart-convention assertions).
            let net: [[popkern_core::geometry::Point3; 5]; 5] =
                core::array::from_fn(|j| core::array::from_fn(|i| patch[i][j]));
            SurfaceGeometry::PolyhedralSpline(Box::new(
                popkern_core::geometry::polyhedral_spline_patch::PolyhedralSplinePatch::new(
                    net,
                    popkern_core::geometry::polyhedral_spline_patch::PolyhedralSplineProvenance {
                        valence: n,
                        sector,
                    },
                ),
            ))
        })
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
    // The poke op expects either 2N half-edges (even, N ≥ 3) OR N
    // half-edges (N ∈ [5, 16]). Cube corners at N=3, 4 where the
    // bridge edges were consumed by a prior kef end up with N HEs
    // outside that range. Each such HE carries a curve that spans
    // (arc_k + bridge_k); we split it at the arc/bridge junction so
    // the resulting 2N HEs correspond to BlendBoundary segments
    // 1-to-1 (poke then pairs them into sectors aligned with the
    // cap's bi-2 strip midpoint).
    let expected_n = boundary.segments.len() / 2;
    normalize_low_valence_face_for_poke(brep, blend_face, expected_n, boundary)?;

    // The committed loop can traverse the rim in the OPPOSITE
    // rotational direction to the boundary ring, and can start at any
    // segment. The poke pairs loop half-edges (2k, 2k+1) into sector
    // k, so a misaligned ring builds every cap sector for a different
    // span than the sub-face it lands on (measured 2026-07-26:
    // half-sector mispairing at n=3, cap_curve_gap 8.3e-2 once the
    // rim itself was exact). Align the boundary to the loop before
    // anything consumes it, then materialize the rim from it — at the
    // low valences the committed loop is a scaffold of stale Lines
    // whose curves and vertices only the serving route rebuilds.
    let aligned = if (3..=4).contains(&expected_n) {
        align_boundary_to_loop(brep, blend_face, boundary)
    } else {
        None
    };
    if let Some(aligned) = aligned.as_ref() {
        materialize_pns_rim_from_boundary(brep, blend_face, aligned);
    }
    let boundary = aligned.as_ref().unwrap_or(boundary);

    let surfaces = promote_blend_boundary_to_pns_surfaces(boundary, fillet_surfaces, brep)?;

    // Derive the pole the same way `promote_...` did, so the
    // topology's central vertex sits at the cap's geometric pole.
    let ring = blend_boundary_to_pns_ring_with_surfaces(boundary, fillet_surfaces, brep)?;
    let pole = pole_from_boundary(&ring);

    let poke = poke_n_gon_face_into_n_quads(brep, blend_face, pole)?;
    if poke.sub_faces.len() != surfaces.len() {
        return Err(PnsPromoteError::SubFaceCountMismatch {
            n_sub_faces: poke.sub_faces.len(),
            expected: surfaces.len(),
        });
    }

    // Assign each sub-face its patch. On the aligned low-valence path
    // the correspondence is verified by measurement (the sector pairing
    // has been observed rotated; the true patch reads ~1e-3 against the
    // face's rim, every other patch at sector scale ~5e-2); the n >= 5
    // route keeps its long-standing positional zip untouched.
    if aligned.is_some() {
        let assignment = measure_sub_face_patch_assignment(brep, &poke.sub_faces, &surfaces);
        for (face_idx, surface_idx) in assignment.into_iter().enumerate() {
            let sub_face = poke.sub_faces[face_idx];
            brep.faces[sub_face].surface = Some(SurfaceOnFace {
                geometry: surfaces[surface_idx].clone(),
                orientation: true,
            });
        }
    } else {
        for (sub_face, geometry) in poke.sub_faces.iter().zip(surfaces) {
            brep.faces[*sub_face].surface = Some(SurfaceOnFace {
                geometry,
                orientation: true,
            });
        }
    }

    // The sub-faces inherited their boundary pcurves from the pre-poke
    // blend face, whose surface had an entirely different
    // parameterization. On the aligned low-valence path the chart
    // correspondence is analytic by construction, and general 2-D
    // projection cannot be trusted to rebuild the pcurves — measured
    // 2026-07-26, it lands interior preimages (cap_uv (0.512, 0.639)
    // for a rim sample) and turns a 5e-3 surface into a 4.4e-2 pcurve
    // reading — so the cap edge geometry is authored from the chart.
    // The n >= 5 route keeps its long-standing projection recompute.
    if aligned.is_some() {
        author_pns_cap_edge_geometry(brep, &poke.sub_faces, poke.central_vertex);
    } else {
        for sub_face in &poke.sub_faces {
            let _ = brep.recompute_face_pcurves(*sub_face);
        }
    }

    Ok(PnsInsertResult {
        sub_faces: poke.sub_faces,
        central_vertex: poke.central_vertex,
    })
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
