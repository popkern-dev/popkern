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
    let sigma = sigma_center_jet_default(n);
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

    Ok(PnsInsertResult {
        sub_faces: poke.sub_faces,
        central_vertex: poke.central_vertex,
    })
}
