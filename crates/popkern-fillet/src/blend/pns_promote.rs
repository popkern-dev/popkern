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

use popkern_core::geometry::gregory::{
    BetaReparam, CapAssembleInput, assemble_pns_bi4_cap, pns_bi4_to_nurbs, pole_from_boundary,
    sigma_center_jet_default,
};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::BRep;

use super::boundary::BlendBoundary;
use super::pns_adapter::{PnsAdapterError, blend_boundary_to_pns_ring_with_surfaces};
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
