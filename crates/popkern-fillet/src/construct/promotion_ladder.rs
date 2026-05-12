//! Promotion ladder for the three-quad decomposition.
//!
//! Walks the rung hierarchy:
//! deg5 → biquintic-full → biseptic → two-stripe → KP13 bi-5
//! (env-gated) → KP20 bi-6 macro (env-gated) → bi-6 polar →
//! KP binonic.
//!
//! Extracted from `three_quad_apply.rs` to keep that file under
//! the 500-line structural limit.

use popkern_core::geometry::gregory::GregoryPatch;

use super::three_quad_promotions::{
    promote_decomp_to_biquintic_full_with_spoke_g1, promote_decomp_to_biseptic_with_spoke_g1,
    promote_decomp_to_two_stripe_with_spoke_g1,
};
use crate::blend::gregory_three_quad::NQuadDecomposition;
use crate::fillet_surface::FilletSurfaceResult;

/// Run the full promotion ladder. Returns `true` if the KP13
/// bi-5 4N path committed (caller uses 4N poke downstream),
/// `false` otherwise.
///
/// Phase 6c-B: `promote_binonic` callback takes both `sub_patches`
/// and `fillet_surfaces` so the binonic optimiser can use outer-
/// fillet G1 residuals (cap-vs-fillet tangent-plane match along
/// each cap's outer boundary).
#[allow(clippy::too_many_arguments)]
pub(super) fn run_promotion_ladder(
    decomp: &mut NQuadDecomposition,
    fillet_surfaces: &[FilletSurfaceResult],
    relaxed_coplanarity: &mut Option<f64>,
    use_deg5: bool,
    use_biquintic_full: bool,
    use_biseptic: bool,
    use_two_stripe: bool,
    use_kp_binonic: bool,
    use_kp_polar: bool,
    promote_binonic: impl FnOnce(&mut [GregoryPatch], &[FilletSurfaceResult]) -> Option<f64>,
) -> bool {
    if !(use_biquintic_full && use_deg5) {
        return false;
    }
    promote_decomp_to_biquintic_full_with_spoke_g1(&mut decomp.sub_patches);
    if !(use_biseptic
        && promote_decomp_to_biseptic_with_spoke_g1(&mut decomp.sub_patches)
        && use_two_stripe)
    {
        return false;
    }
    promote_decomp_to_two_stripe_with_spoke_g1(&mut decomp.sub_patches);
    if use_kp_polar
        && decomp.sub_patches.len() >= 3
        && let Some(patches_4n) =
            super::kp13_bi5_promote::try_kp13_bi5_dispatch(&decomp.sub_patches, fillet_surfaces)
    {
        decomp.sub_patches = patches_4n;
        return true;
    }
    if use_kp_polar
        && decomp.sub_patches.len() >= 5
        && let Some(patches_4n) = super::kp20_polar_dispatch::try_kp20_polar_dispatch(
            &decomp.sub_patches,
            fillet_surfaces,
        )
    {
        decomp.sub_patches = patches_4n;
        return true;
    }
    let polar_committed = use_kp_polar
        && decomp.sub_patches.len() >= 5
        && super::kp_polar_bi6_promote::promote_and_update_relaxed(
            &mut decomp.sub_patches,
            fillet_surfaces,
            relaxed_coplanarity,
        );
    if use_kp_binonic && !polar_committed {
        *relaxed_coplanarity = promote_binonic(&mut decomp.sub_patches, fillet_surfaces);
    }
    false
}
