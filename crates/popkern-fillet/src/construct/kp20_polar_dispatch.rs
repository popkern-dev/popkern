//! Gated KP20 polar macro-cap candidate construction for fillet dispatch.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6CapInput, RectangularGregoryBiNonic, TensorBorderBi3,
    assemble_kp20_bi6_cap, elevate_bi6_grid_to_bi9, measure_kp20_bi6_macro_cap,
};

use super::kp_polar_bi6_promote::SELF_INT_TOL;

pub(super) fn try_kp20_polar_macro_candidate(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Vec<GregoryPatch>, String> {
    let cap = assemble_kp20_bi6_cap(&Kp20Bi6CapInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let report = measure_kp20_bi6_macro_cap(&cap, tensor_border).map_err(|e| e.to_string())?;
    if !report.passes(1e-8) {
        return Err(format!("validation residuals exceed tolerance: {report:?}"));
    }

    Ok(cap
        .patches_bi6
        .into_iter()
        .map(elevate_bi6_grid_to_bi9)
        .map(|grid| GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid })))
        .collect())
}

pub(super) fn first_self_intersecting_patch(candidate: &[GregoryPatch]) -> Option<usize> {
    candidate
        .iter()
        .position(|gp| gp.self_intersects_with_tolerance(SELF_INT_TOL))
}
