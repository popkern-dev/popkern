//! Gated KP20 polar macro-cap candidate construction for fillet dispatch.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6CapInput, PolarBoundaryInput, RectangularGregoryBiNonic, Side, Spoke,
    SpokeEnd, TensorBorder, TensorBorderBi3, assemble_kp20_bi6_cap, elevate_bi6_grid_to_bi9,
    extract_bi3_tensor_border_from_polar_input, measure_kp20_bi6_macro_cap,
};

use super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::kp_polar_bi6_promote::extract_fillet_radius;
use super::three_quad_apply::extract_direct_boundary_input;
use crate::fillet_surface::FilletSurfaceResult;

pub(super) fn try_kp20_polar_dispatch(
    sub_patches: &[GregoryPatch],
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<Vec<GregoryPatch>> {
    let n = sub_patches.len();
    let mut two_stripe = Vec::with_capacity(n);
    for (idx, patch) in sub_patches.iter().enumerate() {
        let GregoryPatch::RectangularTwoStripe(ts) = patch else {
            trace("non-two-stripe input", n, Some(format!("sub-patch {idx}")));
            return None;
        };
        two_stripe.push((**ts).clone());
    }

    let spokes: Vec<Spoke> = (0..n)
        .map(|k| Spoke {
            a: SpokeEnd {
                patch_idx: k,
                side: Side::UOne,
                reversed: false,
            },
            b: SpokeEnd {
                patch_idx: (k + 1) % n,
                side: Side::VZero,
                reversed: true,
            },
        })
        .collect();
    let Some(mut direct) = extract_direct_boundary_input(&two_stripe, &spokes, None) else {
        trace("extract_direct_boundary_input", n, None);
        return None;
    };
    if fillet_surfaces.len() >= n {
        for (k, spoke) in direct.spokes.iter_mut().enumerate() {
            if let Some(radius) = extract_fillet_radius(&fillet_surfaces[k]) {
                spoke.radius_at_c = radius;
                spoke.radius_at_v = radius;
            }
        }
    }
    let polar = match PolarBoundaryInput::from_direct_boundary_input(&direct) {
        Ok(polar) => polar,
        Err(e) => {
            trace("from_direct_boundary_input", n, Some(e));
            return None;
        }
    };
    let tensor_border = match extract_bi3_tensor_border_from_polar_input(&polar) {
        Ok(TensorBorder::Bi3(tensor_border)) => tensor_border,
        Ok(TensorBorder::Bi5(_)) => {
            trace("tensor-border kind", n, Some("expected bi-3".into()));
            return None;
        }
        Err(e) => {
            trace("tensor-border extract", n, Some(e));
            return None;
        }
    };
    let candidate = match try_kp20_polar_macro_candidate(polar.central_vertex, &tensor_border) {
        Ok(candidate) => candidate,
        Err(e) => {
            trace("candidate", n, Some(e));
            return None;
        }
    };
    if candidate.len() != 4 * n {
        trace(
            "output length mismatch",
            n,
            Some(format!("expected {}, got {}", 4 * n, candidate.len())),
        );
        return None;
    }
    if let Some(idx) = first_self_intersecting_patch(&candidate) {
        trace("self-intersection", n, Some(format!("patch index {idx}")));
        return None;
    }
    trace("accepted", n, Some(format!("4N={} bi-9 patches", 4 * n)));
    Some(candidate)
}

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

fn trace(stage: &str, n: usize, detail: Option<String>) {
    if std::env::var_os("POPKERN_KP20_POLAR_TRACE").is_some() {
        match detail {
            Some(d) => eprintln!("[kp20_polar N={n}] {stage}: {d}"),
            None => eprintln!("[kp20_polar N={n}] {stage}"),
        }
    }
}
