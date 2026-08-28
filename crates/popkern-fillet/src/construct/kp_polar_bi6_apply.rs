//! `popkern-fillet` entry point for the K-P 2015
//! bi-6 polar construction (N ≥ 5).
//!
//! This module is a thin wrapper over the `popkern-core` polar
//! pipeline (I.1–I.10d). Callers that already have a
//! [`PolarBoundaryInput`] can invoke [`try_kp_polar_bi6_cap`] to
//! get the N bi-9 CP grids for the polar cap.
//!
//! # Scope
//!
//! This commit lands the **flag** ([`FilletOptions::use_kp_polar_gregory`]
//! in `types.rs`) and the **wrapper entry point**. It does NOT
//! yet include the automatic dispatch from [`BlendBoundary`] at
//! [`insert_vertex_blend`](super::blend_insert::insert_vertex_blend) —
//! that requires a `BlendBoundary` → `PolarBoundaryInput`
//! conversion (planned as I.10f).
//!
//! # Pipeline
//!
//! The wrapper composes the I.10 core-level steps:
//!
//! 1. [`extract_bi3_tensor_border_from_polar_input`] — derives a
//!    bi-3 [`TensorBorder`] from the geometric inputs.
//! 2. [`convert_polar_to_bi6_cap_input`] — packs the data into
//!    [`Bi6CapInput`] (tensor-border + cap-net).
//! 3. [`f5_solve_optimal_bi6_cap`] — F_5-minimises over the free
//!    parameter set, returns both optimal parameters and
//!    fully-assembled patches.
//!
//! Returns [`Err`] on any failure (validation, singular systems,
//! etc.) so the caller's fallback ladder can fire.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    PolarBoundaryInput, convert_polar_to_bi6_cap_input, extract_bi3_tensor_border_from_polar_input,
    f_kappa_solve_optimal_bi6_cap, f5_solve_optimal_bi6_cap,
};

/// Run the K-P 2015 bi-6 polar construction end-to-end on a
/// [`PolarBoundaryInput`]. Returns the `N` patch grids already
/// degree-elevated to bi-9 (10×10 CPs each), ready to wrap as
/// [`RectangularGregoryBiNonic`](popkern_core::geometry::gregory::RectangularGregoryBiNonic).
///
/// Caller is responsible for wrapping the returned grids into
/// `GregoryPatch::RectangularBiNonic` variants for insertion into
/// a [`BRep`](popkern_core::topology::BRep).
///
/// # Errors
///
/// Returns an error string if:
/// - The polar input fails validation (non-finite radii etc.).
/// - `n < 3`.
/// - The tensor-border extraction fails.
/// - The circulant systems are singular without proper pinning
///   (e.g., `n ∈ {3, 6}` without `p0_21_pinned` — though the
///   bridge auto-pins at those valences).
/// - The F_5 SPD solve is singular.
// Scaffold: no crate-internal caller yet (auto-dispatch comes in
// I.10f). Tests exercise the function; silence dead-code.
#[allow(dead_code)]
pub fn try_kp_polar_bi6_cap(
    polar_input: &PolarBoundaryInput,
) -> Result<Vec<[[Point3; 10]; 10]>, String> {
    let n = polar_input.n();
    if n < 3 {
        return Err(format!("polar bi-6 cap needs n ≥ 3, got {n}"));
    }
    let tensor_border = extract_bi3_tensor_border_from_polar_input(polar_input)?;
    let cap_input = convert_polar_to_bi6_cap_input(polar_input, &tensor_border)?;
    let (_optimal, out) = f5_solve_optimal_bi6_cap(&cap_input)?;
    Ok(out.patches_bi9)
}

/// Variant of [`try_kp_polar_bi6_cap`] parameterised on the
/// functional order κ ∈ {1, 3, 4, 5}. Same pipeline otherwise
/// (tensor-border → bi-6 cap input → SPD solve → bi-9). Used by
/// the F_κ fallback ladder (f5-tuning-plan.md option C): the
/// caller retries with κ ∈ {4, 3, 1} when κ = 5 produces self-
/// intersecting output.
#[allow(dead_code)]
pub fn try_kp_polar_bi6_cap_with_kappa(
    polar_input: &PolarBoundaryInput,
    kappa: usize,
) -> Result<Vec<[[Point3; 10]; 10]>, String> {
    let n = polar_input.n();
    if n < 3 {
        return Err(format!("polar bi-6 cap needs n ≥ 3, got {n}"));
    }
    let tensor_border = extract_bi3_tensor_border_from_polar_input(polar_input)?;
    let cap_input = convert_polar_to_bi6_cap_input(polar_input, &tensor_border)?;
    let (_optimal, out) = f_kappa_solve_optimal_bi6_cap(&cap_input, kappa)?;
    Ok(out.patches_bi9)
}

#[cfg(test)]
#[path = "kp_polar_bi6_apply_tests.rs"]
mod tests;
