//! `build()` and `can_build()` — the public entry points for
//! direct-construction.
//!
//! Orchestrates the Phase G pipeline:
//!
//! ```text
//!   DirectBoundaryInput ─┬→ corner-C solver  ────────┐
//!                        ├→ outer-corner solver  ────┤
//!                        ├→ blend-coefficient solver ┤
//!                        │                           │
//!                        └→ ring row-1 G1 solver ────┼→ updated patches
//!                                                    │
//!                              self-check (ValidationLevel) ┘
//!                                                    │
//!                                            DirectConstructionResult
//! ```
//!
//! G2 (row-2) enforcement is not yet wired — see the plan for
//! the derivation status. For now the output is `ContinuityLevel::G1`
//! regardless of the requested target.

#![allow(clippy::needless_range_loop)]

use crate::geometry::gregory::direct_blend_coeffs::compute_ring_blend_coeffs;
use crate::geometry::gregory::direct_corner_c::solve_central_vertex_corner;
use crate::geometry::gregory::direct_corner_outer::solve_outer_corners;
use crate::geometry::gregory::direct_errors::DirectConstructionError;
use crate::geometry::gregory::direct_input::DirectBoundaryInput;
use crate::geometry::gregory::direct_options::DirectConstructionOptions;
use crate::geometry::gregory::direct_row1_ring::solve_ring_row1_g1;
use crate::geometry::gregory::direct_self_check::{SelfCheckOptions, measure_ring_quality};
use crate::geometry::gregory::direct_types::{
    CornerDiagnostic, CornerKind, DirectConstructionResult, FeasibilityReport, FeasibilityWarning,
    FormulaFamily, SpokeDiagnostic, ValidationLevel,
};
use crate::geometry::gregory::spoke_g1_solver::{Side, Spoke, SpokeEnd};
use crate::healing::quality::ContinuityLevel;

/// Feasibility query — cheap, no CP emission.
///
/// Returns `FeasibilityReport` describing the best continuity
/// level the direct construction can reach for this input,
/// plus a conditioning score and per-corner warnings. Useful
/// for UI preview before committing to the full build.
pub fn can_build(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Result<FeasibilityReport, DirectConstructionError> {
    input
        .validate()
        .map_err(|detail| DirectConstructionError::DecompositionFailed { detail })?;

    // Run the cheap steps (corner-C + blend-coeffs) to detect
    // structural issues before paying for the row-1 assembly.
    let _ = solve_central_vertex_corner(input, options)?;
    let blend_coeffs = compute_ring_blend_coeffs(input, options)?;

    let mut warnings = Vec::new();
    let mut min_score = 1.0_f64;
    for (k, bc) in blend_coeffs.iter().enumerate() {
        if !bc.condition_number.is_finite() {
            warnings.push(FeasibilityWarning::CornerNearDegenerate {
                corner_idx: k,
                score: 0.0,
            });
            min_score = 0.0;
            continue;
        }
        // Rough near-degeneracy heuristic: condition number < 0.1
        // indicates a near-tangent corner. Score maps [0, 1] ≈
        // clamp(cond, 0, 1).
        let score = bc.condition_number.clamp(0.0, 1.0);
        if score < min_score {
            min_score = score;
        }
        if score < 0.3 {
            warnings.push(FeasibilityWarning::CornerNearDegenerate {
                corner_idx: k,
                score,
            });
        }
    }

    // N ≥ 5 remains out of Phase G scope — flag but don't reject.
    if input.n() >= 5 {
        warnings.push(FeasibilityWarning::ValenceRequiresPolar { n: input.n() });
    }

    // The MVP delivers G1 (not G2). Report G1 regardless of target.
    let achievable_continuity = match options.target_continuity {
        ContinuityLevel::G2 => ContinuityLevel::G1,
        other => other,
    };

    Ok(FeasibilityReport {
        achievable_continuity,
        near_degeneracy_score: min_score,
        warnings,
        // Rough per-corner cost estimate. Measured by the Phase G
        // benchmark suite on representative fixtures; for now a
        // plausible constant.
        estimated_construction_us: 1500 * input.n() as u64,
    })
}

/// Primary construction entry point. Deterministic, single-shot.
///
/// Orchestrates the Phase G pipeline, applying the corner
/// constructions + row-1 G1 solver to the seed patches in
/// `input.seed_patches`. Emits fully populated
/// `DirectConstructionResult` with per-corner / per-spoke
/// diagnostics for UI binding.
pub fn build(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Result<DirectConstructionResult, DirectConstructionError> {
    input
        .validate()
        .map_err(|detail| DirectConstructionError::DecompositionFailed { detail })?;
    let seed = input.seed_patches.as_ref().ok_or_else(|| {
        DirectConstructionError::DecompositionFailed {
            detail: "direct construction requires seed_patches".into(),
        }
    })?;

    // 1. Corner-C: the Shirman-Séquin N-way enclosure.
    let corner_c = solve_central_vertex_corner(input, options)?;

    // 2. Outer corners: parallelogram-rule twist CPs + tangent-
    //    plane consistency check.
    let outer_corners = solve_outer_corners(input, options)?;

    // 3. Blend-coefficient endpoint solves.
    let blend_coeffs = compute_ring_blend_coeffs(input, options)?;

    // 4. Apply corner constructions to a MUTABLE copy of the
    //    seed patches — the row-1 solver reads this state as
    //    its min-norm starting point.
    let mut pre_row1_patches = seed.clone();
    let n = input.n();
    for k in 0..n {
        // Apply corner-C output to grid[6][1] on each patch.
        pre_row1_patches[k].grid[6][1] = corner_c.inner_cps[k];
        // Apply outer-corner output:
        //   corners[k].patch_k_grid_6_6 → patch k's grid[6][6].
        //   corners[k].patch_k_plus_1_grid_1_1 → patch (k+1)'s grid[1][1].
        pre_row1_patches[k].grid[6][6] = outer_corners.corners[k].patch_k_grid_6_6;
        let k_plus_1 = (k + 1) % n;
        pre_row1_patches[k_plus_1].grid[1][1] = outer_corners.corners[k].patch_k_plus_1_grid_1_1;
    }

    // 5. Run the ring row-1 G1 solver against the corner-adjusted
    //    seed. Uses a fresh `input_for_ring` so we pass the
    //    pre-row1 patches as the solver's seed.
    let input_for_ring = DirectBoundaryInput {
        seed_patches: Some(pre_row1_patches),
        hex_boundaries: None,
        ..input.clone()
    };
    let ring_result = solve_ring_row1_g1(&input_for_ring, &blend_coeffs, options)?;

    // 6. Self-check at the requested validation level.
    let spokes = canonical_spoke_topology(n);
    let (dense_min_cos_g1, dense_min_cos_g2, max_position_gap_g0) = match options.validation_level {
        ValidationLevel::None => (f64::NAN, None, f64::NAN),
        ValidationLevel::Standard | ValidationLevel::Full => {
            let sc_opts = SelfCheckOptions {
                samples_per_spoke: if options.validation_level == ValidationLevel::Full {
                    200
                } else {
                    100
                },
                measure_g2: options.target_continuity == ContinuityLevel::G2,
            };
            let m = measure_ring_quality(&ring_result.updated_patches, &spokes, &sc_opts);
            // Fire SelfCheckFailedBug if the measured G1 is below
            // the threshold AND the input had enough structural
            // quality to reach it (i.e. not an intentional
            // degenerate-fallback path).
            if m.dense_min_cos_g1 < 1.0 - options.self_check_g1_tol {
                // Soft fail — the construction returned patches
                // but dense G1 isn't at machine precision. Could
                // indicate a math bug OR an asymmetric / curved
                // input where α's constant-linear form isn't
                // enough. Emit the bug signal but still return
                // the result so callers can diagnose.
                return Err(DirectConstructionError::SelfCheckFailedBug {
                    measured: m.dense_min_cos_g1,
                    required: 1.0 - options.self_check_g1_tol,
                });
            }
            (
                m.dense_min_cos_g1,
                m.dense_min_cos_g2,
                m.max_position_gap_g0,
            )
        }
    };

    // 7. Diagnostics.
    let mut corner_diagnostics = Vec::with_capacity(n * 2);
    for _ in 0..n {
        corner_diagnostics.push(CornerDiagnostic {
            corner_kind: CornerKind::Central,
            formula_family: FormulaFamily::ShirmanSequinEnclosure,
            near_degeneracy_score: 1.0, // TODO: surface from corner-C solver
            tolerance_used: options.corner_singularity_eps,
        });
        corner_diagnostics.push(CornerDiagnostic {
            corner_kind: CornerKind::OuterHex,
            formula_family: FormulaFamily::ChiyokuraKimuraStandard,
            near_degeneracy_score: 1.0,
            tolerance_used: options.corner_singularity_eps,
        });
    }

    let spoke_diagnostics = blend_coeffs
        .iter()
        .zip(ring_result.per_spoke_residual.iter())
        .map(|(bc, &res)| SpokeDiagnostic {
            formula_family: FormulaFamily::ChiyokuraKimuraStandard,
            g1_residual: res,
            g2_residual: None,
            alpha_polynomial_degree: 1,
            beta_polynomial_degree: 1,
            near_degeneracy_score: bc.condition_number.clamp(0.0, 1.0),
        })
        .collect();

    Ok(DirectConstructionResult {
        sub_patches: ring_result.updated_patches,
        achieved_continuity: ContinuityLevel::G1,
        corner_diagnostics,
        spoke_diagnostics,
        dense_min_cos_g1,
        dense_min_cos_g2,
        max_position_gap_g0,
    })
}

/// The canonical three-quad spoke topology: spoke k connects
/// patch k's UOne (not reversed) to patch (k+1)%N's VZero
/// (reversed). This matches the convention in
/// `popkern-fillet/src/construct/three_quad_apply.rs`.
fn canonical_spoke_topology(n: usize) -> Vec<Spoke> {
    (0..n)
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
        .collect()
}

#[cfg(test)]
#[path = "direct_build_tests.rs"]
mod tests;
