//! Phase I.16 production dispatch — `popkern-fillet` entry point
//! for the KP13 bi-5 macro-patch construction.
//!
//! Mirrors the KP15 bi-6 polar wrapper [`try_kp_polar_bi6_cap`]
//! (`kp_polar_bi6_apply.rs`), but targets the KP13 bi-5 pipeline
//! (`g△ → g → a → p → bi-9`, paper-faithful through I.14.e) from
//! the same [`PolarBoundaryInput`] input type.
//!
//! # Representation adapter
//!
//! The KP15 bi-6 path consumes `kp_polar_boundary_input::TensorBorder`
//! (**bi-2 in v**: 3 rows × 4 CPs per sector = position + 1st
//! cross + 2nd cross derivatives as BB coefs at `v_polar ∈ {0,
//! 1/2, 1}`). The KP13 bi-5 path consumes
//! `kp13_bi5_boundary::TensorBorder` (4 × 4 BB CPs per sector —
//! a **full bi-3 tensor-product grid**).
//!
//! Coordinate reconciliation: polar's `v_polar = 0` is the outer
//! fillet mouth (furthest from the star); KP13's `v = 0` is the
//! star edge (corner `(0, 0)`). The two v-axes run in **opposite
//! directions**: `v_kp13 = 1 − v_polar`.
//!
//! The adapter is therefore a (v-reversal + bi-2 → bi-3 degree
//! elevation) acting on the CPs per sector:
//!
//! ```text
//!   D_0 (KP13 v = 0, star edge)   = C_2 = second_cross_row
//!   D_1 (KP13 v = 1/3)            = (1/3) · C_2 + (2/3) · C_1
//!   D_2 (KP13 v = 2/3)            = (2/3) · C_1 + (1/3) · C_0
//!   D_3 (KP13 v = 1, outer edge)  = C_0 = position_row
//! ```
//!
//! where `C_i = polar row i` (C_0 = position, C_1 = first cross,
//! C_2 = second cross). The middle two CPs are the standard
//! bi-2 → bi-3 degree-elevation formulas (after accounting for
//! the v-reversal).
//!
//! This conversion is **exact** up to the bi-2 expressiveness of
//! the polar input — the returned bi-3 grid represents the same
//! polynomial in v, now in 4-CP form. No extrapolation or
//! information loss.
//!
//! # Usage
//!
//! Call sites that already have a [`PolarBoundaryInput`] (the
//! natural input type for the polar dispatch chain in
//! [`super::kp_polar_bi6_apply`]) can invoke
//! [`try_kp_bi5_macro_cap`] to get **`4n` bi-9 CP grids** — one
//! per `(sector, piece)` in the KP13 2 × 2 sub-patch layout.
//!
//! # Note on output shape
//!
//! Unlike the bi-6 polar path (`N` patches per cap), KP13 bi-5
//! produces `4N` patches. Callers that BRep-wrap the output must
//! handle this 4× multiplicity — the paper's 2 × 2 piece layout
//! is a structural property of the bi-5 macro construction, not
//! a representation detail.

#![allow(dead_code)]

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Bi5CapBi3Patch, Bi5CapInput, Bi5CapTensorBorder, PolarBoundaryInput, RectangularGregoryBiNonic,
    TensorBorder, assemble_bi5_macro_cap, catmull_clark_star_point,
    extract_bi3_tensor_border_from_polar_input,
};

/// Degree-7 Bezier evaluation on a spoke's 8-CP boundary. Used to
/// sample the spokes at `t = 1/3` and `t = 2/3` for the
/// Catmull-Clark cap-net inputs `c5`, `c6`. Mirrors the same
/// helper in the bi-6 polar bridge.
fn eval_degree7_bezier(cps: &[Point3; 8], t: f64) -> Point3 {
    let s = 1.0 - t;
    // Bi-7 Bernstein: B^7_i(t) = C(7, i) t^i (1−t)^(7−i).
    let binom: [f64; 8] = [1.0, 7.0, 21.0, 35.0, 35.0, 21.0, 7.0, 1.0];
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for i in 0..8 {
        let w = binom[i] * t.powi(i as i32) * s.powi(7 - i as i32);
        x += w * cps[i].x;
        y += w * cps[i].y;
        z += w * cps[i].z;
    }
    Point3::new(x, y, z)
}

/// Convert a `kp_polar_boundary_input::TensorBorder` (bi-3,
/// 3 rows × 4 CPs per sector) to a `kp13_bi5_boundary::TensorBorder`
/// (4 × 4 full bi-3 grid per sector).
///
/// See module header for the row-mapping conventions and the
/// "row 3 via linear extrapolation" simplification.
fn tensor_border_kp_polar_to_kp13(polar_tb: &TensorBorder) -> Option<Bi5CapTensorBorder> {
    // Require the polar input to be bi-3 (not bi-5 — that's for
    // the ternary-split variant). Return None if bi-5 is passed.
    use popkern_core::geometry::gregory::TensorBorderKind;
    if polar_tb.kind() != TensorBorderKind::Bi3 {
        return None;
    }

    let n = polar_tb.n_sectors();
    let mut patches = Vec::with_capacity(n);
    for k in 0..n {
        let r0 = polar_tb.position_row(k); // C_0 — polar v = 0 (outer mouth).
        let r1 = polar_tb.first_cross_row(k); // C_1 — polar v = 1/2.
        let r2 = polar_tb.second_cross_row(k); // C_2 — polar v = 1 (deepest inward).
        // Convert bi-2 (polar) → bi-3 (KP13) in v with v_kp13 =
        // 1 − v_polar, per module-header derivation.
        //   D_0 = C_2 (star side)
        //   D_1 = (1/3) C_2 + (2/3) C_1
        //   D_2 = (2/3) C_1 + (1/3) C_0
        //   D_3 = C_0 (outer boundary)
        let mut cps = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
        for i in 0..4 {
            let c0 = r0[i];
            let c1 = r1[i];
            let c2 = r2[i];
            cps[i][0] = c2;
            cps[i][1] = Point3::new(
                c2.x / 3.0 + 2.0 * c1.x / 3.0,
                c2.y / 3.0 + 2.0 * c1.y / 3.0,
                c2.z / 3.0 + 2.0 * c1.z / 3.0,
            );
            cps[i][2] = Point3::new(
                2.0 * c1.x / 3.0 + c0.x / 3.0,
                2.0 * c1.y / 3.0 + c0.y / 3.0,
                2.0 * c1.z / 3.0 + c0.z / 3.0,
            );
            cps[i][3] = c0;
        }
        patches.push(Bi5CapBi3Patch { cps });
    }
    Some(Bi5CapTensorBorder { patches })
}

/// End-to-end KP13 bi-5 macro-patch construction from a
/// [`PolarBoundaryInput`].
///
/// Runs `g△ → g → a → p → bi-9 elevation` and returns the `4N`
/// bi-9 CP grids (paper's 2 × 2 sub-patch layout per sector). See
/// module header for representation-adapter details and the
/// output-shape caveat.
///
/// Errors on: invalid polar input, `n < 3`, tensor-border kind
/// mismatch (bi-5 input rejected), or any singular linear system
/// in the assembly (LSQ or F_1 min).
pub fn try_kp_bi5_macro_cap(
    polar_input: &PolarBoundaryInput,
) -> Result<Vec<[[Point3; 10]; 10]>, String> {
    polar_input.validate()?;
    let n = polar_input.n();
    if n < 3 {
        return Err(format!("KP13 bi-5 cap needs n ≥ 3, got {n}"));
    }

    // CC limit point for the star: standard Catmull-Clark formula
    // with c7 = central vertex, c5s / c6s sampled at t = 1/3, 2/3
    // on each spoke.
    let c7 = polar_input.central_vertex;
    let mut c5s = Vec::with_capacity(n);
    let mut c6s = Vec::with_capacity(n);
    for spoke in &polar_input.spokes {
        c5s.push(eval_degree7_bezier(&spoke.boundary_cps, 1.0 / 3.0));
        c6s.push(eval_degree7_bezier(&spoke.boundary_cps, 2.0 / 3.0));
    }
    let cc_limit = catmull_clark_star_point(n, c7, &c5s, &c6s);

    // Extract bi-3 tensor border from polar input (reuses the
    // bi-6 path's extraction), then adapt to KP13's full-bi-3-grid
    // representation.
    let polar_tb = extract_bi3_tensor_border_from_polar_input(polar_input)?;
    let target = tensor_border_kp_polar_to_kp13(&polar_tb).ok_or_else(|| {
        "kp_bi5_macro_apply: adapter expects bi-3 tensor-border (bi-5 not supported)".to_string()
    })?;

    // Assemble the bi-5 macro cap end-to-end.
    let input = Bi5CapInput { target, cc_limit };
    let out = assemble_bi5_macro_cap(&input).map_err(|e| format!("KP13 bi-5 assembly: {e}"))?;
    Ok(out.patches_bi9)
}

/// Variant of [`try_kp_bi5_macro_cap`] that wraps each of the
/// `4n` bi-9 CP grids as a [`RectangularGregoryBiNonic`] — the
/// same BRep-insertable wrapper type used by the Phase H bi-9
/// direct and KP15 bi-6 polar paths.
///
/// # Output shape
///
/// Returns a flat `Vec<RectangularGregoryBiNonic>` of length
/// `4n`, in canonical `(sector, piece)` row-major order — i.e.,
/// `output[4k + piece_index]` for `piece_index ∈ 0..=3`
/// (`O, L, R, M`).
///
/// # Design choice
///
/// We treat each of the 4 sub-patches per sector as an
/// independent bi-9 Gregory patch rather than introducing a new
/// `RectangularGregoryBiNonicQuadPatchwork` wrapper type. The
/// sub-patches live on their own `[0, 1]²` piece-local domains
/// and are `C²`-connected internally (via I.14.d) and `G²`-
/// connected across sectors (via I.14.e eq (18)); the BRep
/// doesn't need to know about the 2 × 2 grouping.
pub fn try_kp_bi5_macro_cap_as_gregory(
    polar_input: &PolarBoundaryInput,
) -> Result<Vec<RectangularGregoryBiNonic>, String> {
    let grids = try_kp_bi5_macro_cap(polar_input)?;
    Ok(grids
        .into_iter()
        .map(|grid| RectangularGregoryBiNonic { grid })
        .collect())
}

#[cfg(test)]
#[path = "kp_bi5_macro_apply_tests.rs"]
mod tests;
