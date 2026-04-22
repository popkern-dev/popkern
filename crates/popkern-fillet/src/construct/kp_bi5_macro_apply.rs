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
//! (3 rows × 4 CPs per sector: position + 1st cross + 2nd cross
//! derivatives). The KP13 bi-5 path consumes
//! `kp13_bi5_boundary::TensorBorder` (4 × 4 BB CPs per sector —
//! a **full bi-3 tensor-product grid**).
//!
//! The adapter [`tensor_border_kp_polar_to_kp13`] maps the 3-row
//! depth-2 strip to the 4-row full bi-3 grid via:
//!
//! ```text
//!   row 0 (toward star)     = second_cross_row
//!   row 1                   = first_cross_row
//!   row 2                   = position_row
//!   row 3 (outer boundary)  = 2·position_row − first_cross_row
//!                             (linear extrapolation outward)
//! ```
//!
//! This reading treats the polar input's "position" row as the
//! outermost visible boundary (matching the paper's Fig. 2 which
//! shows the tensor-border as a depth-2 ring abutting the cap).
//! Row 3 is a linear continuation of the depth gradient — a
//! deliberate MVP simplification.
//!
//! **Caveat**: paper-faithful tensor-border extraction for the
//! KP13 pipeline would re-derive all 4 rows from the underlying
//! bi-3 boundary curve + face tangent planes without the
//! extrapolation. This is additive paper-fidelity work.
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
    Bi5CapBi3Patch, Bi5CapInput, Bi5CapTensorBorder, PolarBoundaryInput, TensorBorder,
    assemble_bi5_macro_cap, catmull_clark_star_point, extract_bi3_tensor_border_from_polar_input,
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
        let r0 = polar_tb.position_row(k);
        let r1 = polar_tb.first_cross_row(k);
        let r2 = polar_tb.second_cross_row(k);
        // Row layout: row 0 (star side) = r2, row 1 = r1, row 2 =
        // r0 (outermost position), row 3 = 2·r0 − r1 (linear
        // outward extrapolation).
        let mut cps = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
        for i in 0..4 {
            cps[i][0] = r2[i];
            cps[i][1] = r1[i];
            cps[i][2] = r0[i];
            cps[i][3] = Point3::new(
                2.0 * r0[i].x - r1[i].x,
                2.0 * r0[i].y - r1[i].y,
                2.0 * r0[i].z - r1[i].z,
            );
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

#[cfg(test)]
#[path = "kp_bi5_macro_apply_tests.rs"]
mod tests;
