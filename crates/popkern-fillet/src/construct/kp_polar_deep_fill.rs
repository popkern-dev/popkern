//! Boundary-owned deep-interior candidate for polar bi-6 caps.
//!
//! The KP15 scaffold leaves most `j >= 3` control points at the
//! cap star unless the caller supplies `interior_filler`. That is
//! structurally too little ownership for high-valence caps: the
//! open interface edge `U=1` must be the next sector's row-0
//! curve, not a collapsed star row. This module builds a bounded
//! candidate that populates those cells before the nonlinear L-M
//! fallback runs.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Bi6CapFreeParams, Bi6CapInput, GregoryPatch, RectangularGregoryBiNonic, assemble_bi6_polar_cap,
    rule_s1_star_propagation,
};

use super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::three_quad_apply::worst_ring_coplanarity;

pub(super) enum DeepFillOutcome {
    Accepted(Vec<GregoryPatch>),
    Rejected(&'static str),
}

pub(super) fn try_boundary_owned_deep_fill(
    base: &Bi6CapInput,
    free: &Bi6CapFreeParams,
    commercial_gate: f64,
) -> DeepFillOutcome {
    let mut last_reason = "deep_fill_no_candidate";
    for scale in [1.0, 0.5, 0.25, 0.05, 0.0] {
        let candidate = contract_free_params(free, base.c7, scale);
        for boundary_weight in [0.0, 0.1, 0.25, 0.5, 0.75, 1.0] {
            match try_single_boundary_owned_deep_fill(
                base,
                &candidate,
                commercial_gate,
                boundary_weight,
            ) {
                DeepFillOutcome::Accepted(patches) => return DeepFillOutcome::Accepted(patches),
                DeepFillOutcome::Rejected(reason) => last_reason = reason,
            }
        }
    }
    DeepFillOutcome::Rejected(last_reason)
}

fn try_single_boundary_owned_deep_fill(
    base: &Bi6CapInput,
    free: &Bi6CapFreeParams,
    commercial_gate: f64,
    boundary_weight: f64,
) -> DeepFillOutcome {
    let Some(input) = boundary_owned_deep_fill_input(base, free, boundary_weight) else {
        return DeepFillOutcome::Rejected("deep_fill_input_shape");
    };
    let Ok(out) = assemble_bi6_polar_cap(&input) else {
        return DeepFillOutcome::Rejected("deep_fill_final_assembly");
    };
    let patches: Vec<GregoryPatch> = out
        .patches_bi9
        .into_iter()
        .map(|grid| GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid })))
        .collect();
    if patches
        .iter()
        .any(|patch| patch.self_intersects_with_tolerance(SELF_INT_TOL))
    {
        return DeepFillOutcome::Rejected("deep_fill_self_intersection");
    }
    if worst_ring_coplanarity(&patches) > commercial_gate {
        return DeepFillOutcome::Rejected("deep_fill_g1_exceeded");
    }
    DeepFillOutcome::Accepted(patches)
}

pub(super) fn boundary_owned_deep_fill_input(
    base: &Bi6CapInput,
    free: &Bi6CapFreeParams,
    boundary_weight: f64,
) -> Option<Bi6CapInput> {
    let mut input = apply_near_star_and_free_cells(base, free)?;
    let Ok(seed) = assemble_bi6_polar_cap(&input) else {
        return None;
    };

    let mut filler = input
        .interior_filler
        .take()
        .unwrap_or_else(|| vec![[[None; 7]; 7]; input.n]);
    for (k, filler_sector) in filler.iter_mut().enumerate().take(input.n) {
        let next = &seed.patches_bi6[(k + 1) % input.n];
        let sector = &seed.patches_bi6[k];
        for j in 3..=6 {
            let v = boundary_weight * ((j - 2) as f64) / 4.0;
            let right = next[6 - j][0];
            let left = sector[0][j];
            for i in 0..=6 {
                if i == 6 {
                    filler_sector[i][j] = Some(right);
                    continue;
                }
                let u = (i as f64) / 6.0;
                let row2 = sector[i][2];
                let owned_edge = lerp(left, right, u);
                filler_sector[i][j] = Some(lerp(row2, owned_edge, v));
            }
        }
        filler_sector[2][3] = Some(free.p_23[k]);
        filler_sector[3][3] = Some(free.p_33[k]);
    }
    input.interior_filler = Some(filler);
    Some(input)
}

fn apply_near_star_and_free_cells(
    base: &Bi6CapInput,
    free: &Bi6CapFreeParams,
) -> Option<Bi6CapInput> {
    if free.p_23.len() != base.n || free.p_33.len() != base.n {
        return None;
    }
    let mut input = base.clone();
    input.p0_22 = free.p0_22;
    if !input.sectors.is_empty() {
        input.sectors[0].p_01 = free.p0_01;
        input.sectors[0].p_02 = free.p0_02;
        input.sectors[0].p_11 = free.p0_11;
    }
    for k in 1..input.sectors.len() {
        let prev = input.sectors[k - 1];
        let s1 = rule_s1_star_propagation(
            &input.params,
            prev.p_00,
            prev.p_01,
            prev.p_02,
            prev.p_10,
            prev.p_11,
            prev.p_20,
        );
        input.sectors[k].p_01 = s1.p01;
        input.sectors[k].p_11 = s1.p11;
        input.sectors[k].p_02 = s1.p02;
    }
    let mut filler = vec![[[None; 7]; 7]; input.n];
    for (k, filler_sector) in filler.iter_mut().enumerate().take(input.n) {
        filler_sector[2][3] = Some(free.p_23[k]);
        filler_sector[3][3] = Some(free.p_33[k]);
    }
    input.interior_filler = Some(filler);
    Some(input)
}

fn lerp(a: Point3, b: Point3, t: f64) -> Point3 {
    a + (b - a) * t
}

fn contract_free_params(free: &Bi6CapFreeParams, anchor: Point3, scale: f64) -> Bi6CapFreeParams {
    let contract = |p: Point3| anchor + (p - anchor) * scale;
    Bi6CapFreeParams {
        p0_01: contract(free.p0_01),
        p0_02: contract(free.p0_02),
        p0_11: contract(free.p0_11),
        p0_22: contract(free.p0_22),
        p_23: free.p_23.iter().map(|p| contract(*p)).collect(),
        p_33: free.p_33.iter().map(|p| contract(*p)).collect(),
    }
}
