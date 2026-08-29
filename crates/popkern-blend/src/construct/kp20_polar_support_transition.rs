//! Bounded support-transition sheet primitive for KP20 support-owned folds.

use popkern_core::geometry::gregory::{
    Kp20Bi6ValidationReport, Kp20MacroPatch, Kp20PieceId, TensorBorderBi3,
};
use popkern_core::geometry::{Point3, Vec3};

use super::super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_cap, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, is_better_summary, validation_is_clean,
};
use super::super::kp20_polar_support_reparameterization::Kp20SupportTransitionReport;

const TARGET_PIECE: Kp20PieceId = Kp20PieceId::NearMidpointLeft;

pub(super) fn refine_support_transition_sheet(
    tensor_border: &TensorBorderBi3,
    macro_patches: Vec<Kp20MacroPatch>,
) -> Result<(Vec<Kp20MacroPatch>, Kp20SupportTransitionReport), String> {
    let baseline = candidate_summary_from_cap(tensor_border, macro_patches.clone())?;
    let mut report = Kp20SupportTransitionReport {
        evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds: baseline.0.fold_count,
        end_folds: baseline.0.fold_count,
        target_piece: TARGET_PIECE,
        alpha: 0.0,
        normal_scale: 0.0,
        skew_scale: 0.0,
        decline_reason: "support_transition_no_candidate",
    };
    let mut best: Option<(
        Vec<Kp20MacroPatch>,
        (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
        Kp20SupportTransitionParams,
    )> = None;
    for params in Kp20SupportTransitionParams::family() {
        report.evaluated += 1;
        let mut candidate = macro_patches.clone();
        let changed = apply_support_transition(&mut candidate, &baseline.0, params);
        if !changed.iter().any(|pieces| pieces[TARGET_PIECE.index()]) {
            report.rejected += 1;
            report.decline_reason = "support_transition_no_target_piece_folds";
            continue;
        }
        let fast = fast_candidate_summary_from_changed_patches(
            tensor_border,
            candidate.clone(),
            &baseline.0,
            &changed,
            baseline.0.fold_count,
        )?;
        if !validation_is_clean(fast.1) {
            report.rejected += 1;
            report.decline_reason = "support_transition_fast_validation_not_clean";
            continue;
        }
        if !is_better_summary(&fast, Some(&baseline)) {
            report.rejected += 1;
            report.decline_reason = "support_transition_no_fast_fold_reduction";
            continue;
        }
        let exact = candidate_summary_from_changed_patches(
            tensor_border,
            candidate.clone(),
            &baseline.0,
            &changed,
            baseline.0.fold_count,
        )?;
        if !validation_is_clean(exact.1) {
            report.rejected += 1;
            report.decline_reason = "support_transition_exact_validation_not_clean";
            continue;
        }
        if !is_better_summary(&exact, Some(&baseline)) {
            report.rejected += 1;
            report.decline_reason = "support_transition_no_exact_fold_reduction";
            continue;
        }
        if best
            .as_ref()
            .is_none_or(|(_, summary, _)| is_better_summary(&exact, Some(summary)))
        {
            best = Some((candidate, exact, params));
        } else {
            report.rejected += 1;
        }
    }
    let Some((candidate, summary, params)) = best else {
        report.rejected = report.evaluated;
        return Ok((macro_patches, report));
    };
    report.accepted = 1;
    report.rejected = report.evaluated - report.accepted;
    report.end_folds = summary.0.fold_count;
    report.alpha = params.alpha;
    report.normal_scale = params.normal_scale;
    report.skew_scale = params.skew_scale;
    report.decline_reason = "none";
    Ok((candidate, report))
}

#[derive(Debug, Clone, Copy)]
struct Kp20SupportTransitionParams {
    alpha: f64,
    normal_scale: f64,
    skew_scale: f64,
}

impl Kp20SupportTransitionParams {
    const fn family() -> [Self; 8] {
        [
            Self::new(0.25, 0.00, -0.08),
            Self::new(0.25, 0.00, 0.08),
            Self::new(0.40, 0.00, -0.08),
            Self::new(0.40, 0.00, 0.08),
            Self::new(0.25, -0.03, 0.00),
            Self::new(0.25, 0.03, 0.00),
            Self::new(0.40, -0.03, 0.00),
            Self::new(0.40, 0.03, 0.00),
        ]
    }

    const fn new(alpha: f64, normal_scale: f64, skew_scale: f64) -> Self {
        Self {
            alpha,
            normal_scale,
            skew_scale,
        }
    }
}

fn apply_support_transition(
    macro_patches: &mut [Kp20MacroPatch],
    baseline: &Kp20SamplerFoldSummary,
    params: Kp20SupportTransitionParams,
) -> Vec<[bool; 4]> {
    let mut changed = vec![[false; 4]; macro_patches.len()];
    for site in baseline
        .fold_sites
        .iter()
        .filter(|site| site.piece == TARGET_PIECE)
    {
        let patch = macro_patches[site.sector].piece_mut(TARGET_PIECE);
        transition_piece(patch, params);
        changed[site.sector][TARGET_PIECE.index()] = true;
    }
    changed
}

fn transition_piece(
    patch: &mut popkern_core::geometry::gregory::Kp20Bi6Patch,
    params: Kp20SupportTransitionParams,
) {
    let original = patch.cps;
    let frame = support_transition_frame(&original);
    let span = frame
        .boundary_tangent
        .length()
        .max(frame.support_inward.length());
    for i in 2..=4 {
        for j in 2..=4 {
            let u = i as f64 / 6.0;
            let v = j as f64 / 6.0;
            let coons = coons_point(&original, u, v, i, j);
            let weight = transition_weight(u, v);
            let displacement = frame.normal * (params.normal_scale * span * weight)
                + frame.cross_sector * (params.skew_scale * span * weight);
            patch.cps[i][j] =
                original[i][j] + (coons - original[i][j]) * params.alpha + displacement;
        }
    }
}

struct Kp20SupportTransitionFrame {
    boundary_tangent: Vec3,
    support_inward: Vec3,
    cross_sector: Vec3,
    normal: Vec3,
}

fn support_transition_frame(cps: &[[Point3; 7]; 7]) -> Kp20SupportTransitionFrame {
    let boundary_tangent = cps[0][6] - cps[0][0];
    let support_inward = cps[6][3] - cps[0][3];
    let cross_sector = cps[3][6] - cps[3][0];
    let normal = boundary_tangent
        .cross(&support_inward)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    Kp20SupportTransitionFrame {
        boundary_tangent,
        support_inward,
        cross_sector,
        normal,
    }
}

fn coons_point(cps: &[[Point3; 7]; 7], u: f64, v: f64, i: usize, j: usize) -> Point3 {
    let uv_edges = lerp_point(cps[0][j], cps[6][j], u);
    let vu_edges = lerp_point(cps[i][0], cps[i][6], v);
    let bottom = lerp_point(cps[0][0], cps[6][0], u);
    let top = lerp_point(cps[0][6], cps[6][6], u);
    let corners = lerp_point(bottom, top, v);
    uv_edges + (vu_edges - corners)
}

fn lerp_point(a: Point3, b: Point3, t: f64) -> Point3 {
    a + (b - a) * t
}

fn transition_weight(u: f64, v: f64) -> f64 {
    let du = (u - 0.5).abs() * 2.0;
    let dv = (v - 0.5).abs() * 2.0;
    (1.0 - du.max(dv)).max(0.0).powi(2)
}
