//! Promotion ladder for the three-quad decomposition.
//!
//! Walks the rung hierarchy:
//! deg5 → biquintic-full → biseptic → two-stripe → KP13 bi-5
//! (env-gated) → KP20 bi-6 macro (env-gated) → bi-6 polar →
//! KP binonic.
//!
//! Extracted from `three_quad_apply.rs` to keep that file under
//! the 500-line structural limit.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::gregory::GregoryPatch;

use super::three_quad_kp_diagnostics;
use super::three_quad_promotions::{
    promote_decomp_to_biquintic_full_with_spoke_g1, promote_decomp_to_biseptic_with_spoke_g1,
    promote_decomp_to_two_stripe_with_spoke_g1,
};
use crate::blend::boundary::{BlendBoundary, BoundarySegment};
use crate::blend::gregory_three_quad::NQuadDecomposition;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ApexConstructionDiagnostic, CornerSubPatchFamily, CornerSubPatchProvenance};

/// Result of running the promotion ladder.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(super) enum PromotionResult {
    /// Caller should use the ordinary N-way poke.
    #[default]
    Standard,
    /// Caller should use the 4N macro-cap poke with this construction family.
    Macro4n(CornerSubPatchFamily),
    /// Caller should use the expanded KP20 shared-trim support topology.
    ExpandedKp20Support,
    /// Caller should use the KP20 support-cell topology.
    Kp20SupportCell,
}

#[derive(Debug, Clone, Default)]
pub(crate) struct KpBinonicPromotionReport {
    pub relaxed_coplanarity: Option<f64>,
    pub diagnostics: Vec<ApexConstructionDiagnostic>,
}

#[derive(Debug, Clone, Default)]
pub(super) struct PromotionLadderReport {
    pub result: PromotionResult,
    pub kp20_support_cell_provenance: Option<Vec<CornerSubPatchProvenance>>,
    pub kp20_support_cell_floor_patches: Option<Vec<GregoryPatch>>,
    pub diagnostics: Vec<ApexConstructionDiagnostic>,
}

/// Exact bi-3 (cubic) BB control points of a boundary segment curve:
/// endpoints interpolated, the two interior CPs least-squares fit over
/// 40 samples. The fillet's `FilletEnd`/`FaceEdge` segments are
/// cubic-exact (PK-087 stage-6 measurement: best-cubic residual
/// ~1e-16), so this recovers the curve's control polygon rather than
/// approximating it — the accurate replacement for the spoke extractor's
/// crude on-curve degree reduction.
fn segment_bi3_cps(curve: &popkern_core::geometry::NurbsCurve) -> [Point3; 4] {
    let (t0, t1) = curve.parameter_range();
    let m = 40usize;
    let pts: Vec<Point3> = (0..=m)
        .map(|i| curve.point_at(t0 + (i as f64 / m as f64) * (t1 - t0)))
        .collect();
    let p0 = pts[0];
    let p3 = pts[m];
    let basis = |s: f64| {
        let u = 1.0 - s;
        [u * u * u, 3.0 * u * u * s, 3.0 * u * s * s, s * s * s]
    };
    let mut ata = [[0.0f64; 2]; 2];
    let mut atb = [[0.0f64; 3]; 2];
    for (i, pt) in pts.iter().enumerate() {
        let s = i as f64 / m as f64;
        let b = basis(s);
        let bi = [b[1], b[2]];
        let resid = [
            pt.x - p0.x * b[0] - p3.x * b[3],
            pt.y - p0.y * b[0] - p3.y * b[3],
            pt.z - p0.z * b[0] - p3.z * b[3],
        ];
        for a in 0..2 {
            for c in 0..2 {
                ata[a][c] += bi[a] * bi[c];
            }
            for k in 0..3 {
                atb[a][k] += bi[a] * resid[k];
            }
        }
    }
    let det = ata[0][0] * ata[1][1] - ata[0][1] * ata[1][0];
    if det.abs() < 1e-18 {
        return [p0, p0, p3, p3];
    }
    let solve = |k: usize| {
        (
            (atb[0][k] * ata[1][1] - atb[1][k] * ata[0][1]) / det,
            (ata[0][0] * atb[1][k] - ata[1][0] * atb[0][k]) / det,
        )
    };
    let (c1x, c2x) = solve(0);
    let (c1y, c2y) = solve(1);
    let (c1z, c2z) = solve(2);
    [
        p0,
        Point3::new(c1x, c1y, c1z),
        Point3::new(c2x, c2y, c2z),
        p3,
    ]
}

/// Per-sector two-outer-edge fit targets from the real blend boundary
/// (PK-087 stage 6). For cap sector `k` the two outer edges are
/// `FilletEnd[k]` (`segments[2k]`) and `FaceEdge[k]` (`segments[2k+1]`),
/// which meet at the mid-sector junction `J_k`. Both edges are oriented
/// so CP0 is `J_k` (the shared corner where the cap's `tad`/`tak` edges
/// meet), matching the fit's `b[0][0]` corner. Returns `None` unless the
/// boundary has exactly `2n` segments.
fn sector_outer_edge_targets(boundary: &BlendBoundary, n: usize) -> Option<Vec<[[Point3; 4]; 2]>> {
    if boundary.segments.len() != 2 * n {
        return None;
    }
    let curve_of = |seg: &BoundarySegment| -> popkern_core::geometry::NurbsCurve {
        match seg {
            BoundarySegment::FilletEnd { curve, .. } | BoundarySegment::FaceEdge { curve, .. } => {
                curve.clone()
            }
        }
    };
    let rev = |cps: [Point3; 4]| [cps[3], cps[2], cps[1], cps[0]];
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let fe = segment_bi3_cps(&curve_of(&boundary.segments[2 * k])); // W_{k-1} -> J_k
        let face = segment_bi3_cps(&curve_of(&boundary.segments[2 * k + 1])); // J_k -> W_k
        // tad edge (script g[s][i][0]) ← FaceEdge (CP0 = J_k already);
        // tak edge (script g[s][0][i]) ← FilletEnd reversed (CP0 = J_k).
        out.push([face, rev(fe)]);
    }
    Some(out)
}

/// Run the full promotion ladder and return its selected topology
/// family plus structured diagnostics for declined promotion rungs.
///
/// Phase 6c-B: `promote_binonic` callback takes both `sub_patches`
/// and `fillet_surfaces` so the binonic optimiser can use outer-
/// fillet G1 residuals (cap-vs-fillet tangent-plane match along
/// each cap's outer boundary).
///
/// PK-090: `apex_boundary` is the apex blend boundary (plus the BRep
/// that owns its face-edge surfaces) for the honest 2n-segment
/// blend-boundary paper candidate in the KP20 polar dispatch; call
/// sites without one pass `None`.
#[allow(clippy::too_many_arguments)]
pub(super) fn run_promotion_ladder(
    decomp: &mut NQuadDecomposition,
    fillet_surfaces: &[FilletSurfaceResult],
    apex_boundary: Option<(&crate::blend::BlendBoundary, &popkern_core::topology::BRep)>,
    relaxed_coplanarity: &mut Option<f64>,
    use_deg5: bool,
    use_biquintic_full: bool,
    use_biseptic: bool,
    use_two_stripe: bool,
    use_kp_binonic: bool,
    use_kp_polar: bool,
    promote_binonic: impl FnOnce(
        &mut [GregoryPatch],
        &[FilletSurfaceResult],
    ) -> KpBinonicPromotionReport,
) -> PromotionLadderReport {
    let n_sides = decomp.sub_patches.len();
    if !(use_biquintic_full && use_deg5) {
        return standard_report_with_prerequisite_decline(
            use_kp_binonic,
            use_kp_polar,
            decomp.sub_patches.len(),
            "pre_biquintic_disabled",
            Some("deg5 or biquintic-full promotion disabled".to_string()),
        );
    }
    promote_decomp_to_biquintic_full_with_spoke_g1(&mut decomp.sub_patches);
    if !use_biseptic || !use_two_stripe {
        return standard_report_with_prerequisite_decline(
            use_kp_binonic,
            use_kp_polar,
            decomp.sub_patches.len(),
            "pre_two_stripe_disabled",
            Some("biseptic or two-stripe promotion disabled".to_string()),
        );
    }
    if !promote_decomp_to_biseptic_with_spoke_g1(&mut decomp.sub_patches) {
        return standard_report_with_prerequisite_decline(
            use_kp_binonic,
            use_kp_polar,
            decomp.sub_patches.len(),
            "biseptic_unavailable",
            None,
        );
    }
    if !promote_decomp_to_two_stripe_with_spoke_g1(&mut decomp.sub_patches) {
        return standard_report_with_prerequisite_decline(
            use_kp_binonic,
            use_kp_polar,
            decomp.sub_patches.len(),
            "two_stripe_unavailable",
            None,
        );
    }
    if use_kp_polar
        && decomp.sub_patches.len() >= 3
        && let Some(patches_4n) =
            super::kp13_bi5_promote::try_kp13_bi5_dispatch(&decomp.sub_patches, fillet_surfaces)
    {
        decomp.sub_patches = patches_4n;
        return PromotionLadderReport {
            result: PromotionResult::Macro4n(CornerSubPatchFamily::Kp13Bi5),
            kp20_support_cell_provenance: None,
            kp20_support_cell_floor_patches: None,
            diagnostics: Vec::new(),
        };
    }
    let mut diagnostics = Vec::new();
    if use_kp_polar && decomp.sub_patches.len() >= 3 {
        let mut report = super::kp20_polar_dispatch::try_kp20_polar_dispatch_with_diagnostics(
            &decomp.sub_patches,
            fillet_surfaces,
            apex_boundary,
        );
        diagnostics.append(&mut report.diagnostics);
        if let Some(patches_4n) = report.patches {
            // Commercial-tier paper accepts record their measured
            // relative g2 residual as a relaxed-continuity fact — the
            // corner ships smooth-but-not-paper-exact and the
            // CornerContinuityRelaxed warning makes that visible.
            if let Some(g2_relative) = report.commercial_g2_relative {
                *relaxed_coplanarity =
                    Some(relaxed_coplanarity.map_or(g2_relative, |worst| worst.max(g2_relative)));
            }
            decomp.sub_patches = patches_4n;
            let result = if report.expanded_support_topology_accepted {
                PromotionResult::ExpandedKp20Support
            } else if report.support_cell_topology_accepted {
                PromotionResult::Kp20SupportCell
            } else {
                PromotionResult::Macro4n(CornerSubPatchFamily::Kp20Bi6)
            };
            return PromotionLadderReport {
                result,
                kp20_support_cell_provenance: report.support_cell_provenance,
                kp20_support_cell_floor_patches: report.support_cell_floor_patches,
                diagnostics,
            };
        }
    }
    let (polar_committed, mut polar_diagnostics) = if use_kp_polar && decomp.sub_patches.len() >= 5
    {
        let boundary_target =
            apex_boundary.and_then(|(b, _)| sector_outer_edge_targets(b, decomp.sub_patches.len()));
        super::kp_polar_bi6_promote::promote_and_update_relaxed_with_diagnostics(
            &mut decomp.sub_patches,
            fillet_surfaces,
            relaxed_coplanarity,
            boundary_target.as_deref(),
        )
    } else {
        (false, Vec::new())
    };
    diagnostics.append(&mut polar_diagnostics);
    if use_kp_binonic && !polar_committed {
        let report = promote_binonic(&mut decomp.sub_patches, fillet_surfaces);
        *relaxed_coplanarity = report.relaxed_coplanarity;
        diagnostics.extend(report.diagnostics);
    }
    if use_kp_polar && n_sides >= 5 {
        diagnostics.push(
            three_quad_kp_diagnostics::support_cell_candidate_nonselection(
                "promotion_ladder_standard",
                n_sides,
                decomp.sub_patches.len(),
            ),
        );
    }
    PromotionLadderReport {
        result: PromotionResult::Standard,
        kp20_support_cell_provenance: None,
        kp20_support_cell_floor_patches: None,
        diagnostics,
    }
}

fn standard_report_with_prerequisite_decline(
    use_kp_binonic: bool,
    use_kp_polar: bool,
    n_sides: usize,
    reason_code: &'static str,
    detail: Option<String>,
) -> PromotionLadderReport {
    let mut diagnostics = Vec::new();
    if use_kp_binonic {
        diagnostics.push(three_quad_kp_diagnostics::promotion_declined(
            reason_code,
            n_sides,
            detail,
            None,
            None,
        ));
    }
    if use_kp_polar && n_sides >= 5 {
        diagnostics.push(
            three_quad_kp_diagnostics::support_cell_candidate_nonselection(
                reason_code,
                n_sides,
                n_sides,
            ),
        );
    }
    PromotionLadderReport {
        result: PromotionResult::Standard,
        kp20_support_cell_provenance: None,
        kp20_support_cell_floor_patches: None,
        diagnostics,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;

    #[test]
    fn binonic_prerequisite_decline_is_structured_diagnostic() {
        let mut decomp = NQuadDecomposition {
            central_point: Point3::new(0.0, 0.0, 0.0),
            sub_patches: Vec::new(),
        };
        let mut relaxed_coplanarity = None;
        let report = run_promotion_ladder(
            &mut decomp,
            &[],
            None,
            &mut relaxed_coplanarity,
            false,
            false,
            true,
            true,
            true,
            false,
            |_, _| KpBinonicPromotionReport::default(),
        );
        assert_eq!(report.result, PromotionResult::Standard);
        assert_eq!(report.diagnostics.len(), 1);
        let crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(cap) =
            &report.diagnostics[0]
        else {
            panic!("expected exact-boundary cap diagnostic");
        };
        assert_eq!(cap.code, "kp_binonic.promotion.declined");
        assert_eq!(cap.stage, "kp_binonic.promotion.pre_biquintic_disabled");
        assert_eq!(cap.sector_name, "kp_binonic_ring_n0");
    }
}
