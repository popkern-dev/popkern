//! K-P bi-nonic promotion helpers for three-quad decompositions.

use super::*;
use crate::construct::kp_binonic_optimizer;
use crate::construct::kp_biquintic_optimizer::optimize_biquintic_and_build_binonic;
use crate::construct::kp_twist_optimizer::optimize_c_twists_and_build_binonic;
use crate::construct::promotion_ladder::KpBinonicPromotionReport;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;

#[path = "three_quad_kp_candidate_selection.rs"]
mod three_quad_kp_candidate_selection;
#[path = "three_quad_kp_ring_metrics.rs"]
mod three_quad_kp_ring_metrics;
pub(crate) use three_quad_kp_ring_metrics::{worst_ring_coplanarity, worst_ring_g0_gap};

/// Promote sub-patches from `RectangularTwoStripe` to
/// `RectangularBiNonic` via the Karciauskas-Peters guided-surface
/// construction. Returns a report containing structured decline
/// diagnostics when the promotion is rejected; accepted commercial
/// grade promotion mutates `sub_patches` in place.
///
/// Phase 6c-B: `fillet_surfaces` (one entry per sub-patch, ordered
/// by ring index) are forwarded to `optimize_binonic_with_outer_fillets`
/// so the bi-9 L-M solver includes per-sample cap-vs-fillet G1
/// residuals along each cap's outer boundary, in addition to the
/// inner-spoke G1 residuals. Closes the cap-fillet seam tangent-
/// plane mismatch documented in Task #144.
///
/// `POPKERN_OUTER_FILLET_G1=1` env var gates the new path during
/// the spike. Without it, behaviour matches the pre-Phase-6c-B
/// optimizer (inner-spoke G1 only).
pub(super) fn promote_decomp_to_kp_binonic(
    sub_patches: &mut [GregoryPatch],
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    hex_boundaries: Option<&[([Point3; 10], [Point3; 10])]>,
) -> KpBinonicPromotionReport {
    let n = sub_patches.len();
    if n == 0 {
        return KpBinonicPromotionReport::default();
    }
    // Every sub-patch must already be two-stripe — we read the
    // grid as the source of the K-P DirectBoundaryInput.
    let mut two_stripe: Vec<RectangularGregoryTwoStripe> = Vec::with_capacity(n);
    for patch in sub_patches.iter() {
        let GregoryPatch::RectangularTwoStripe(ts) = patch else {
            return declined("non_two_stripe_input", n, None, None, None);
        };
        two_stripe.push((**ts).clone());
    }

    // Ring-topology spokes: consistent with the other promotion
    // rungs — spoke `k` shares sub-patch `k`'s UOne with sub-patch
    // `(k+1) % n`'s VZero (reversed).
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

    let Some(input) = extract_direct_boundary_input(&two_stripe, &spokes, hex_boundaries) else {
        return declined("boundary_input_unavailable", n, None, None, None);
    };
    let char_len = estimate_characteristic_length(&two_stripe);
    let options = DirectConstructionOptions::with_characteristic_length(char_len);
    // Phase H.7 gap #2 closure — K-P 2004 §4 iterative twist
    // optimization refines the N central-vertex twists against the
    // ring-wide G1 coplanarity residual (the one-shot Shirman-Séquin
    // twist is under-conditioned at extreme radius asymmetry).
    // Option C step 3: three-phase optimizer (bicubic → bi-quintic
    // → bi-nonic); phase 3 engages only when phases 1+2 fail to
    // converge, and the ladder falls through to the bi-quintic-only
    // then bicubic C-twist optimizers if a phase is unavailable.
    // Phase 6c-B: env-gate on POPKERN_OUTER_FILLET_G1=1. When set,
    // pass the fillet surfaces to the binonic optimiser so its
    // L-M residual includes per-sample cap-vs-fillet G1 along
    // each cap's outer boundary.
    let outer_fillets: Option<Vec<popkern_core::geometry::surface::SurfaceGeometry>> =
        if std::env::var("POPKERN_OUTER_FILLET_G1").is_ok_and(|v| v == "1") {
            Some(fillet_surfaces.iter().map(|f| f.surface.clone()).collect())
        } else {
            None
        };

    let pre_worst = worst_ring_coplanarity(sub_patches);
    let mut optimization_input = input.clone();
    optimization_input.hex_boundaries = None;
    let (guide_binonic, guide_rung_outcome) = match optimize_kp_binonic_candidate(
        &optimization_input,
        &options,
        outer_fillets.as_deref(),
    ) {
        Some(candidate) => candidate,
        None => {
            return declined("candidate_unavailable", n, None, None, None);
        }
    };
    // PK-275 rung observability: record the armed outer-fillet
    // rung's outcome as a structured fact, whatever it was — the
    // rung previously declined at six optimizer sites with no
    // record, which made its (zero) coverage unverifiable.
    let mut rung_facts: Vec<crate::types::ApexConstructionDiagnostic> = Vec::new();
    if let Some(outcome) = &guide_rung_outcome {
        rung_facts.push(
            crate::construct::three_quad_kp_diagnostics::outer_fillet_rung_outcome(
                "guide", outcome, n,
            ),
        );
    }
    let gate_fillets: Vec<popkern_core::geometry::surface::SurfaceGeometry> =
        fillet_surfaces.iter().map(|f| f.surface.clone()).collect();
    let selected = three_quad_kp_candidate_selection::select_boundary_candidate(
        guide_binonic,
        &input,
        &options,
        outer_fillets.as_deref(),
        &gate_fillets,
        pre_worst,
        n,
    );
    let mut binonic = selected.binonic;
    let candidate_provenance = selected.provenance;
    let selection_fact = selected.selection_fact;
    if let Some(outcome) = &selected.exact_rung_outcome {
        rung_facts.push(
            crate::construct::three_quad_kp_diagnostics::outer_fillet_rung_outcome(
                "exact", outcome, n,
            ),
        );
    }

    // Self-intersection gate — never commit a fold.
    let binonic_gp: Vec<GregoryPatch> = binonic
        .iter()
        .cloned()
        .map(|p| GregoryPatch::RectangularBiNonic(Box::new(p)))
        .collect();
    for (i, gp) in binonic_gp.iter().enumerate() {
        if gp.self_intersects_with_tolerance(1e-6) {
            // PK-240: the candidate provenance discriminates the fold
            // classes — "guide" means every optimizer arm's product
            // folds (the extreme-asymmetry flat-sector class), while
            // "exact"/"relocked" would implicate the boundary stages.
            return declined_with_selection(
                "self_intersection",
                n,
                Some(format!("sub-patch {i}; candidate={candidate_provenance}")),
                None,
                Some(1e-6),
                selection_fact.into_iter().chain(rung_facts).collect(),
            );
        }
    }
    // PK-283 ring-level crossing gate (validate_all's face-pair class the
    // per-patch check cannot see). Single choke point for every candidate arm.
    if let Some(crossing) = popkern_core::geometry::gregory::binonic_ring_crossing(
        &binonic,
        crate::construct::kp_exact_boundary_binonic_optimizer::RING_CROSSING_TOLERANCE,
    ) {
        let detail = format!(
            "{}; candidate={candidate_provenance}",
            crossing.diagnostic_label()
        );
        if std::env::var_os("POPKERN_KP_TRACE").is_some() {
            eprintln!("[kp_binonic.promotion] ring_crossing decline: {detail}");
        }
        return declined_with_selection(
            "ring_crossing",
            n,
            Some(detail),
            Some(crossing.separation),
            Some(crossing.exemption_band),
            selection_fact.into_iter().chain(rung_facts).collect(),
        );
    }
    const COMMERCIAL_G0_GATE: f64 = RING_G0_GATE;
    let post_g0 = worst_ring_g0_gap(&binonic_gp);
    if post_g0 > COMMERCIAL_G0_GATE {
        return declined_with_selection(
            "g0_exceeded",
            n,
            Some(format!(
                "worst_g0={post_g0:.3e} gate={COMMERCIAL_G0_GATE:.1e}"
            )),
            Some(post_g0),
            Some(COMMERCIAL_G0_GATE),
            selection_fact.into_iter().chain(rung_facts).collect(),
        );
    }

    // Quality guard — revert to two-stripe when binonic is
    // catastrophically worse than the predecessor rung. Threshold
    // 1.2× absorbs floating-point wobble while catching real
    // regressions.
    let post_worst = worst_ring_coplanarity(&binonic_gp);
    if pre_worst > 1e-12 && post_worst > 1.2 * pre_worst {
        return declined_with_selection(
            "quality_regression",
            n,
            Some(format!("pre={pre_worst:.3e} post={post_worst:.3e}")),
            Some(post_worst),
            Some(1.2 * pre_worst),
            selection_fact.into_iter().chain(rung_facts).collect(),
        );
    }

    const COMMERCIAL_GATE: f64 = RING_G1_GATE;
    if post_worst > COMMERCIAL_GATE {
        return declined_with_selection(
            "g1_exceeded",
            n,
            Some(format!(
                "worst_coplanarity={post_worst:.3e} gate={COMMERCIAL_GATE:.1e}"
            )),
            Some(post_worst),
            Some(COMMERCIAL_GATE),
            selection_fact.into_iter().chain(rung_facts).collect(),
        );
    }

    // PK-249 chart honesty: accepted candidates may carry chart
    // handedness-reversal bands the arc-length fold gate cannot see
    // (measured hugging the shared-spoke rails at 1:20–1:50
    // asymmetry). Attempt the transactional orientation repair;
    // on decline, today's candidate serves and the fact records
    // the residual bands.
    let chart_fact = three_quad_kp_candidate_selection::clear_committed_candidate_bands(
        &mut binonic,
        &gate_fillets,
        pre_worst,
        post_worst,
    );

    for (idx, patch) in binonic.into_iter().enumerate() {
        sub_patches[idx] = GregoryPatch::RectangularBiNonic(Box::new(patch));
    }
    KpBinonicPromotionReport {
        relaxed_coplanarity: None,
        diagnostics: selection_fact
            .into_iter()
            .chain(chart_fact)
            .chain(rung_facts)
            .collect(),
    }
}

const PROMOTED_BOUNDARY_RESIDUAL_GATE: f64 = 5.0e-3;
const RING_G0_GATE: f64 = 1.0e-10;
const RING_G1_GATE: f64 = 5.0e-2;

/// Optimizer ladder for one candidate. The second tuple slot is the
/// PK-275 outer-fillet rung outcome (populated only when the rung is
/// armed AND the bi-nonic optimizer produced this candidate — the
/// bi-quintic/C-twist fallbacks never run the rung).
fn optimize_kp_binonic_candidate(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
    outer_fillets: Option<&[popkern_core::geometry::surface::SurfaceGeometry]>,
) -> Option<(Vec<RectangularGregoryBiNonic>, Option<String>)> {
    kp_binonic_optimizer::optimize_binonic_with_outer_fillets(input, options, outer_fillets)
        .map(|result| (result.binonic, result.outer_rung_outcome))
        .or_else(|| {
            optimize_biquintic_and_build_binonic(input, options)
                .map(|result| (result.binonic, None))
        })
        .or_else(|| {
            optimize_c_twists_and_build_binonic(input, options).map(|result| (result.binonic, None))
        })
}

fn declined(
    reason_code: &'static str,
    n_sides: usize,
    detail: Option<String>,
    residual: Option<f64>,
    tolerance: Option<f64>,
) -> KpBinonicPromotionReport {
    declined_with_selection(
        reason_code,
        n_sides,
        detail,
        residual,
        tolerance,
        Vec::new(),
    )
}

/// Decline report carrying the boundary-candidates selection fact
/// (and any PK-275 outer-fillet rung facts) alongside the decline
/// itself, so a declined promotion still records which candidate
/// paths were tried and why each fell short (PK-249 observability).
fn declined_with_selection(
    reason_code: &'static str,
    n_sides: usize,
    detail: Option<String>,
    residual: Option<f64>,
    tolerance: Option<f64>,
    selection_fact: Vec<crate::types::ApexConstructionDiagnostic>,
) -> KpBinonicPromotionReport {
    let mut diagnostics = vec![
        crate::construct::three_quad_kp_diagnostics::promotion_declined(
            reason_code,
            n_sides,
            detail,
            residual,
            tolerance,
        ),
    ];
    diagnostics.extend(selection_fact);
    KpBinonicPromotionReport {
        relaxed_coplanarity: None,
        diagnostics,
    }
}

/// Build a [`DirectBoundaryInput`] from the biseptic-promoted
/// two-stripe seed patches + ring topology. Extracts the per-
/// spoke boundary CPs, corner xtan values, and input-face
/// tangent planes (approximated from the patch normals at C
/// and V — upstream biseptic promotion already wrote these
/// consistent with the true input face planes, so using the
/// patch's computed normal is equivalent).
///
/// `pub(super)` so the sibling `kp_polar_bi6_promote` module
/// can share the same extraction (I.10f).
pub(crate) fn extract_direct_boundary_input(
    patches: &[RectangularGregoryTwoStripe],
    spokes: &[Spoke],
    hex_boundaries: Option<&[([Point3; 10], [Point3; 10])]>,
) -> Option<DirectBoundaryInput> {
    let n = patches.len();
    if n == 0 || spokes.len() != n {
        return None;
    }
    // Central vertex C: patch 0's grid[7][0] (u=1, v=0) per the
    // three-quad convention.
    let central_vertex = patches[0].grid[7][0];

    let mut spoke_data = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    let mut central_vertex_tps = Vec::with_capacity(n);

    for (k, spoke) in spokes.iter().enumerate() {
        let patch_a = &patches[spoke.a.patch_idx];
        let patch_b = &patches[spoke.b.patch_idx];

        // Spoke boundary CPs: for the standard three-quad (A UOne
        // not reversed), these are patch_a.grid[7][0..=7].
        let mut boundary_cps = [Point3::new(0.0, 0.0, 0.0); 8];
        for (j, cp) in boundary_cps.iter_mut().enumerate() {
            *cp = match spoke.a.side {
                Side::UOne => patch_a.grid[7][j],
                Side::UZero => patch_a.grid[0][j],
                Side::VOne => patch_a.grid[j][7],
                Side::VZero => patch_a.grid[j][0],
            };
        }

        // xtan at corner C (spoke-t = 0) and V (spoke-t = 1) on
        // each side. Side A: partial in the cross-boundary
        // direction at the spoke's side. Side B: same, but with
        // reversal mapping spoke-t → B's local parameter via
        // t_local = 1 - t when b.reversed is true.
        let (ua_c, va_c) = uv_at_spoke_t(spoke.a.side, 0.0, spoke.a.reversed);
        let (ua_v, va_v) = uv_at_spoke_t(spoke.a.side, 1.0, spoke.a.reversed);
        let xtan_a_at_c = cross_tangent_on_side(patch_a, spoke.a.side, ua_c, va_c);
        let xtan_a_at_v = cross_tangent_on_side(patch_a, spoke.a.side, ua_v, va_v);

        let (ub_c, vb_c) = uv_at_spoke_t(spoke.b.side, 0.0, spoke.b.reversed);
        let (ub_v, vb_v) = uv_at_spoke_t(spoke.b.side, 1.0, spoke.b.reversed);
        let xtan_b_at_c = cross_tangent_on_side(patch_b, spoke.b.side, ub_c, vb_c);
        let xtan_b_at_v = cross_tangent_on_side(patch_b, spoke.b.side, ub_v, vb_v);

        spoke_data.push(SpokeBoundaryData {
            boundary_cps,
            xtan_a_at_c,
            xtan_b_at_c,
            xtan_a_at_v,
            xtan_b_at_v,
            // Radii not used by the current blend-coeffs formula
            // (constant-α form). Populate with a plausible placeholder
            // so the struct validates; the K-K-C 1991 variable-
            // radius extension will consume these later.
            radius_at_c: 1.0,
            radius_at_v: 1.0,
            edge_id: format!("spoke_{k}"),
        });

        // Outer corner position: the shared (u=1, v=1) of patch A.
        let v_pos = patch_a.grid[7][7];
        // Tangent plane at V from patch A's normal there.
        let du_v = patch_a.partial_u(1.0, 1.0);
        let dv_v = patch_a.partial_v(1.0, 1.0);
        let n_v = du_v.cross(&dv_v).normalized()?;
        let tps_v = TangentPlane::from_normal(n_v)?;
        outer_corners.push(OuterCornerData {
            position: v_pos,
            tangent_plane: tps_v,
        });

        // Central-vertex tangent plane from patch A's normal at C.
        let du_c = patch_a.partial_u(1.0, 0.0);
        let dv_c = patch_a.partial_v(1.0, 0.0);
        let n_c = du_c.cross(&dv_c).normalized()?;
        let tps_c = TangentPlane::from_normal(n_c)?;
        central_vertex_tps.push(tps_c);
    }

    // Phase 6c-D Stage 1d: when the caller provides
    // hex_boundaries (deg-9 Bernstein CPs of the substituted
    // FilletEnd / FaceEdge curves per sub-patch), thread them
    // into the new `HexBoundaryData` field. The K-P bi-nonic
    // ring builder uses these to bypass `k_p_boundary` Hermite
    // synthesis and trace the analytic `sphere ∩ cylinder` arc
    // exactly. `None` falls through to existing behaviour.
    let hex_boundaries_data = hex_boundaries.map(|src| {
        src.iter()
            .map(
                |(row0, col9)| popkern_core::geometry::gregory::HexBoundaryData {
                    row0: *row0,
                    col9: *col9,
                },
            )
            .collect::<Vec<_>>()
    });

    Some(DirectBoundaryInput {
        central_vertex,
        spokes: spoke_data,
        outer_corners,
        central_vertex_tps,
        seed_patches: Some(patches.to_vec()),
        hex_boundaries: hex_boundaries_data,
    })
}

/// Map a spoke parameter `t ∈ [0, 1]` to the `(u, v)` on a
/// side-indexed patch boundary, honouring the `reversed` flag
/// that three-quad uses for its VZero B-ends (spoke-t = 0 ↔
/// patch-u = 1 when reversed).
fn uv_at_spoke_t(side: Side, t: f64, reversed: bool) -> (f64, f64) {
    let t_local = if reversed { 1.0 - t } else { t };
    match side {
        Side::UZero => (0.0, t_local),
        Side::UOne => (1.0, t_local),
        Side::VZero => (t_local, 0.0),
        Side::VOne => (t_local, 1.0),
    }
}

/// Partial derivative in the cross-boundary direction at the
/// given side. For U-sides that's `∂S/∂u`; for V-sides `∂S/∂v`.
fn cross_tangent_on_side(
    patch: &RectangularGregoryTwoStripe,
    side: Side,
    u: f64,
    v: f64,
) -> popkern_core::geometry::Vec3 {
    match side {
        Side::UZero | Side::UOne => patch.partial_u(u, v),
        Side::VZero | Side::VOne => patch.partial_v(u, v),
    }
}

/// Estimate a characteristic length for tolerance scaling — the
/// max bounding-box extent across the seed patches. Populates
/// `DirectConstructionOptions::characteristic_length` so all
/// ratio-based tolerances scale with the geometry's natural
/// size.
fn estimate_characteristic_length(patches: &[RectangularGregoryTwoStripe]) -> f64 {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for patch in patches {
        for i in 0..8 {
            for j in 0..8 {
                let p = patch.grid[i][j];
                for (k, v) in [p.x, p.y, p.z].iter().enumerate() {
                    if *v < min[k] {
                        min[k] = *v;
                    }
                    if *v > max[k] {
                        max[k] = *v;
                    }
                }
            }
        }
    }
    let dx = max[0] - min[0];
    let dy = max[1] - min[1];
    let dz = max[2] - min[2];
    (dx * dx + dy * dy + dz * dz).sqrt().max(1e-9)
}
