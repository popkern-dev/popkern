//! K-P bi-nonic promotion helpers for three-quad decompositions.

use super::*;
use crate::construct::kp_binonic_optimizer;
use crate::construct::kp_biquintic_optimizer::optimize_biquintic_and_build_binonic;
use crate::construct::kp_twist_optimizer::optimize_c_twists_and_build_binonic;

/// Promote sub-patches from `RectangularTwoStripe` to
/// `RectangularBiNonic` via the Karciauskas-Peters guided-surface
/// construction. Returns `Some(worst_coplanarity)` when the
/// binonic ships with relaxed continuity (worst > 5e-2 commercial
/// gate — Phase H.8 step-D graceful degradation), `None` when the
/// promotion is accepted at commercial grade or reverted.
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
) -> Option<f64> {
    let n = sub_patches.len();
    if n == 0 {
        return None;
    }
    // Every sub-patch must already be two-stripe — we read the
    // grid as the source of the K-P DirectBoundaryInput.
    let mut two_stripe: Vec<RectangularGregoryTwoStripe> = Vec::with_capacity(n);
    for patch in sub_patches.iter() {
        let GregoryPatch::RectangularTwoStripe(ts) = patch else {
            trace_kp("non-two-stripe input", n, None);
            return None;
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
        trace_kp("extract_direct_boundary_input → None", n, None);
        return None;
    };
    let char_len = estimate_characteristic_length(&two_stripe);
    let options = DirectConstructionOptions::with_characteristic_length(char_len);
    // Phase H.7 gap #2 closure — K-P 2004 §4 iterative twist
    // optimization. The one-shot Shirman-Séquin twist is under-
    // conditioned for extreme radius asymmetry (e.g. cube_extreme
    // 1:10); this loop refines the N central-vertex twists by
    // Gauss-Newton-style gradient descent against the ring-wide
    // G1 coplanarity residual. Converges in 3-5 iterations
    // typically; always monotonically improves the baseline
    // (the iterate tracks the best-seen residual).
    // Option C step 3: three-phase optimizer (bicubic → bi-quintic
    // → bi-nonic). Phase 3 only engages when phases 1+2 fail to
    // converge; for the vast majority of fixtures, behaviour is
    // identical to the bi-quintic path. Falls through to the bi-
    // quintic-only optimizer, then to the bicubic C-twist
    // optimizer (Option B) if both fail.
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

    let binonic = match kp_binonic_optimizer::optimize_binonic_with_outer_fillets(
        &input,
        &options,
        outer_fillets.as_deref(),
    ) {
        Some(r) => r.binonic,
        None => match optimize_biquintic_and_build_binonic(&input, &options) {
            Some(r) => r.binonic,
            None => {
                let Some(result) = optimize_c_twists_and_build_binonic(&input, &options) else {
                    trace_kp("optimize_binonic_and_build_binonic → None", n, None);
                    return None;
                };
                result.binonic
            }
        },
    };

    // Self-intersection gate — never commit a fold.
    let binonic_gp: Vec<GregoryPatch> = binonic
        .iter()
        .cloned()
        .map(|p| GregoryPatch::RectangularBiNonic(Box::new(p)))
        .collect();
    for (i, gp) in binonic_gp.iter().enumerate() {
        if gp.self_intersects_with_tolerance(1e-6) {
            trace_kp(
                "self_intersects_with_tolerance",
                n,
                Some(format!("sub-patch {i}")),
            );
            return None;
        }
    }

    // Quality guard — revert to two-stripe when binonic is
    // catastrophically worse than the predecessor rung. Threshold
    // 1.2× absorbs floating-point wobble while catching real
    // regressions.
    let pre_worst = worst_ring_coplanarity(sub_patches);
    let post_worst = worst_ring_coplanarity(&binonic_gp);
    // Phase 6c-D Stage 3: under boundary injection, the binonic
    // patch is built from a boundary-consistent frame + Coons
    // interior, NOT from the K-P guide-driven L-M optimum. Spoke
    // G1 is determined entirely by the (analytic) row0/col9 +
    // (input) row9/col0 geometry — the L-M can't refine it. The
    // 1.2× regression guard was tuned for the guide-driven path
    // and rejects legitimate-but-different geometry under
    // injection. Skip the guard when injection is active; the
    // self-intersection gate (already checked above) is the
    // remaining quality floor under injection.
    let injection_active = hex_boundaries.is_some();
    if !injection_active && pre_worst > 1e-12 && post_worst > 1.2 * pre_worst {
        trace_kp(
            "quality regression",
            n,
            Some(format!("pre={pre_worst:.3e} post={post_worst:.3e}")),
        );
        return None;
    }

    for (idx, patch) in binonic.into_iter().enumerate() {
        sub_patches[idx] = GregoryPatch::RectangularBiNonic(Box::new(patch));
    }

    // Phase H.8 step-D — `Some(worst)` signals relaxed continuity.
    const COMMERCIAL_GATE: f64 = 5e-2;
    if post_worst > COMMERCIAL_GATE {
        Some(post_worst)
    } else {
        None
    }
}

/// Worst-case G1 coplanarity across the full ring's shared spokes,
/// using 50 uniform samples on `t ∈ [0.05, 0.95]` per spoke. The
/// endpoint margin matches the commercial-gate test's convention
/// (avoids the rational-Gregory corner-singularity finite-angle
/// transition zone). Returns 0 if normals are ill-conditioned at
/// every sample (can't happen on a real ring, but guards against
/// log-scale division).
pub(crate) fn worst_ring_coplanarity(patches: &[GregoryPatch]) -> f64 {
    let n = patches.len();
    let mut worst = 0.0_f64;
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        for sample in 0..50 {
            let t = 0.05 + 0.90 * (sample as f64) / 49.0;
            let t_b = 1.0 - t;
            let xtan_a = a.partial_u(1.0, t);
            let xtan_b = b.partial_v(t_b, 0.0);
            let t_spoke = a.partial_v(1.0, t);
            let cross = xtan_a.cross(&xtan_b);
            let numer = (cross.x * t_spoke.x + cross.y * t_spoke.y + cross.z * t_spoke.z).abs();
            let na = (xtan_a.x * xtan_a.x + xtan_a.y * xtan_a.y + xtan_a.z * xtan_a.z).sqrt();
            let nb = (xtan_b.x * xtan_b.x + xtan_b.y * xtan_b.y + xtan_b.z * xtan_b.z).sqrt();
            let ns = (t_spoke.x * t_spoke.x + t_spoke.y * t_spoke.y + t_spoke.z * t_spoke.z).sqrt();
            let denom = na * nb * ns;
            if denom > 1e-300 {
                let coplan = numer / denom;
                if coplan > worst {
                    worst = coplan;
                }
            }
        }
    }
    worst
}

/// Emit a diagnostic trace when the K-P promotion silently
/// falls back to two-stripe. Gated on the `POPKERN_KP_TRACE`
/// env var so production builds stay silent; set `=1` before
/// running a fillet test to see which stage rejected the input.
fn trace_kp(stage: &str, n: usize, detail: Option<String>) {
    if std::env::var_os("POPKERN_KP_TRACE").is_some() {
        match detail {
            Some(d) => eprintln!("[kp_binonic N={n}] fallback @ {stage}: {d}"),
            None => eprintln!("[kp_binonic N={n}] fallback @ {stage}"),
        }
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
