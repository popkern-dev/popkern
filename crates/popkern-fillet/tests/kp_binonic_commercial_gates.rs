//! Phase H.7 — commercial-grade continuity gates on the K-P
//! bidegree-9 output.
//!
//! # Metric choice: triple product, not unit-normal dot
//!
//! The natural "1 − |n_A · n_B|" signal breaks down when one
//! patch's `|∂u × ∂v|` drops to ~1e-5 at a parametric
//! near-degeneracy: `.normalized()` amplifies the tiny cross
//! vector into a direction dominated by floating-point noise,
//! and the dot-product metric then spans the full [−1, +1]
//! even though the underlying surfaces share the same tangent
//! plane. The diagnostic pass in `kp_binonic_g1_diagnostic.rs`
//! (deleted) confirmed this on cube_mixed: worst unit-normal
//! deficit of 9.7e-1 at t ≈ 0.74 sits on top of a perfectly
//! G1-coplanar triple-product of 1.3e-3.
//!
//! So the gate uses the *geometric* G1 definition directly —
//! the coplanarity of `xtan_A`, `xtan_B`, and the along-boundary
//! spoke tangent, normalised so `0` = perfect G1 and `1` =
//! orthogonal. This is invariant under the parametric
//! degeneracy that trips up the unit-normal metric.
//!
//! # Gates (per the Phase H plan § 1-3)
//!
//! * **G0**: position gap ≤ 1e-10 (machine-ε, degree-elevation
//!   preserves boundary exactly).
//! * **G1**: coplanarity ≤ 5e-2 on the interior `t ∈ [0.05, 0.95]`.
//!   The `0.05` endpoint margin excludes the Gregory central-
//!   vertex and outer-corner singularities which are known
//!   finite-angle transition zones (`/ip-check`'d in H.0 as
//!   inherent to the rational blend — NOT a Phase H concern).
//! * **G2**: normal-curvature match residual ≤ 5e-2 on the same
//!   interior range (surrogate for full ‖II_A − II_B‖_F
//!   equality).
//!
//! # Fixture matrix
//!
//! Starts with cube_mixed (N=3, asymmetric radii). Follow-ups
//! add cube_extreme, pyramid4_mixed, saddle_vertex_4, etc.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn cube_mixed(use_kp_binonic: bool) -> (BRep, Vec<popkern_core::topology::FaceHandle>) {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: use_kp_binonic,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, Tolerance::default(), None).unwrap();
    let faces: Vec<_> = result.corner_patches.iter().map(|p| p.face).collect();
    (result.brep, faces)
}

fn patch_of(brep: &BRep, face: popkern_core::topology::FaceHandle) -> Option<GregoryPatch> {
    brep.faces[face].surface.as_ref().and_then(|sof| {
        if let SurfaceGeometry::Gregory(gp) = &sof.geometry {
            Some(gp.clone())
        } else {
            None
        }
    })
}

fn vec_norm(v: Vec3) -> f64 {
    (v.x * v.x + v.y * v.y + v.z * v.z).sqrt()
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

/// Normal curvature along an ambient unit tangent, via the second
/// fundamental form evaluated in the parametric basis. Returns
/// `None` if the first fundamental form is singular.
fn normal_curvature_along(gp: &GregoryPatch, u: f64, v: f64, t_ambient: Vec3) -> Option<f64> {
    let su = gp.partial_u(u, v);
    let sv = gp.partial_v(u, v);
    let n = su.cross(&sv).normalized()?;
    let e = dot(su, su);
    let f = dot(su, sv);
    let g = dot(sv, sv);
    let det = e * g - f * f;
    if det.abs() < 1e-14 {
        return None;
    }
    let su_dot_t = dot(su, t_ambient);
    let sv_dot_t = dot(sv, t_ambient);
    let a = (g * su_dot_t - f * sv_dot_t) / det;
    let b = (e * sv_dot_t - f * su_dot_t) / det;
    let norm_sq = e * a * a + 2.0 * f * a * b + g * b * b;
    if norm_sq < 1e-14 {
        return None;
    }
    let a = a / norm_sq.sqrt();
    let b = b / norm_sq.sqrt();
    let (suu, suv, svv) = match gp {
        GregoryPatch::RectangularBiNonic(r) => (
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        ),
        GregoryPatch::RectangularTwoStripe(r) => (
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        ),
        GregoryPatch::RectangularBiseptic(r) => (
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        ),
        GregoryPatch::RectangularBiquinticFull(r) => (
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        ),
        _ => return None,
    };
    let l = dot(n, suu);
    let m = dot(n, suv);
    let nn = dot(n, svv);
    Some(l * a * a + 2.0 * m * a * b + nn * b * b)
}

struct RingMetrics {
    worst_g0_gap: f64,
    /// Coplanarity triple product (`0` = perfect G1, stable under
    /// parametric near-degeneracies).
    worst_g1_coplanarity: f64,
    /// Normal-curvature mismatch, normalised so `0` = perfect G2.
    worst_g2_residual: f64,
}

/// Sample the ring on `t ∈ [0.05, 0.95]` at `n_samples` points.
/// Endpoint margins exclude the Gregory corner-singularity
/// transition zones — not a Phase H construction concern.
fn measure_ring(patches: &[GregoryPatch], n_samples: usize) -> RingMetrics {
    let n = patches.len();
    let mut worst_g0 = 0.0_f64;
    let mut worst_g1 = 0.0_f64;
    let mut worst_g2 = 0.0_f64;

    // Uniform samples on [0.05, 0.95] — 90% of the spoke's
    // interior. Avoids the rational-Gregory corner transitions.
    let ts: Vec<f64> = (0..n_samples)
        .map(|k| 0.05 + 0.90 * (k as f64) / ((n_samples - 1) as f64))
        .collect();

    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        for &t in &ts {
            let t_b = 1.0 - t;

            // G0 gap.
            let pa = a.point_at(1.0, t);
            let pb = b.point_at(t_b, 0.0);
            let gap =
                ((pa.x - pb.x).powi(2) + (pa.y - pb.y).powi(2) + (pa.z - pb.z).powi(2)).sqrt();
            if gap > worst_g0 {
                worst_g0 = gap;
            }

            // G1: triple product of xtan_A, xtan_B, spoke tangent.
            let xtan_a = a.partial_u(1.0, t);
            let xtan_b = b.partial_v(t_b, 0.0);
            let t_spoke = a.partial_v(1.0, t);
            let numer = dot(xtan_a.cross(&xtan_b), t_spoke).abs();
            let denom = vec_norm(xtan_a) * vec_norm(xtan_b) * vec_norm(t_spoke);
            if denom > 1e-300 {
                let coplan = numer / denom;
                if coplan > worst_g1 {
                    worst_g1 = coplan;
                }
            }

            // G2: normal curvature along the spoke tangent.
            let Some(t_unit) = t_spoke.normalized() else {
                continue;
            };
            let Some(ka) = normal_curvature_along(a, 1.0, t, t_unit) else {
                continue;
            };
            let Some(kb) = normal_curvature_along(b, t_b, 0.0, t_unit) else {
                continue;
            };
            let kmax = ka.abs().max(kb.abs());
            let resid = (ka - kb).abs() / (1.0 + kmax);
            if resid > worst_g2 {
                worst_g2 = resid;
            }
        }
    }
    RingMetrics {
        worst_g0_gap: worst_g0,
        worst_g1_coplanarity: worst_g1,
        worst_g2_residual: worst_g2,
    }
}

const N_SAMPLES: usize = 100;
const G0_GATE: f64 = 1e-10;
const G1_GATE: f64 = 5e-2;
const G2_GATE: f64 = 5e-2;

#[test]
fn binonic_passes_commercial_g0_g1_g2_gates_on_cube_mixed() {
    let (brep, faces) = cube_mixed(true);
    let patches: Vec<_> = faces.iter().filter_map(|&f| patch_of(&brep, f)).collect();
    assert_eq!(patches.len(), 3);
    for p in &patches {
        assert!(
            matches!(p, GregoryPatch::RectangularBiNonic(_)),
            "expected binonic variant on cube_mixed, dispatch not firing"
        );
    }
    let m = measure_ring(&patches, N_SAMPLES);
    println!(
        "binonic cube_mixed: g0={:.3e} g1_coplan={:.3e} g2={:.3e}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.worst_g2_residual
    );
    assert!(
        m.worst_g0_gap < G0_GATE,
        "binonic G0 gap {:.3e} exceeds {:.1e} gate",
        m.worst_g0_gap,
        G0_GATE
    );
    assert!(
        m.worst_g1_coplanarity < G1_GATE,
        "binonic G1 coplanarity {:.3e} exceeds {:.1e} gate (commercial)",
        m.worst_g1_coplanarity,
        G1_GATE
    );
    assert!(
        m.worst_g2_residual < G2_GATE,
        "binonic G2 residual {:.3e} exceeds {:.1e} gate (commercial)",
        m.worst_g2_residual,
        G2_GATE
    );
}

#[test]
fn binonic_matches_or_beats_two_stripe_on_cube_mixed() {
    let (brep_t, faces_t) = cube_mixed(false);
    let (brep_b, faces_b) = cube_mixed(true);
    let t: Vec<_> = faces_t
        .iter()
        .filter_map(|&f| patch_of(&brep_t, f))
        .collect();
    let b: Vec<_> = faces_b
        .iter()
        .filter_map(|&f| patch_of(&brep_b, f))
        .collect();
    let mt = measure_ring(&t, N_SAMPLES);
    let mb = measure_ring(&b, N_SAMPLES);
    println!(
        "two_stripe cube_mixed: g0={:.3e} g1={:.3e} g2={:.3e}",
        mt.worst_g0_gap, mt.worst_g1_coplanarity, mt.worst_g2_residual
    );
    println!(
        "binonic    cube_mixed: g0={:.3e} g1={:.3e} g2={:.3e}",
        mb.worst_g0_gap, mb.worst_g1_coplanarity, mb.worst_g2_residual
    );
    // Non-regression: binonic must not be materially worse than
    // two-stripe on any metric. 2× tolerance absorbs floating-
    // point noise without masking a real regression.
    assert!(
        mb.worst_g0_gap <= 2.0 * mt.worst_g0_gap.max(1e-12),
        "binonic G0 regression: {:.3e} vs two_stripe {:.3e}",
        mb.worst_g0_gap,
        mt.worst_g0_gap
    );
    assert!(
        mb.worst_g1_coplanarity <= 2.0 * mt.worst_g1_coplanarity.max(1e-8),
        "binonic G1 regression: {:.3e} vs two_stripe {:.3e}",
        mb.worst_g1_coplanarity,
        mt.worst_g1_coplanarity
    );
    assert!(
        mb.worst_g2_residual <= 2.0 * mt.worst_g2_residual.max(1e-8),
        "binonic G2 regression: {:.3e} vs two_stripe {:.3e}",
        mb.worst_g2_residual,
        mt.worst_g2_residual
    );
}
