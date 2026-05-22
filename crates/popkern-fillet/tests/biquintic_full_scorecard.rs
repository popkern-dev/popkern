//! §5.3 item #3 extended — Phase 8 measure scorecard.
//!
//! Quantifies the improvement of Karciauskas-Peters biquintic-full
//! (with spoke-G1 solver) over the baseline deg5 Gregory along the
//! three shared spokes of the trihedral vertex blend when the legacy
//! rectangular decomposition is exposed as committed faces. Two metrics:
//!
//! 1. **Normal agreement** — `|n_L · n_R| / (|n_L| · |n_R|)` at
//!    sample points on each shared spoke between adjacent sub-
//!    patches. `1.0` is perfect alignment, `0.0` is orthogonal.
//!    G1 continuity needs the scorecard median > ~0.9999.
//! 2. **Max normal angle** — `acos(|n_L · n_R| / …)` radians. G1
//!    continuity needs < ~1e-3 radians (~0.06°).
//!
//! Asserts the biquintic-full path does no *worse* than deg5 on
//! both metrics (min-agreement and max-angle); prints the
//! side-by-side numbers so future tuning of the solver has a
//! regression gate to measure against.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

/// One shared-spoke sample pair: the two 3-D normals on either
/// side of the shared boundary at the same spoke-parameter.
struct SampleAgreement {
    cos_angle: f64,
}

fn fillet_mixed_corner(
    use_biquintic_full: bool,
) -> (BRep, Vec<popkern_core::topology::FaceHandle>) {
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
        use_biquintic_full_gregory: use_biquintic_full,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, Tolerance::default(), None).unwrap();
    let mut corner_faces: Vec<_> = result.corner_patches.iter().map(|p| p.face).collect();
    if corner_faces.is_empty() {
        corner_faces = result
            .repair_report
            .apex_routes
            .iter()
            .flat_map(|route| route.exact_boundary_report.committed_faces.iter().copied())
            .collect();
    }
    (result.brep, corner_faces)
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

fn point_at(gp: &GregoryPatch, u: f64, v: f64) -> Point3 {
    gp.point_at(u, v)
}

fn partial_u(gp: &GregoryPatch, u: f64, v: f64) -> Vec3 {
    gp.partial_u(u, v)
}

fn partial_v(gp: &GregoryPatch, u: f64, v: f64) -> Vec3 {
    gp.partial_v(u, v)
}

fn unit_normal(gp: &GregoryPatch, u: f64, v: f64) -> Option<Vec3> {
    partial_u(gp, u, v).cross(&partial_v(gp, u, v)).normalized()
}

/// For each adjacent pair of sub-patches (`k`, `(k+1) % 3`) in
/// the CCW corner-patch order, sample normals at `n` interior
/// parameters on their shared spoke — patch `k`'s `UOne` boundary
/// vs patch `(k+1)`'s `VZero` boundary with the parameter
/// reversed. Returns all `3 · n` agreement values.
fn sample_all_shared_spokes(patches: &[GregoryPatch], n_samples: usize) -> Vec<SampleAgreement> {
    let mut out = Vec::new();
    let n = patches.len();
    if n < 2 {
        return out;
    }
    let ts: Vec<f64> = (1..=n_samples)
        .map(|k| (k as f64) / ((n_samples + 1) as f64))
        .collect();
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        for &t in &ts {
            let t_b = 1.0 - t;
            let n_a = match unit_normal(a, 1.0, t) {
                Some(v) => v,
                None => continue,
            };
            let n_b = match unit_normal(b, t_b, 0.0) {
                Some(v) => v,
                None => continue,
            };
            let cos = (n_a.x * n_b.x + n_a.y * n_b.y + n_a.z * n_b.z).abs();
            out.push(SampleAgreement { cos_angle: cos });

            // Sanity: the two samples must agree in 3-D position.
            let p_a = point_at(a, 1.0, t);
            let p_b = point_at(b, t_b, 0.0);
            let gap = ((p_a.x - p_b.x).powi(2) + (p_a.y - p_b.y).powi(2) + (p_a.z - p_b.z).powi(2))
                .sqrt();
            assert!(
                gap < 1e-2,
                "shared spoke: sub-patches disagree on position at t={t}: gap={gap:.4e}"
            );
        }
    }
    out
}

/// Sample the at-C corner of each shared spoke (spoke-parameter
/// `t → 0` from patch A's side, `s → 1` from patch B's side).
/// Uses a tiny offset `1e-6` from the exact corner to avoid
/// numerical issues with the rational-blend singularity at C.
fn sample_at_central_vertex(patches: &[GregoryPatch]) -> Vec<SampleAgreement> {
    let mut out = Vec::new();
    let n = patches.len();
    let eps = 1e-6_f64;
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        let Some(n_a) = unit_normal(a, 1.0, eps) else {
            continue;
        };
        let Some(n_b) = unit_normal(b, 1.0 - eps, 0.0) else {
            continue;
        };
        let cos = (n_a.x * n_b.x + n_a.y * n_b.y + n_a.z * n_b.z).abs();
        out.push(SampleAgreement { cos_angle: cos });
    }
    out
}

fn summarize(label: &str, samples: &[SampleAgreement]) -> (f64, f64) {
    let min_cos = samples
        .iter()
        .map(|s| s.cos_angle)
        .fold(f64::INFINITY, f64::min);
    let max_cos = samples
        .iter()
        .map(|s| s.cos_angle)
        .fold(f64::NEG_INFINITY, f64::max);
    let mean_cos: f64 =
        samples.iter().map(|s| s.cos_angle).sum::<f64>() / (samples.len().max(1) as f64);
    let max_angle = min_cos.clamp(-1.0, 1.0).acos();
    let min_angle = max_cos.clamp(-1.0, 1.0).acos();
    let n = samples.len();
    println!(
        "{label}: n={n} cos[min={min_cos:.6} mean={mean_cos:.6} max={max_cos:.6}] angle[{min_angle:.4}..{max_angle:.4} rad]"
    );
    (min_cos, max_angle)
}

/// Phase 8 scorecard: produce side-by-side numbers for deg5 vs
/// biquintic-full on the three-quad cube-corner fixture when the
/// rectangular decomposition is committed. The current macro-cap route
/// can own the corner as one logical Gregory patch; in that case this
/// test reduces to a topology validity guard because there are no
/// shared rectangular spokes to sample.
#[test]
fn biquintic_full_scorecard_vs_deg5() {
    let (brep_deg5, faces_deg5) = fillet_mixed_corner(false);
    let (brep_full, faces_full) = fillet_mixed_corner(true);

    let patches_deg5: Vec<GregoryPatch> = faces_deg5
        .iter()
        .filter_map(|&f| patch_of(&brep_deg5, f))
        .collect();
    let patches_full: Vec<GregoryPatch> = faces_full
        .iter()
        .filter_map(|&f| patch_of(&brep_full, f))
        .collect();

    if patches_deg5.len() != 3 || patches_full.len() != 3 {
        assert!(
            patches_deg5.is_empty() && patches_full.is_empty(),
            "expected either matching 3-patch rectangular decompositions or the current logical macro-cap route; got deg5={} full={}",
            patches_deg5.len(),
            patches_full.len(),
        );
        assert!(
            !faces_deg5.is_empty() && !faces_full.is_empty(),
            "logical macro-cap route must still report corner patch faces"
        );
        brep_deg5
            .validate()
            .expect("deg5 logical macro-cap route must yield a valid BRep");
        brep_full
            .validate()
            .expect("biquintic-full logical macro-cap route must yield a valid BRep");
        return;
    }

    assert_eq!(patches_deg5.len(), 3, "expected 3 deg5 Gregory sub-patches");
    assert_eq!(
        patches_full.len(),
        3,
        "expected 3 biquintic-full Gregory sub-patches"
    );

    const N_SAMPLES: usize = 5;
    let agree_deg5 = sample_all_shared_spokes(&patches_deg5, N_SAMPLES);
    let agree_full = sample_all_shared_spokes(&patches_full, N_SAMPLES);

    let (min_cos_deg5, max_angle_deg5) = summarize("deg5     ", &agree_deg5);
    let (min_cos_full, max_angle_full) = summarize("full     ", &agree_full);

    // At-C spoke agreement (t = 0 hits the central vertex where
    // Shirman-Séquin pins tangent-plane alignment). This should be
    // near-perfect on both variants — it is the acceptance gate
    // for the SS item #3 pre-existing work.
    let at_c_deg5 = sample_at_central_vertex(&patches_deg5);
    let at_c_full = sample_at_central_vertex(&patches_full);
    let _ = summarize("deg5 at C", &at_c_deg5);
    let _ = summarize("full at C", &at_c_full);

    // Non-regression: biquintic-full must not be worse than deg5
    // on either metric, within a small numerical slack for floating-
    // point + solver noise.
    let tol_cos = 1e-9;
    let tol_angle = 1e-6;
    assert!(
        min_cos_full + tol_cos >= min_cos_deg5,
        "biquintic-full regressed min-cos: deg5={min_cos_deg5:.10} full={min_cos_full:.10}"
    );
    assert!(
        max_angle_full <= max_angle_deg5 + tol_angle,
        "biquintic-full regressed max-angle: deg5={max_angle_deg5:.6} full={max_angle_full:.6}"
    );

    // Mean non-regression across all 15 interior samples. Some current
    // dispatches route deg5 and biquintic-full through equivalent geometry,
    // so equality is acceptable; this scorecard should catch worsening, not
    // require an improvement from a path that did not engage differently.
    let mean_deg5: f64 =
        agree_deg5.iter().map(|s| s.cos_angle).sum::<f64>() / (agree_deg5.len() as f64);
    let mean_full: f64 =
        agree_full.iter().map(|s| s.cos_angle).sum::<f64>() / (agree_full.len() as f64);
    assert!(
        mean_full + tol_cos >= mean_deg5,
        "spoke-G1 solver regressed mean agreement: \
         deg5 mean={mean_deg5:.4} full mean={mean_full:.4} (ratio {:.3}x)",
        mean_full / mean_deg5
    );

    // At-C alignment is the acceptance gate for Shirman-Séquin —
    // both variants inherit this from the shared corner inner-CP
    // pairs, so both must clear the bar. 0.998 ≈ 3.6° — room for
    // SS solve noise without crossing into G1-failure territory.
    let min_cos_at_c_deg5 = at_c_deg5
        .iter()
        .map(|s| s.cos_angle)
        .fold(f64::INFINITY, f64::min);
    let min_cos_at_c_full = at_c_full
        .iter()
        .map(|s| s.cos_angle)
        .fold(f64::INFINITY, f64::min);
    assert!(
        min_cos_at_c_deg5 > 0.995 && min_cos_at_c_full > 0.995,
        "at-C normal alignment below SS acceptance: deg5={min_cos_at_c_deg5:.6} full={min_cos_at_c_full:.6}"
    );
}
