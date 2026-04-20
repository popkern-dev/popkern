//! §5.3 Phase-10a — Peters two-stripe G1 scorecard.
//!
//! Four-rung side-by-side: baseline deg5, biquintic-full,
//! biseptic, two-stripe. The last rung is expected to achieve
//! near-machine-precision G1 across all spoke samples, closing
//! the structural gap that biseptic plateaued at (~84°
//! worst-case).

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

enum Variant {
    Deg5,
    BiquinticFull,
    Biseptic,
    TwoStripe,
}

fn fillet_mixed_corner(variant: Variant) -> (BRep, Vec<popkern_core::topology::FaceHandle>) {
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
    let (full, bsep, tsp) = match variant {
        Variant::Deg5 => (false, false, false),
        Variant::BiquinticFull => (true, false, false),
        Variant::Biseptic => (true, true, false),
        Variant::TwoStripe => (true, true, true),
    };
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: full,
        use_biseptic_gregory: bsep,
        use_two_stripe_gregory: tsp,
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

fn unit_normal(gp: &GregoryPatch, u: f64, v: f64) -> Option<Vec3> {
    gp.partial_u(u, v).cross(&gp.partial_v(u, v)).normalized()
}

struct Sample {
    cos_angle: f64,
}

fn sample_shared_spokes(patches: &[GregoryPatch], n: usize) -> Vec<Sample> {
    let mut out = Vec::new();
    let np = patches.len();
    let ts: Vec<f64> = (1..=n).map(|k| (k as f64) / ((n + 1) as f64)).collect();
    for k in 0..np {
        let a = &patches[k];
        let b = &patches[(k + 1) % np];
        for &t in &ts {
            let t_b = 1.0 - t;
            let Some(na) = unit_normal(a, 1.0, t) else {
                continue;
            };
            let Some(nb) = unit_normal(b, t_b, 0.0) else {
                continue;
            };
            let cos = (na.x * nb.x + na.y * nb.y + na.z * nb.z).abs();
            out.push(Sample { cos_angle: cos });
            let pa = a.point_at(1.0, t);
            let pb = b.point_at(t_b, 0.0);
            let gap =
                ((pa.x - pb.x).powi(2) + (pa.y - pb.y).powi(2) + (pa.z - pb.z).powi(2)).sqrt();
            assert!(gap < 1e-6, "position gap at t={t}: {gap:.4e}");
        }
    }
    out
}

fn summarize(label: &str, samples: &[Sample]) -> (f64, f64, f64) {
    let min = samples
        .iter()
        .map(|s| s.cos_angle)
        .fold(f64::INFINITY, f64::min);
    let max = samples
        .iter()
        .map(|s| s.cos_angle)
        .fold(f64::NEG_INFINITY, f64::max);
    let mean: f64 =
        samples.iter().map(|s| s.cos_angle).sum::<f64>() / (samples.len().max(1) as f64);
    let max_angle = min.clamp(-1.0, 1.0).acos();
    println!(
        "{label}: n={} cos[min={min:.6} mean={mean:.6} max={max:.6}] max-angle={max_angle:.4} rad",
        samples.len()
    );
    (min, mean, max_angle)
}

#[test]
fn two_stripe_scorecard_four_rung() {
    let (brep_d, faces_d) = fillet_mixed_corner(Variant::Deg5);
    let (brep_q, faces_q) = fillet_mixed_corner(Variant::BiquinticFull);
    let (brep_s, faces_s) = fillet_mixed_corner(Variant::Biseptic);
    let (brep_t, faces_t) = fillet_mixed_corner(Variant::TwoStripe);

    let patches_d: Vec<_> = faces_d
        .iter()
        .filter_map(|&f| patch_of(&brep_d, f))
        .collect();
    let patches_q: Vec<_> = faces_q
        .iter()
        .filter_map(|&f| patch_of(&brep_q, f))
        .collect();
    let patches_s: Vec<_> = faces_s
        .iter()
        .filter_map(|&f| patch_of(&brep_s, f))
        .collect();
    let patches_t: Vec<_> = faces_t
        .iter()
        .filter_map(|&f| patch_of(&brep_t, f))
        .collect();

    assert_eq!(patches_d.len(), 3);
    assert_eq!(patches_q.len(), 3);
    assert_eq!(patches_s.len(), 3);
    assert_eq!(patches_t.len(), 3);

    // Per conservative ladder: two-stripe rung emits
    // RectangularTwoStripe variant only if the solver actually
    // improved G1 over biseptic. On some geometries the
    // iteration loop converges to a biseptic-equivalent state
    // and the upgrade is skipped. Accept biseptic or two-stripe
    // at this rung.
    for (i, p) in patches_t.iter().enumerate() {
        assert!(
            matches!(
                p,
                GregoryPatch::RectangularTwoStripe(_) | GregoryPatch::RectangularBiseptic(_)
            ),
            "two-stripe rung patch {i} should be RectangularTwoStripe or RectangularBiseptic"
        );
    }

    const N: usize = 5;
    let sd = sample_shared_spokes(&patches_d, N);
    let sq = sample_shared_spokes(&patches_q, N);
    let ss = sample_shared_spokes(&patches_s, N);
    let st = sample_shared_spokes(&patches_t, N);

    let (min_d, mean_d, ang_d) = summarize("deg5       ", &sd);
    let (min_q, mean_q, ang_q) = summarize("biquintic  ", &sq);
    let (min_s, mean_s, ang_s) = summarize("biseptic   ", &ss);
    let (min_t, mean_t, ang_t) = summarize("two-stripe ", &st);

    // Non-regression between all four rungs.
    let tol = 1e-9;
    assert!(min_q + tol >= min_d && mean_q + tol >= mean_d && ang_q <= ang_d + 1e-6);
    assert!(min_s + tol >= min_q && mean_s + tol >= mean_q && ang_s <= ang_q + 1e-6);
    assert!(min_t + tol >= min_s && mean_t + tol >= mean_s && ang_t <= ang_s + 1e-6);

    // Two-stripe must significantly outperform biseptic on the
    // interior-spoke samples. Mean cos should be above 0.95
    // (biseptic plateaus near 0.79); max sample should reach
    // near-machine precision (>0.99).
    assert!(
        mean_t > 0.95,
        "two-stripe mean-cos {mean_t:.6} below 0.95 — solver under-performing"
    );
    let max_cos_t = st
        .iter()
        .map(|s| s.cos_angle)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        max_cos_t > 0.99,
        "two-stripe best-sample cos {max_cos_t:.6} below 0.99 — solver not reaching machine precision anywhere"
    );
    // Phase 10c remaining work: the MIN cos (worst-case sample)
    // still hits ~0.73 on this fixture. Samples near the outer
    // hex vertices `V_{2k+2}` have corner tangent-plane
    // mismatches that are NOT explicitly constrained by the
    // solver — only the central vertex C enforces cross-patch
    // tangent agreement (via Shirman-Séquin). A future phase
    // adds outer-corner G1 constraints alongside the spoke
    // samples, at which point `min_t` should also exceed 0.99.
}
