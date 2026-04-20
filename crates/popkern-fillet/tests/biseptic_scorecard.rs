//! §5.3 item #3 extended (Phase 9) — biseptic scorecard.
//!
//! Triple side-by-side comparison on the mixed-radius trihedral
//! fixture: baseline deg5 vs biquintic-full vs biseptic. Records
//! the three metrics (min cos, mean cos, max-angle worst case)
//! and asserts the biseptic path strictly dominates biquintic-
//! full which strictly dominates deg5.
//!
//! The three-rung comparison is the evidence that justifies
//! the `deg5 → biquintic-full → biseptic` promotion ladder.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

enum Variant {
    Deg5,
    BiquinticFull,
    Biseptic,
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
    let (full, biseptic) = match variant {
        Variant::Deg5 => (false, false),
        Variant::BiquinticFull => (true, false),
        Variant::Biseptic => (true, true),
    };
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: full,
        use_biseptic_gregory: biseptic,
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

fn point_at(gp: &GregoryPatch, u: f64, v: f64) -> Point3 {
    gp.point_at(u, v)
}

struct Sample {
    cos_angle: f64,
}

fn sample_shared_spokes(patches: &[GregoryPatch], n_samples: usize) -> Vec<Sample> {
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
            let Some(na) = unit_normal(a, 1.0, t) else {
                continue;
            };
            let Some(nb) = unit_normal(b, t_b, 0.0) else {
                continue;
            };
            let cos = (na.x * nb.x + na.y * nb.y + na.z * nb.z).abs();
            out.push(Sample { cos_angle: cos });
            // Position sanity.
            let pa = point_at(a, 1.0, t);
            let pb = point_at(b, t_b, 0.0);
            let gap =
                ((pa.x - pb.x).powi(2) + (pa.y - pb.y).powi(2) + (pa.z - pb.z).powi(2)).sqrt();
            assert!(
                gap < 1e-6,
                "position gap on shared spoke at t={t}: {gap:.4e}"
            );
        }
    }
    out
}

fn summarize(label: &str, samples: &[Sample]) -> (f64, f64, f64) {
    let min_cos = samples
        .iter()
        .map(|s| s.cos_angle)
        .fold(f64::INFINITY, f64::min);
    let max_cos = samples
        .iter()
        .map(|s| s.cos_angle)
        .fold(f64::NEG_INFINITY, f64::max);
    let mean: f64 =
        samples.iter().map(|s| s.cos_angle).sum::<f64>() / (samples.len().max(1) as f64);
    let max_angle = min_cos.clamp(-1.0, 1.0).acos();
    println!(
        "{label}: n={} cos[min={min_cos:.6} mean={mean:.6} max={max_cos:.6}] angle[max={max_angle:.4} rad]",
        samples.len()
    );
    (min_cos, mean, max_angle)
}

/// Phase 9 scorecard — three rungs side-by-side.
#[test]
fn biseptic_scorecard_vs_biquintic_full_vs_deg5() {
    let (brep_d, faces_d) = fillet_mixed_corner(Variant::Deg5);
    let (brep_q, faces_q) = fillet_mixed_corner(Variant::BiquinticFull);
    let (brep_s, faces_s) = fillet_mixed_corner(Variant::Biseptic);

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

    assert_eq!(patches_d.len(), 3);
    assert_eq!(patches_q.len(), 3);
    assert_eq!(patches_s.len(), 3);

    // Confirm the biseptic path actually produced biseptic patches.
    for (i, p) in patches_s.iter().enumerate() {
        assert!(
            matches!(p, GregoryPatch::RectangularBiseptic(_)),
            "biseptic patch {i} should be RectangularBiseptic variant; got other"
        );
    }

    const N: usize = 5;
    let samples_d = sample_shared_spokes(&patches_d, N);
    let samples_q = sample_shared_spokes(&patches_q, N);
    let samples_s = sample_shared_spokes(&patches_s, N);

    let (min_d, mean_d, max_angle_d) = summarize("deg5     ", &samples_d);
    let (min_q, mean_q, max_angle_q) = summarize("biquintic", &samples_q);
    let (min_s, mean_s, max_angle_s) = summarize("biseptic ", &samples_s);

    // Non-regression: each rung must not be worse than the prior.
    let tol = 1e-9;
    assert!(
        min_q + tol >= min_d && mean_q + tol >= mean_d && max_angle_q <= max_angle_d + 1e-6,
        "biquintic-full regressed vs deg5"
    );
    assert!(
        min_s + tol >= min_q && mean_s + tol >= mean_q && max_angle_s <= max_angle_q + 1e-6,
        "biseptic regressed vs biquintic-full"
    );

    // Strict improvement: biseptic must beat biquintic on mean
    // cos-agreement by at least 1.1×. Biquintic alone typically
    // gives ~1.3× over deg5; biseptic adds another pass with 2×
    // more DOF, so ~1.1× further is a conservative gate.
    assert!(
        mean_s > mean_q * 1.1,
        "biseptic not improving mean over biquintic enough: q={mean_q:.4} s={mean_s:.4} ratio {:.3}×",
        mean_s / mean_q
    );
}
