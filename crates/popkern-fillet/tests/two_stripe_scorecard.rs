//! §5.3 Phase-10a — Peters two-stripe G1 scorecard.
//!
//! Rung ladder side-by-side on the fixture the rectangular sub-patch
//! ladder actually reaches: the uniform square-pyramid apex
//! (r = 0.05), whose uniform promotion arm stands down on explicit
//! legacy-family flags (PK-229 `explicit_legacy_surface_family`
//! carve-out).
//!
//! PK-230: migrated off the mixed-radius cube corners. cube_mixed has
//! been setback-served since PK-226 (the old scorecard passed through
//! a macro-cap escape hatch with zero flag coverage), and cube_extreme
//! reached the ladder only while the pre-PK-230 grazing false
//! positive made the serving retry re-route its promoted direct-apex
//! serve back to legacy — the honest validator retired that
//! accident. On mixed corners the deg5-only baseline IS the default
//! flag set, and defaults promote mixed N∈{3,4} apexes (PK-099).
//!
//! Measured rungs (validation sweep 2026-08-11, min/mean cos over
//! shared-spoke samples):
//!
//! | rung        | min    | mean   |
//! |-------------|--------|--------|
//! | deg5        | (promoted direct-apex serve — no sub-patches)   |
//! | biquintic   | 0.9424 | 0.9752 |
//! | biseptic    | 0.9424 | 0.9752 | (coupled-ladder gate: see below)
//! | two-stripe  | 0.9999 | 1.0000 |
//!
//! This pin must FAIL, not skip, if the ladder stops committing
//! rectangular sub-patches on this fixture (PK-230 acceptance: no
//! vacuous flag pins).

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{BlendResult, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

enum Variant {
    Deg5,
    BiquinticFull,
    Biseptic,
    TwoStripe,
}

fn fillet_pyramid4_apex(variant: Variant) -> BlendResult {
    let brep = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let incident_edges = |vertex: popkern_core::topology::VertexHandle| {
        brep.half_edges
            .values()
            .filter(|he| he.origin == vertex)
            .map(|he| he.edge)
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>()
    };
    let vertex = brep
        .vertices
        .keys()
        .find(|&v| incident_edges(v).len() == 4)
        .expect("square_pyramid valence-4 apex");
    let edges = incident_edges(vertex);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.05),
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
    fillet_edges(&brep, &specs, &opts, Tolerance::default(), None).unwrap()
}

fn corner_patches(result: &BlendResult) -> Vec<GregoryPatch> {
    result
        .corner_patches
        .iter()
        .filter_map(|p| result.brep.faces.get(p.face))
        .filter_map(|f| f.surface.as_ref())
        .filter_map(|sof| {
            if let SurfaceGeometry::Gregory(gp) = &sof.geometry {
                Some(gp.clone())
            } else {
                None
            }
        })
        .collect()
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
fn two_stripe_scorecard_ladder_rungs() {
    let result_d = fillet_pyramid4_apex(Variant::Deg5);
    let result_q = fillet_pyramid4_apex(Variant::BiquinticFull);
    let result_s = fillet_pyramid4_apex(Variant::Biseptic);
    let result_t = fillet_pyramid4_apex(Variant::TwoStripe);

    // The deg5-only set IS the default configuration: the uniform
    // apex promotes to the direct-apex route and no rectangular
    // sub-patches ship. This is the documented reaching condition —
    // the LADDER rungs below are the flag-effect claims.
    let patches_d = corner_patches(&result_d);
    assert!(
        patches_d.is_empty(),
        "deg5-only (default) flag set should serve the promoted direct-apex route \
         with no rectangular Gregory sub-patches; got {}",
        patches_d.len()
    );

    let patches_q = corner_patches(&result_q);
    let patches_s = corner_patches(&result_s);
    let patches_t = corner_patches(&result_t);

    // The explicit-family rungs must COMMIT rectangular sub-patches —
    // an empty set means the fixture no longer reaches the ladder and
    // the scorecard must migrate (PK-230: no vacuous flag pins).
    assert_eq!(
        patches_q.len(),
        4,
        "biquintic rung must commit 4 sub-patches"
    );
    assert_eq!(
        patches_s.len(),
        4,
        "biseptic rung must commit 4 sub-patches"
    );
    assert_eq!(
        patches_t.len(),
        4,
        "two-stripe rung must commit 4 sub-patches"
    );

    for (i, p) in patches_q.iter().enumerate() {
        assert!(
            matches!(p, GregoryPatch::RectangularBiquinticFull(_)),
            "biquintic rung patch {i} should be RectangularBiquinticFull"
        );
    }
    // Biseptic and two-stripe are a coupled ladder unit: with
    // two-stripe off, run_promotion_ladder short-circuits at the
    // pre-two-stripe gate ("pre_two_stripe_disabled") BEFORE biseptic
    // promotion, so this rung ships the biquintic state by design.
    // Biseptic-typed patches are an intermediate the two-stripe rung
    // consumes; they surface only if the two-stripe promote declines.
    for (i, p) in patches_s.iter().enumerate() {
        assert!(
            matches!(
                p,
                GregoryPatch::RectangularBiseptic(_) | GregoryPatch::RectangularBiquinticFull(_)
            ),
            "biseptic rung patch {i} should be RectangularBiseptic or RectangularBiquinticFull"
        );
    }
    // Two-stripe measurably improves on this fixture and must ship
    // its own variant (flag effect — the PK-230 acceptance).
    for (i, p) in patches_t.iter().enumerate() {
        assert!(
            matches!(p, GregoryPatch::RectangularTwoStripe(_)),
            "two-stripe rung patch {i} should be RectangularTwoStripe"
        );
    }

    const N: usize = 5;
    let sq = sample_shared_spokes(&patches_q, N);
    let ss = sample_shared_spokes(&patches_s, N);
    let st = sample_shared_spokes(&patches_t, N);

    let (min_q, mean_q, ang_q) = summarize("biquintic  ", &sq);
    let (min_s, mean_s, ang_s) = summarize("biseptic   ", &ss);
    let (min_t, mean_t, ang_t) = summarize("two-stripe ", &st);

    // Non-regression between the explicit-family rungs.
    let tol = 1e-9;
    assert!(min_s + tol >= min_q && mean_s + tol >= mean_q && ang_s <= ang_q + 1e-6);
    assert!(min_t + tol >= min_s && mean_t + tol >= mean_s && ang_t <= ang_s + 1e-6);

    // Quality at the measured class (margins below the 2026-08-11
    // measurements: biquintic 0.9424/0.9752, two-stripe
    // 0.9999/1.0000).
    assert!(
        mean_q > 0.9,
        "biquintic mean-cos {mean_q:.6} below measured class (0.9752)"
    );
    assert!(
        min_t > 0.99 && mean_t > 0.995,
        "two-stripe min/mean-cos {min_t:.6}/{mean_t:.6} below measured class \
         (0.9999/1.0000)"
    );
    // Two-stripe must measurably outperform biquintic — the
    // end-to-end justification for the upper ladder.
    assert!(
        min_t > min_q + 0.01,
        "two-stripe min-cos {min_t:.6} not measurably above biquintic's {min_q:.6}"
    );
}
