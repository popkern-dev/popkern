//! PK-192 / PK-193 — PnS cap coverage at the low valences.
//!
//! The PnS bi-4 cap was dispatched only at `n ≥ 5` for most of its
//! life. `n = 4` was believed structurally impossible (σ's
//! `1/cos(2π/n)` pole) and `n = 3` was never wired up. Both beliefs
//! were wrong, and both valences are the *common* ones in mechanical
//! parts: `n = 3` is the cube corner, `n = 4` the pyramid apex.
//!
//! These gates hold the low-valence caps to the same live-seam
//! continuity bar as the rest of the family. Continuity is measured
//! across the committed seam topology with finite differences in UV,
//! so it is independent of which surface kind the cap shipped.
//!
//! The A/B against the incumbent route is
//! `pns_low_valence_route_comparison`, which reads the PnS dispatch
//! switch from the environment and is therefore ignored by default —
//! see its doc comment for the run command.

mod kp_gate_helpers;

use helpers::{
    G0_GATE, G1_GATE, N_SAMPLES, incident_edges, measure_live_cap_seams, pick_vertex_with_valence,
};
use kp_gate_helpers as helpers;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::types::{CornerMethod, FilletOptions, FilletRadius, FilletSpec};
use popkern_fillet::{ContinuityLevel, fillet_edges};

/// What a low-valence apex fillet produced: the continuity metrics
/// plus the construction route that actually served the corner.
struct ApexOutcome {
    metrics: helpers::LiveSeamMetrics,
    /// Distinct `CornerMethod`s across the committed corner patches.
    methods: Vec<CornerMethod>,
    /// Distinct achieved continuity levels across those patches.
    achieved: Vec<ContinuityLevel>,
    cap_faces: usize,
    valid: bool,
    /// Structured `pns.*` diagnostics the pipeline recorded for this
    /// corner — the record of why PnS did or did not serve it.
    pns_diagnostics: Vec<String>,
    /// Per-seam breakdown: which surface kinds meet, and how badly.
    seam_detail: Vec<String>,
}

impl ApexOutcome {
    fn served_by_pns(&self) -> bool {
        self.methods.contains(&CornerMethod::PolyhedralNetSpline)
    }
}

/// Fillet an n-gonal pyramid apex and measure the committed cap
/// across its live seam topology, recording which route served it.
///
/// Route identity is read from the pipeline's own
/// `CornerPatchInfo::method`, not re-derived here — a gate that
/// guessed the route could pass while measuring a cap the route under
/// test never built.
fn run_pyramid_apex_outcome(n: usize, radii: &[f64]) -> ApexOutcome {
    assert_eq!(radii.len(), n);
    let pyr = BRep::n_gonal_pyramid(n, 1.0, 3.0, Tolerance::default())
        .expect("n_gonal_pyramid fixture must build");
    let vertex = pick_vertex_with_valence(&pyr, n).expect("n_gonal_pyramid apex");
    let edges = incident_edges(&pyr, vertex);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        use_kp_polar_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&pyr, &specs, &opts, Tolerance::default(), None)
        .expect("low-valence apex fillet must not error");
    let mut methods: Vec<CornerMethod> = Vec::new();
    let mut achieved: Vec<ContinuityLevel> = Vec::new();
    for patch in &result.corner_patches {
        if !methods.contains(&patch.method) {
            methods.push(patch.method);
        }
        if !achieved.contains(&patch.achieved_continuity) {
            achieved.push(patch.achieved_continuity);
        }
    }
    let faces = helpers::cap_outer_faces_from_result(&result);
    let valid = result.brep.validate().is_ok();
    let metrics = measure_live_cap_seams(&result.brep, &faces, N_SAMPLES);
    let pns_diagnostics: Vec<String> = result
        .repair_report
        .apex_diagnostics
        .iter()
        .filter_map(|d| match d {
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(e)
                if e.stage.starts_with("pns.") =>
            {
                Some(format!(
                    "{}[{}]{}",
                    e.stage,
                    e.code,
                    e.detail_code
                        .as_ref()
                        .map(|d| format!(" {d}"))
                        .unwrap_or_default()
                ))
            }
            _ => None,
        })
        .collect();
    let seam_detail = describe_cap_seams(&result.brep, &faces);
    ApexOutcome {
        metrics,
        methods,
        achieved,
        cap_faces: faces.len(),
        valid,
        pns_diagnostics,
        seam_detail,
    }
}

/// Name each measured cap-cap seam by the surface kinds it joins, so
/// a bad number can be attributed to a real pair of faces rather than
/// blamed on the measurement.
fn describe_cap_seams(
    brep: &popkern_core::topology::BRep,
    faces: &[popkern_core::topology::FaceHandle],
) -> Vec<String> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let kind = |f: popkern_core::topology::FaceHandle| -> String {
        match brep.faces.get(f).and_then(|d| d.surface.as_ref()) {
            Some(s) => match &s.geometry {
                SurfaceGeometry::Plane { .. } => "Plane",
                SurfaceGeometry::Cylinder { .. } => "Cylinder",
                SurfaceGeometry::Cone { .. } => "Cone",
                SurfaceGeometry::Sphere { .. } => "Sphere",
                SurfaceGeometry::Torus { .. } => "Torus",
                SurfaceGeometry::Nurbs(_) => "Nurbs",
                SurfaceGeometry::Gregory(_) => "Gregory",
                SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline",
            }
            .to_string(),
            None => "none".to_string(),
        }
    };
    let set: std::collections::BTreeSet<_> = faces.iter().copied().collect();
    let mut out = Vec::new();
    for (_eh, edge) in brep.edges.iter() {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let (Some(la), Some(lb)) = (brep.loops.get(he.loop_), brep.loops.get(twin.loop_)) else {
            continue;
        };
        let (fa, fb) = (la.face, lb.face);
        if fa == fb || !set.contains(&fa) || !set.contains(&fb) {
            continue;
        }
        let m = measure_live_cap_seams(brep, &[fa, fb], N_SAMPLES);
        let profile =
            helpers::seam_cos_profile(brep, fa, fb, &[0.05, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95])
                .map(|p| {
                    p.iter()
                        .map(|(f, c)| format!("{f:.2}:{c:.4}"))
                        .collect::<Vec<_>>()
                        .join(" ")
                })
                .unwrap_or_else(|| "<no profile>".to_string());
        out.push(format!(
            "{} x {} -> g0={:.2e} g1={:.2e} cos_min={:.5}  cos(f)= {profile}",
            kind(fa),
            kind(fb),
            m.worst_g0_gap,
            m.worst_g1_normal,
            m.cos_min
        ));
    }
    out
}

fn run_pyramid_apex(n: usize, radii: &[f64]) -> helpers::LiveSeamMetrics {
    let outcome = run_pyramid_apex_outcome(n, radii);
    assert!(outcome.valid, "n={n} apex fillet must leave a valid BRep");
    outcome.metrics
}

/// `n = 3` — the cube/tetra corner. Previously below the dispatch
/// floor despite the adapter supporting it all along.
#[test]
#[ignore = "PK-192 target gate: PnS does not yet serve n=3 corners — they are built by the exact-boundary Coons route, which is G0-only across its internal seams (PK-197). Flip when PnS dispatch reaches this valence."]
fn pns_serves_valence_3_apex_within_gates() {
    let m = run_pyramid_apex(3, &[0.05; 3]);
    assert!(
        m.seam_count > 0,
        "n=3 cap must commit at least one live seam to measure"
    );
    assert!(
        m.worst_g0_gap < G0_GATE,
        "n=3 cap worst G0 seam gap {:.3e} exceeds {:.0e}",
        m.worst_g0_gap,
        G0_GATE
    );
    assert!(
        m.worst_g1_normal < G1_GATE,
        "n=3 cap worst G1 normal misalignment {:.3e} exceeds {:.0e}",
        m.worst_g1_normal,
        G1_GATE
    );
}

/// `n = 4` — the pyramid apex and cross junction. This is the
/// valence the kernel refused outright until the σ pole was shown
/// removable.
#[test]
#[ignore = "PK-192 target gate: PnS does not yet serve n=4 corners — they are built by the exact-boundary Coons route, which is G0-only across its internal seams (PK-197). Flip when PnS dispatch reaches this valence."]
fn pns_serves_valence_4_apex_within_gates() {
    let m = run_pyramid_apex(4, &[0.05; 4]);
    assert!(
        m.seam_count > 0,
        "n=4 cap must commit at least one live seam to measure"
    );
    assert!(
        m.worst_g0_gap < G0_GATE,
        "n=4 cap worst G0 seam gap {:.3e} exceeds {:.0e}",
        m.worst_g0_gap,
        G0_GATE
    );
    assert!(
        m.worst_g1_normal < G1_GATE,
        "n=4 cap worst G1 normal misalignment {:.3e} exceeds {:.0e}",
        m.worst_g1_normal,
        G1_GATE
    );
}

/// Mixed radii at the low valences — the asymmetric case that
/// distinguishes a genuine cap construction from one that only works
/// on a symmetric fixture.
#[test]
#[ignore = "PK-192 target gate: PnS does not yet serve low-valence mixed-radius corners — they are built by the exact-boundary Coons route, which is G0-only across its internal seams (PK-197). Flip when PnS dispatch reaches this valence."]
fn pns_serves_low_valence_mixed_radii() {
    for (n, radii) in [
        (3_usize, vec![0.04, 0.06, 0.05]),
        (4, vec![0.04, 0.06, 0.05, 0.055]),
    ] {
        let m = run_pyramid_apex(n, &radii);
        assert!(
            m.worst_g0_gap < G0_GATE,
            "n={n} mixed-radius cap worst G0 {:.3e} exceeds {:.0e}",
            m.worst_g0_gap,
            G0_GATE
        );
        assert!(
            m.worst_g1_normal < G1_GATE,
            "n={n} mixed-radius cap worst G1 {:.3e} exceeds {:.0e}",
            m.worst_g1_normal,
            G1_GATE
        );
    }
}

/// A/B of the PnS route against the incumbent route at `n ∈ {3, 4}`.
///
/// Reads `POPKERN_PNS`, which is process-wide, so it must run in its
/// own process and is ignored by default. This is the measurement
/// that decides whether lowering the dispatch floor was justified:
///
/// ```text
/// POPKERN_PNS=1 cargo test -p popkern-fillet --test pns_low_valence_gates \
///     pns_low_valence_route_comparison -- --ignored --nocapture
/// POPKERN_PNS=0 cargo test -p popkern-fillet --test pns_low_valence_gates \
///     pns_low_valence_route_comparison -- --ignored --nocapture
/// ```
#[test]
#[ignore = "PK-192 route A/B: reads process-wide POPKERN_PNS; see doc comment for the two run commands"]
fn pns_low_valence_route_comparison() {
    let route = if std::env::var("POPKERN_PNS")
        .map(|v| v != "0")
        .unwrap_or(true)
    {
        "pns"
    } else {
        "incumbent"
    };
    for (n, radii) in [
        (3_usize, vec![0.05, 0.05, 0.05]),
        (3, vec![0.04, 0.06, 0.05]),
        (4, vec![0.05, 0.05, 0.05, 0.05]),
        (4, vec![0.04, 0.06, 0.05, 0.055]),
    ] {
        let o = run_pyramid_apex_outcome(n, &radii);
        let m = &o.metrics;
        println!(
            "[low-valence-ab] switch={route} n={n} radii={radii:?} served_by_pns={} methods={:?} achieved={:?} cap_faces={} valid={} seams={} g0={:.3e} g1={:.3e} g2={:.3e} cos_min={:.6}",
            o.served_by_pns(),
            o.methods,
            o.achieved,
            o.cap_faces,
            o.valid,
            m.seam_count,
            m.worst_g0_gap,
            m.worst_g1_normal,
            m.worst_g2_residual,
            m.cos_min
        );
        for d in &o.pns_diagnostics {
            println!("[low-valence-ab]     diag: {d}");
        }
        for line in &o.seam_detail {
            println!("[low-valence-ab]     seam: {line}");
        }
    }
}

/// Regression anchor for PK-197: the incumbent low-valence route is
/// G0-only across its own cap seams while reporting G1.
///
/// This gate does not assert the defect is acceptable — it pins the
/// measurement so the repair is visible when it lands, and so the
/// claim cannot quietly change without a test moving. It asserts
/// three things that are true today and must stay coupled: the cap
/// commits measurable internal seams, those seams are positionally
/// exact, and the normal misalignment is large. If a fix lands, this
/// test fails and should be replaced by the real continuity gate.
#[test]
fn incumbent_low_valence_cap_is_g0_only_across_its_own_seams() {
    for (n, radii) in [
        (3_usize, vec![0.04, 0.06, 0.05]),
        (4, vec![0.04, 0.06, 0.05, 0.055]),
    ] {
        let o = run_pyramid_apex_outcome(n, &radii);
        assert!(
            o.metrics.seam_count > 0,
            "n={n} mixed must commit measurable internal cap seams"
        );
        assert!(
            o.metrics.worst_g0_gap < G0_GATE,
            "n={n} cap seams must stay positionally exact (got {:.3e})",
            o.metrics.worst_g0_gap
        );
        assert!(
            o.metrics.worst_g1_normal > 0.5,
            "n={n} PK-197 anchor: incumbent cap was G0-only (worst normal {:.3e}); if this is now small the defect is FIXED and this anchor should be replaced by a real G1 gate",
            o.metrics.worst_g1_normal
        );
        assert!(
            o.achieved.contains(&ContinuityLevel::G1),
            "n={n} PK-197 anchor: route still claims G1 while measuring G0-only"
        );
    }
}
