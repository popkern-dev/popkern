//! PK-219 welded shared-corner quality pins: commercial cap gates
//! for the uniform rows and the suite's first mid-seam assertions.

use crate::direct_apex_validate_helpers::{self, build_under_direct_apex, native_cap_faces};
use crate::helpers::measure_cap_outer_boundary;

#[test]
fn direct_apex_pyramid3_uniform_cap_outer_commercial_gate() {
    // PK-219: the welded N=3 Coons cap measures g0 1.23e-6 /
    // g1 2.34e-7 — pin at the commercial cap gates. (The dispatch
    // G1 gate once fired for this fixture, but it was measuring the
    // restore-leak body; the welded cap itself is exact.)
    let result = build_under_direct_apex(true, &[0.10, 0.10, 0.10]);
    let cap_faces = native_cap_faces(&result).into_iter().collect::<Vec<_>>();
    let cap_outer = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
    eprintln!(
        "[direct-apex validation cap_outer] pyramid3_uniform (welded native): g0={:.3e} g1={:.3e} edges={}",
        cap_outer.worst_g0_gap, cap_outer.worst_g1_coplanarity, cap_outer.edges_measured,
    );
    assert!(
        cap_outer.edges_measured >= 3,
        "pyramid3_uniform must expose all 3 cap-to-cylinder seams"
    );
    assert!(
        cap_outer.worst_g0_gap <= 5.0e-3,
        "pyramid3_uniform cap_outer_g0 {:.3e} exceeds the commercial gate 5e-3",
        cap_outer.worst_g0_gap,
    );
    assert!(
        cap_outer.worst_g1_coplanarity <= 5.0e-2,
        "pyramid3_uniform cap_outer_g1 {:.3e} exceeds the commercial gate 5e-2",
        cap_outer.worst_g1_coplanarity,
    );
}

/// Mid-seam G1 pins per class (PK-219, tightened by PK-220). Two
/// numbers are recorded per committed multi-sector cap:
///
/// - `worst_mid_seam_normal` — measured AT COMMIT, mid-pipeline,
///   through pcurves that are not yet final. On welded polar wedges
///   this over-reads the pole half of each spoke (profile clean
///   below fraction 0.5, tens of degrees above it) — the historic
///   2.8e-1 "pole vertex-enclosure class" was exactly this artifact.
///   Pinned at the recorded ceiling as a construction-regression
///   canary only.
/// - `served_mid_seam_normal` — re-measured at the pipeline tail on
///   the final body (PK-220). This is the truth a user sees, and
///   every class meets the crate's 5e-2 G1 bar: measured pyramid3
///   3.8e-2 (pole-end sample; mid-seam 0.0), pyramid4 2.1e-3,
///   mixed 1.5e-8.
///
/// A missing measurement on a multi-face commit is an assertion
/// failure, never a pass — for BOTH numbers.
#[test]
fn direct_apex_pyramid_mid_seam_normals_pinned() {
    for (name, radii, commit_ceiling) in [
        ("pyramid3_uniform", &[0.10, 0.10, 0.10][..], 3.5e-1),
        ("pyramid4_uniform", &[0.05, 0.05, 0.05, 0.05][..], 3.5e-1),
        ("pyramid4_mixed", &[0.04, 0.06, 0.05, 0.07][..], 5.0e-2),
    ] {
        const SERVED_GATE: f64 = 5.0e-2;
        let result = build_under_direct_apex(true, radii);
        let mut measured = false;
        for route in &result.repair_report.apex_routes {
            let committed = route.exact_boundary_report.committed_faces.len();
            for diagnostic in &route.exact_boundary_report.diagnostics {
                if diagnostic.code != "exact_boundary_cap.mid_seam_normal" {
                    continue;
                }
                let detail = diagnostic.detail_code.as_deref().unwrap_or_default();
                let field = |key: &str| -> f64 {
                    detail
                        .split(';')
                        .find_map(|part| part.strip_prefix(key)?.strip_prefix('=')?.parse().ok())
                        .unwrap_or(f64::NAN)
                };
                let commit_worst = field("worst_mid_seam_normal");
                let served = field("served_mid_seam_normal");
                eprintln!("[direct-apex validation mid_seam] {name}: {detail}");
                assert!(
                    commit_worst.is_finite() && commit_worst <= commit_ceiling,
                    "{name} commit-time worst_mid_seam_normal {commit_worst:.3e} exceeds \
                     its recorded ceiling {commit_ceiling:.1e} ({detail})"
                );
                assert!(
                    served.is_finite(),
                    "{name} recorded no served_mid_seam_normal — the serving-boundary \
                     refresh must measure every surviving committed cap ({detail})"
                );
                assert!(
                    served <= SERVED_GATE,
                    "{name} served_mid_seam_normal {served:.3e} exceeds the crate G1 \
                     bar {SERVED_GATE:.1e} ({detail})"
                );
                measured = true;
            }
            if committed >= 2 {
                assert!(
                    measured,
                    "{name} committed {committed} cap faces but recorded no mid-seam \
                     measurement — a missing measurement is not a pass"
                );
            }
        }
    }
}

/// PK-219 acceptance 339: nearly-uniform radii (corner separations
/// inside the short-bridge band, 10·tol .. 2%·r_min) serve through
/// the SAME welded contract — no preflight-refusal fallback. The
/// weld point is the closest approach of the adjacent contact
/// tangent lines (on both cylinders to second order), so the
/// re-terminated arcs stay inside the seam gates. Radii outside the
/// band stay honestly bridged (the second fixture mixes both kinds
/// in one apex).
#[test]
fn direct_apex_nearly_uniform_band_serves_welded_native() {
    for (name, radii) in [
        (
            "pyramid4_near_uniform_1e4",
            &[0.05, 0.0501, 0.0502, 0.0503][..],
        ),
        (
            "pyramid4_short_bridge_mix",
            &[0.05, 0.0505, 0.051, 0.0495][..],
        ),
    ] {
        let result = build_under_direct_apex(true, radii);
        direct_apex_validate_helpers::assert_welded_native_postconditions(name, &result);
        let cap_faces = native_cap_faces(&result).into_iter().collect::<Vec<_>>();
        let cap_outer = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
        eprintln!(
            "[direct-apex validation cap_outer] {name} (nearly-uniform welded): g0={:.3e} g1={:.3e} edges={}",
            cap_outer.worst_g0_gap, cap_outer.worst_g1_coplanarity, cap_outer.edges_measured,
        );
        assert!(
            cap_outer.edges_measured >= 4,
            "{name} must expose all 4 cap-to-cylinder seams"
        );
        assert!(
            cap_outer.worst_g0_gap <= 5.0e-3,
            "{name} cap_outer_g0 {:.3e} exceeds the commercial gate 5e-3",
            cap_outer.worst_g0_gap,
        );
        assert!(
            cap_outer.worst_g1_coplanarity <= 5.0e-2,
            "{name} cap_outer_g1 {:.3e} exceeds the commercial gate 5e-2",
            cap_outer.worst_g1_coplanarity,
        );
    }
}
