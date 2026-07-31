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

/// First mid-seam assertions in the suite (PK-219): the mid-seam
/// normal profile is measured between adjacent committed cap faces
/// (transverse continuity the along-seam metric cannot see). It has
/// been informational since the facade landed; these pins make it a
/// gated number per class. Bridged multi-sector caps (pyramid4_mixed)
/// meet the crate's 5e-2 G1 bar outright. WELDED multi-sector caps
/// are polar wedges whose inner spokes meet at a shared pole: the
/// welded LeftColumnToColumn cross-spoke pairing (this task) cut the
/// crease 8.4e-1 -> 2.8e-1, and the remainder is the pole
/// vertex-enclosure class — pairwise midpoint offsets around a pole
/// cannot satisfy all sectors simultaneously (the PK-192/PnS shared-
/// spoke frontier; follow-on filed from PK-219). Pinned at the
/// measured class ceiling so regressions still fail. A missing
/// measurement on a multi-face commit is an assertion failure, never
/// a pass.
#[test]
fn direct_apex_pyramid_mid_seam_normals_pinned() {
    for (name, radii, gate) in [
        ("pyramid3_uniform", &[0.10, 0.10, 0.10][..], 3.5e-1),
        ("pyramid4_uniform", &[0.05, 0.05, 0.05, 0.05][..], 3.5e-1),
        ("pyramid4_mixed", &[0.04, 0.06, 0.05, 0.07][..], 5.0e-2),
    ] {
        let result = build_under_direct_apex(true, radii);
        let mut measured = false;
        for route in &result.repair_report.apex_routes {
            let committed = route.exact_boundary_report.committed_faces.len();
            for diagnostic in &route.exact_boundary_report.diagnostics {
                if diagnostic.code != "exact_boundary_cap.mid_seam_normal" {
                    continue;
                }
                let detail = diagnostic.detail_code.as_deref().unwrap_or_default();
                let worst: f64 = detail
                    .split("worst_mid_seam_normal=")
                    .nth(1)
                    .and_then(|value| value.parse().ok())
                    .unwrap_or(f64::NAN);
                eprintln!("[direct-apex validation mid_seam] {name}: {detail}");
                assert!(
                    worst.is_finite() && worst <= gate,
                    "{name} worst_mid_seam_normal {worst:.3e} exceeds its pinned gate \
                     {gate:.1e} ({detail})"
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
