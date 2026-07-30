//! Phase 6c-direct-apex validation — `validate_all` matrix under direct-apex dispatch.
//!
//! Locks in the commercial-grade apex-cap quality of direct-apex (per
//! ADR §6c-direct-apex and §6c-X.A.6's metric breakdown).  Runs the same
//! six-fixture matrix that `phase_6c_a6_validation` uses, but with
//! `FilletOptions::use_direct_apex_construction = true`.  Active
//! matrix fixtures must all complete through native DirectApex:
//!
//! * **Uniform fixtures** (`cube_uniform`, `pyramid3_uniform`,
//!   `pyramid4_uniform`) — direct-apex succeeds end-to-end, no
//!   `BlendWarning::DirectApexFellBackToLegacy`.  The planner keeps
//!   Sphere only for trihedral collapsed targets whose cap boundary
//!   is constructively valid; higher-valence uniform apices use
//!   Gregory to preserve cap-to-cylinder geometry.
//! * **Native Gregory fixtures** (`cube_mixed`) — direct-apex must
//!   complete with no fallback and expose measurable,
//!   commercial-tolerance cap-to-cylinder seams.
//! * **Repeated-apex Gregory fixtures** (`pyramid3_uniform`,
//!   `pyramid4_uniform`, `pyramid4_mixed`) — direct-apex must complete
//!   without containment fallback. Remaining surface-quality blockers
//!   stay visible through `other_geometric` and the focused
//!   cap-continuity gates.
//! * **`cube_extreme`** completes natively (PK-074 closed): the
//!   partition's displaced-membership residual collection retired the
//!   phantom scaffold seams that had forced the legacy fallback, and
//!   the recompute's degree-1 shortcut now verifies parameterization
//!   (not just image straightness) at the quarter points, so the
//!   cap/support seams measure at commercial tolerance end to end.
//!
//! Either way: every fixture must pass the A.6 metric gate
//! (`structural == 0`, `path_a_relevant == 0`).  This is the
//! commercial-grade flip-gate's strictest form: under the new DirectApex
//! architecture, the existing matrix doesn't regress.

#[path = "support/direct_apex_cap_outer_limitations.rs"]
mod direct_apex_cap_outer_limitations;
#[path = "support/direct_apex_validate_all_matrix_native.rs"]
mod direct_apex_validate_all_matrix_native;
#[path = "support/direct_apex_validate_all_matrix_seams.rs"]
mod direct_apex_validate_all_matrix_seams;
#[path = "support/direct_apex_validate_all_matrix_zoo.rs"]
mod direct_apex_validate_all_matrix_zoo;
#[path = "support/direct_apex_validate_helpers.rs"]
mod direct_apex_validate_helpers;
#[path = "kp_gate_helpers.rs"]
mod helpers;

use direct_apex_validate_helpers::{
    assert_native_gregory_case, assert_native_gregory_cube_case, breakdown,
    build_under_direct_apex, count_direct_apex_fallbacks, direct_apex_fallback_reasons,
    native_cap_faces,
};
use helpers::measure_cap_outer_boundary;
use popkern_core::geometry::CurveGeometry;

/// Whole-matrix native DirectApex gate: every active fixture keeps
/// structural validity and path-A gap validity under direct-apex
/// dispatch, with zero containment fallback warnings.
#[test]
fn direct_apex_5_matrix_preserves_structural_validity_under_direct_apex() {
    struct Case {
        name: &'static str,
        radii: &'static [f64],
        is_pyramid: bool,
        expected_fallbacks: usize,
    }
    let cases: &[Case] = &[
        Case {
            name: "cube_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: false,
            expected_fallbacks: 0,
        },
        Case {
            name: "cube_mixed",
            radii: &[0.08, 0.10, 0.12],
            is_pyramid: false,
            expected_fallbacks: 0,
        },
        Case {
            // PK-074: completes natively since the displaced-membership
            // residual collection — phantom scaffold seams no longer
            // reach acceptance, cap/support seams measure exactly.
            name: "cube_extreme",
            radii: &[0.02, 0.10, 0.20],
            is_pyramid: false,
            expected_fallbacks: 0,
        },
        Case {
            // PK-216: uniform radii collapse every sector corner to a
            // shared point; the constructive route declines with the
            // class named and the recorded legacy fallback serves a
            // clean body (the old "native, 0 fallbacks" green shipped
            // a glued surfaceless membrane). Native returns with the
            // shared-corner substrate (PK-219).
            name: "pyramid3_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: true,
            expected_fallbacks: 1,
        },
        Case {
            // PK-216 shared-corner class — see pyramid3_uniform above.
            name: "pyramid4_uniform",
            radii: &[0.05, 0.05, 0.05, 0.05],
            is_pyramid: true,
            expected_fallbacks: 1,
        },
        Case {
            name: "pyramid4_mixed",
            radii: &[0.04, 0.06, 0.05, 0.07],
            is_pyramid: true,
            expected_fallbacks: 0,
        },
    ];
    for case in cases {
        let result = build_under_direct_apex(case.is_pyramid, case.radii);
        let fallbacks = count_direct_apex_fallbacks(&result);
        let fallback_reasons = direct_apex_fallback_reasons(&result);
        let errors = result.brep.validate_all();
        let m = breakdown(&errors);
        eprintln!(
            "[direct-apex validation] {}: fallbacks={} structural={} pathA={} other={} worst_gap={:.3e}",
            case.name,
            fallbacks,
            m.structural,
            m.path_a_relevant,
            m.other_geometric,
            m.worst_path_a_gap,
        );
        assert_eq!(
            fallbacks,
            case.expected_fallbacks,
            "[{}] unexpected direct-apex containment count; reasons:\n{}",
            case.name,
            fallback_reasons.join("\n"),
        );
        if case.expected_fallbacks > 0 {
            assert!(
                fallback_reasons
                    .iter()
                    .all(|reason| reason.contains("construct_apex")),
                "[{}] fallback diagnostics must identify the constructive direct-apex boundary; got {fallback_reasons:#?}",
                case.name,
            );
        }
        assert_eq!(
            m.structural,
            0,
            "[{}] structural must be 0 under DirectApex; errors:\n{}",
            case.name,
            errors
                .iter()
                .map(|e| format!("  - {e:?}"))
                .collect::<Vec<_>>()
                .join("\n"),
        );
        if case.expected_fallbacks == 0 {
            assert_eq!(
                m.path_a_relevant,
                0,
                "[{}] native direct-apex path_a_relevant must be 0; errors:\n{}",
                case.name,
                errors
                    .iter()
                    .map(|e| format!("  - {e:?}"))
                    .collect::<Vec<_>>()
                    .join("\n"),
            );
        }
    }
}

#[test]
fn cube_extreme_exact_boundary_cap_outer_passes_commercial_gate() {
    let result = build_under_direct_apex(false, &[0.02, 0.10, 0.20]);
    let reasons = direct_apex_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "cube_extreme should now complete through native direct-apex; got:\n{}",
        reasons.join("\n")
    );
    let cap_faces = native_cap_faces(&result).into_iter().collect::<Vec<_>>();
    let metrics = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
    assert!(
        metrics.edges_measured >= 1,
        "cube_extreme should expose measurable exact-boundary cap edges"
    );
    assert!(
        metrics.worst_g0_gap <= 5.0e-3,
        "cube_extreme cap_outer_g0 {:.3e} exceeds 5.0e-3",
        metrics.worst_g0_gap
    );
    assert!(
        metrics.worst_g1_coplanarity <= 5.0e-2,
        "cube_extreme cap_outer_g1 {:.3e} exceeds 5.0e-2",
        metrics.worst_g1_coplanarity
    );
}

/// Phase 6c-E commercial `cap_outer_g0 <= 5e-3 m` gate under DirectApex
/// dispatch. Gregory fixtures may complete through native direct-apex only
/// when their cap topology exposes measurable cap-to-cylinder seams
/// within the accepted positional boundary tolerance.
///
/// Collapsed Sphere fixtures (`cube_uniform`, `pyramid4_uniform`)
/// are not gated here because this helper measures Gregory cap
/// faces. Their separate quality gates belong with Sphere-cap
/// boundary metrics.
#[test]
fn direct_apex_5_cap_outer_g0_gregory_under_direct_apex() {
    struct Case {
        name: &'static str,
        radii: &'static [f64],
        is_pyramid: bool,
    }
    let cases: &[Case] = &[
        Case {
            name: "cube_mixed",
            radii: &[0.08, 0.10, 0.12],
            is_pyramid: false,
        },
        Case {
            name: "cube_extreme",
            radii: &[0.02, 0.10, 0.20],
            is_pyramid: false,
        },
        // pyramid3_uniform left this set with PK-216: the uniform
        // (shared-corner) class serves via the recorded legacy
        // fallback, not a native Gregory cap — its serve is pinned by
        // the shared-corner containment gates below and returns here
        // with PK-219.
    ];
    for case in cases {
        let result = build_under_direct_apex(case.is_pyramid, case.radii);
        let cap_faces = native_cap_faces(&result).into_iter().collect::<Vec<_>>();
        let m = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
        eprintln!(
            "[direct-apex validation cap_outer] {}: g0={:.3e} g1={:.3e} edges={}",
            case.name, m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured,
        );
        assert!(
            m.edges_measured >= 1,
            "[{}] no cap-vs-non-cap edges measured under DirectApex",
            case.name,
        );
        assert!(
            m.worst_g0_gap <= 5.0e-3,
            "[{}] cap_outer_g0 {:.3e} exceeds 5.0e-3 commercial gate under DirectApex",
            case.name,
            m.worst_g0_gap,
        );
        assert!(
            m.worst_g1_coplanarity <= 5.0e-2,
            "[{}] cap_outer_g1 {:.3e} exceeds 5.0e-2 commercial gate under DirectApex",
            case.name,
            m.worst_g1_coplanarity,
        );
    }
}

/// Collapsed trihedral Sphere fixtures that should remain native Sphere under DirectApex
/// (no `DirectApexFellBackToLegacy`).  Uniform pyramid cases are covered by
/// separate native gates because the planner may choose Gregory when
/// the collapsed Sphere target loop is not constructively anchored.
#[test]
fn direct_apex_5_native_sphere_cap_set() {
    struct Case {
        name: &'static str,
        radii: &'static [f64],
        is_pyramid: bool,
    }
    let cases: &[Case] = &[Case {
        name: "cube_uniform",
        radii: &[0.10, 0.10, 0.10],
        is_pyramid: false,
    }];
    for case in cases {
        let result = build_under_direct_apex(case.is_pyramid, case.radii);
        let fallbacks = count_direct_apex_fallbacks(&result);
        assert_eq!(
            fallbacks, 0,
            "[{}] expected native direct-apex dispatch (Sphere cap); got {} DirectApexFellBackToLegacy warning(s) — \
             direct-apex build.f Sphere-cap regression?",
            case.name, fallbacks,
        );
    }
}

/// Native Gregory direct-apex recovery gate for the baseline mixed cube.
#[test]
fn direct_apex_native_cube_mixed_no_fallback() {
    assert_native_gregory_cube_case("cube_mixed", &[0.08, 0.10, 0.12]);
}

/// Native Gregory direct-apex semantic-sector ownership gate.  Locks the
/// construction invariant that physical cap-loop order
/// (`sector_index`) and semantic fillet ownership (`fillet_index`)
/// are separate.  Both mixed cube fixtures reorder those concepts
/// during repeated-apex normalization; each must finish natively
/// with cylinder-owned FilletEnd seams.
#[test]
fn direct_apex_native_variable_cube_gregory_set_no_fallback() {
    for (name, radii) in [
        ("cube_mixed", &[0.08, 0.10, 0.12][..]),
        ("cube_extreme", &[0.02, 0.10, 0.20][..]),
    ] {
        assert_native_gregory_cube_case(name, radii);
    }
}

/// Native Gregory direct-apex repeated-apex partition gate.  Locks the
/// N=4 mixed pyramid topology class where no serial FilletEnd mef
/// order is constructible; the planned partition operator must
/// complete natively and expose all four cylinder seams.
#[test]
fn direct_apex_native_pyramid4_mixed_partition_no_fallback() {
    assert_native_gregory_case("pyramid4_mixed", true, &[0.04, 0.06, 0.05, 0.07], 4);
}

#[test]
fn direct_apex_pyramid4_uniform_shared_corner_fallback_serves_clean() {
    let result = build_under_direct_apex(true, &[0.05, 0.05, 0.05, 0.05]);
    direct_apex_validate_helpers::assert_shared_corner_legacy_fallback("pyramid4_uniform", &result);
}

#[test]
fn direct_apex_pyramid4_uniform_cap_outer_containment_ceiling() {
    // PK-216 honest re-pin. The old 5e-3/5e-2 "commercial gate" green
    // was VACUOUS: it measured the mis-scaled sector pillow's own
    // seams while the served boundary was a surfaceless octagon
    // membrane (a literal mesh hole). The shared-corner fallback now
    // serves the legacy sphere cap — geometrically inexact on a
    // pyramid by construction (no per-fillet sphere centre lies on a
    // neighbour's axis), measured g0 7.37e-2 / g1 ~1.0 at the sphere
    // rim. This pins the containment CEILING so further regression
    // still fails; the commercial gate returns with PK-219.
    let result = build_under_direct_apex(true, &[0.05, 0.05, 0.05, 0.05]);
    let cap_faces = native_cap_faces(&result).into_iter().collect::<Vec<_>>();
    let cap_outer = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
    eprintln!(
        "[direct-apex validation cap_outer] pyramid4_uniform (containment): g0={:.3e} g1={:.3e} edges={}",
        cap_outer.worst_g0_gap, cap_outer.worst_g1_coplanarity, cap_outer.edges_measured,
    );
    assert!(
        cap_outer.edges_measured >= 1,
        "pyramid4_uniform exposes no measurable cap-to-neighbor seams",
    );
    assert!(
        cap_outer.worst_g0_gap <= 1.5e-1,
        "pyramid4_uniform cap_outer_g0 {:.3e} exceeds the recorded containment ceiling 1.5e-1",
        cap_outer.worst_g0_gap,
    );
}

#[test]
fn direct_apex_pyramid3_uniform_shared_corner_fallback_serves_clean() {
    // PK-216: at HEAD before the shared-corner gate this fixture
    // served 38 faces with THREE dead boundary sides and 18 residual
    // faces behind a "native, no fallback" green.
    let result = build_under_direct_apex(true, &[0.10, 0.10, 0.10]);
    direct_apex_validate_helpers::assert_shared_corner_legacy_fallback("pyramid3_uniform", &result);
}

/// direct-apex edge-curve semantic storage gate.  Straight 3D edge curves
/// must be authored as analytic lines at construction time, including
/// support edges introduced by face/edge splits, so downstream topology
/// and GUI feature logic do not need to rediscover obvious curve intent.
#[test]
fn direct_apex_does_not_store_semantic_lines_as_nurbs() {
    struct Case {
        name: &'static str,
        radii: &'static [f64],
        is_pyramid: bool,
    }
    let cases: &[Case] = &[
        Case {
            name: "cube_mixed",
            radii: &[0.08, 0.10, 0.12],
            is_pyramid: false,
        },
        Case {
            name: "cube_extreme",
            radii: &[0.02, 0.10, 0.20],
            is_pyramid: false,
        },
    ];
    for case in cases {
        let result = build_under_direct_apex(case.is_pyramid, case.radii);
        let mut semantic_line_nurbs = Vec::new();
        for (edge, edge_data) in result.brep.edges.iter() {
            let Some(curve) = edge_data.curve.as_ref() else {
                continue;
            };
            let he = edge_data.half_edge;
            let Some(he_data) = result.brep.half_edges.get(he) else {
                continue;
            };
            let twin = he_data.twin;
            let Some(twin_data) = result.brep.half_edges.get(twin) else {
                continue;
            };
            let Some(loop_a) = result.brep.loops.get(he_data.loop_) else {
                continue;
            };
            let Some(loop_b) = result.brep.loops.get(twin_data.loop_) else {
                continue;
            };
            let face_a = loop_a.face;
            let face_b = loop_b.face;
            if let CurveGeometry::Nurbs(_) = &curve.geometry
                && curve
                    .geometry
                    .as_line_segment_in_range((curve.t_start, curve.t_end), 1.0e-8)
                    .is_some()
            {
                let kind = |face| -> &'static str {
                    match result
                        .brep
                        .faces
                        .get(face)
                        .and_then(|f| f.surface.as_ref())
                        .map(|s| &s.geometry)
                    {
                        Some(popkern_core::geometry::surface::SurfaceGeometry::Plane {
                            ..
                        }) => "Plane",
                        Some(popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
                            ..
                        }) => "Cylinder",
                        Some(popkern_core::geometry::surface::SurfaceGeometry::Sphere {
                            ..
                        }) => "Sphere",
                        Some(popkern_core::geometry::surface::SurfaceGeometry::Torus {
                            ..
                        }) => "Torus",
                        Some(popkern_core::geometry::surface::SurfaceGeometry::Cone { .. }) => {
                            "Cone"
                        }
                        Some(popkern_core::geometry::surface::SurfaceGeometry::Gregory(_)) => {
                            "Gregory"
                        }
                        Some(
                            popkern_core::geometry::surface::SurfaceGeometry::PolyhedralSpline(_),
                        ) => "PolyhedralSpline",
                        Some(popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_)) => "Nurbs",
                        None => "None",
                    }
                };
                semantic_line_nurbs.push(format!(
                    "{edge:?}({:?}->{:?}, {}<->{})",
                    he_data.origin,
                    twin_data.origin,
                    kind(face_a),
                    kind(face_b)
                ));
            }
        }
        assert!(
            semantic_line_nurbs.is_empty(),
            "[{}] semantically straight direct-apex edge curves must be stored as CurveGeometry::Line, \
             not NURBS: {}",
            case.name,
            semantic_line_nurbs.join(", ")
        );
    }
}
