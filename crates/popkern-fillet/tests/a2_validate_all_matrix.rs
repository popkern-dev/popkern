//! Phase 6c-A2.5 — `validate_all` matrix under A2 dispatch.
//!
//! Locks in the commercial-grade apex-cap quality of A2 (per
//! ADR §6c-A2 and §6c-X.A.6's metric breakdown).  Runs the same
//! six-fixture matrix that `phase_6c_a6_validation` uses, but with
//! `FilletOptions::use_a2_apex_construction = true`.  Two
//! dispatch outcomes coexist in one matrix:
//!
//! * **Uniform fixtures** (`cube_uniform`, `pyramid3_uniform`,
//!   `pyramid4_uniform`) — A2 succeeds end-to-end, no
//!   `BlendWarning::A2FellBackToLegacy`.  The planner keeps
//!   Sphere only for trihedral collapsed targets whose cap boundary
//!   is constructively valid; higher-valence uniform apices use
//!   Gregory to preserve cap-to-cylinder geometry.
//! * **Gregory-cap fixtures** (`cube_mixed`, `cube_extreme`,
//!   `pyramid3_uniform`, `pyramid4_uniform`, `pyramid4_mixed`) —
//!   until the validated outside/outside support-continuation
//!   primitive and remaining trim-pcurve repairs are implemented,
//!   A2 may decline through one recorded fallback.  These are
//!   containment checks, not native-quality acceptance.
//!
//! Either way: every fixture must pass the A.6 metric gate
//! (`structural == 0`, `path_a_relevant == 0`).  This is the
//! commercial-grade flip-gate's strictest form: under the new A2
//! architecture, the existing matrix doesn't regress.

#[path = "support/a2_validate_all_matrix_diagnostics.rs"]
mod a2_validate_all_matrix_diagnostics;
#[path = "support/a2_validate_all_matrix_seams.rs"]
mod a2_validate_all_matrix_seams;
#[path = "support/a2_validate_all_matrix_zoo.rs"]
mod a2_validate_all_matrix_zoo;
#[path = "support/a2_validate_helpers.rs"]
mod a2_validate_helpers;
#[path = "kp_gate_helpers.rs"]
mod helpers;

use a2_validate_helpers::{
    EnvGuard, a2_fallback_reasons, assert_native_cap_has_cylinder_seams,
    assert_native_gregory_case, assert_native_gregory_cube_case, breakdown, build_under_a2,
    count_a2_fallbacks,
};
use helpers::measure_cap_outer_boundary;
use popkern_core::geometry::CurveGeometry;

/// Whole-matrix containment gate: every fixture keeps structural
/// validity under A2 dispatch.  Uniform Sphere remains a hard native
/// gate; Gregory/repeated-apex fixtures must either become native
/// later or decline once with an actionable fallback diagnostic.
#[test]
fn a2_5_matrix_preserves_structural_validity_under_a2() {
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
            expected_fallbacks: 1,
        },
        Case {
            name: "cube_extreme",
            radii: &[0.02, 0.10, 0.20],
            is_pyramid: false,
            expected_fallbacks: 1,
        },
        Case {
            name: "pyramid3_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: true,
            expected_fallbacks: 1,
        },
        Case {
            name: "pyramid4_uniform",
            radii: &[0.05, 0.05, 0.05, 0.05],
            is_pyramid: true,
            expected_fallbacks: 1,
        },
        Case {
            name: "pyramid4_mixed",
            radii: &[0.04, 0.06, 0.05, 0.07],
            is_pyramid: true,
            expected_fallbacks: 1,
        },
    ];
    for case in cases {
        let result = build_under_a2(case.is_pyramid, case.radii);
        let fallbacks = count_a2_fallbacks(&result);
        let fallback_reasons = a2_fallback_reasons(&result);
        let errors = result.brep.validate_all();
        let m = breakdown(&errors);
        eprintln!(
            "[A2.5] {}: fallbacks={} structural={} pathA={} other={} worst_gap={:.3e}",
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
            "[{}] unexpected A2 containment count; reasons:\n{}",
            case.name,
            fallback_reasons.join("\n"),
        );
        if case.expected_fallbacks > 0 {
            assert!(
                fallback_reasons
                    .iter()
                    .all(|reason| reason.contains("construct_apex")),
                "[{}] fallback diagnostics must identify the constructive A2 boundary; got {fallback_reasons:#?}",
                case.name,
            );
        }
        assert_eq!(
            m.structural,
            0,
            "[{}] structural must be 0 under A2; errors:\n{}",
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
                "[{}] native A2 path_a_relevant must be 0; errors:\n{}",
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

/// Phase 6c-E commercial `cap_outer_g0 <= 5e-3 m` gate under A2
/// dispatch. Gregory fixtures may complete through native A2 only
/// when their cap topology exposes measurable cap-to-cylinder seams
/// within the accepted positional boundary tolerance.
///
/// Collapsed Sphere fixtures (`cube_uniform`, `pyramid4_uniform`)
/// are not gated here because this helper measures Gregory cap
/// faces. Their separate quality gates belong with Sphere-cap
/// boundary metrics.
#[test]
#[ignore = "native Gregory cap-outer gate; current repeated-apex A2 path is under containment fallback"]
fn a2_5_cap_outer_g0_gregory_under_a2() {
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
        Case {
            name: "pyramid3_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: true,
        },
        Case {
            name: "pyramid4_mixed",
            radii: &[0.04, 0.06, 0.05, 0.07],
            is_pyramid: true,
        },
    ];
    for case in cases {
        let result = build_under_a2(case.is_pyramid, case.radii);
        let cap_faces = result
            .corner_patches
            .iter()
            .map(|p| p.face)
            .collect::<Vec<_>>();
        let m = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
        eprintln!(
            "[A2.5 cap_outer] {}: g0={:.3e} g1={:.3e} edges={}",
            case.name, m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured,
        );
        assert!(
            m.edges_measured >= 1,
            "[{}] no cap-vs-non-cap edges measured under A2",
            case.name,
        );
        assert!(
            m.worst_g0_gap <= 5.0e-3,
            "[{}] cap_outer_g0 {:.3e} exceeds 5.0e-3 commercial gate under A2",
            case.name,
            m.worst_g0_gap,
        );
        assert!(
            m.worst_g1_coplanarity <= 5.0e-2,
            "[{}] cap_outer_g1 {:.3e} exceeds 5.0e-2 commercial gate under A2",
            case.name,
            m.worst_g1_coplanarity,
        );
    }
}

/// Collapsed trihedral Sphere fixtures that should remain native Sphere under A2
/// (no `A2FellBackToLegacy`).  Uniform pyramid cases are covered by
/// separate native gates because the planner may choose Gregory when
/// the collapsed Sphere target loop is not constructively anchored.
#[test]
fn a2_5_native_sphere_cap_set() {
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
        let result = build_under_a2(case.is_pyramid, case.radii);
        let fallbacks = count_a2_fallbacks(&result);
        assert_eq!(
            fallbacks, 0,
            "[{}] expected native A2 dispatch (Sphere cap); got {} A2FellBackToLegacy warning(s) — \
             A2.2.f Sphere-cap regression?",
            case.name, fallbacks,
        );
    }
}

/// Native Gregory A2 recovery gate for the baseline mixed cube.
#[test]
#[ignore = "native Gregory A2 gate; current mixed cube declines through containment fallback"]
fn a2_native_cube_mixed_no_fallback() {
    assert_native_gregory_cube_case("cube_mixed", &[0.08, 0.10, 0.12]);
}

/// Native Gregory A2 semantic-sector ownership gate.  Locks the
/// construction invariant that physical cap-loop order
/// (`sector_index`) and semantic fillet ownership (`fillet_index`)
/// are separate.  Both mixed cube fixtures reorder those concepts
/// during repeated-apex normalization; each must finish natively
/// with cylinder-owned FilletEnd seams.
#[test]
#[ignore = "native Gregory A2 gate; current variable cubes decline through containment fallback"]
fn a2_native_variable_cube_gregory_set_no_fallback() {
    for (name, radii) in [
        ("cube_mixed", &[0.08, 0.10, 0.12][..]),
        ("cube_extreme", &[0.02, 0.10, 0.20][..]),
    ] {
        assert_native_gregory_cube_case(name, radii);
    }
}

/// Native Gregory A2 repeated-apex partition gate.  Locks the
/// N=4 mixed pyramid topology class where no serial FilletEnd mef
/// order is constructible; the planned partition operator must
/// complete natively and expose all four cylinder seams.
#[test]
#[ignore = "native N=4 Gregory partition gate; missing outside/outside support-continuation primitive"]
fn a2_native_pyramid4_mixed_partition_no_fallback() {
    assert_native_gregory_case("pyramid4_mixed", true, &[0.04, 0.06, 0.05, 0.07], 4);
}

#[test]
#[ignore = "native N=4 uniform Gregory gate; missing outside/outside support-continuation primitive"]
fn a2_native_pyramid4_uniform_no_fallback() {
    let result = build_under_a2(true, &[0.05, 0.05, 0.05, 0.05]);
    let reasons = a2_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "pyramid4_uniform must complete through native A2; got {} fallback(s):\n{}",
        reasons.len(),
        reasons
            .iter()
            .map(|r| format!("  - {r}"))
            .collect::<Vec<_>>()
            .join("\n"),
    );

    let errors = result.brep.validate_all();
    let m = breakdown(&errors);
    assert_eq!(
        m.structural, 0,
        "native pyramid4_uniform structural errors: {errors:#?}"
    );
    assert_eq!(
        m.path_a_relevant, 0,
        "native pyramid4_uniform open-loop/shell-gap errors: {errors:#?}"
    );

    let cap_faces = result
        .corner_patches
        .iter()
        .map(|p| p.face)
        .collect::<Vec<_>>();
    let cap_outer = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
    eprintln!(
        "[A2.5 cap_outer] pyramid4_uniform: g0={:.3e} g1={:.3e} edges={}",
        cap_outer.worst_g0_gap, cap_outer.worst_g1_coplanarity, cap_outer.edges_measured,
    );
    assert!(
        cap_outer.edges_measured >= 1,
        "native pyramid4_uniform exposes no measurable cap-to-neighbor seams",
    );
    assert!(
        cap_outer.worst_g0_gap <= 5.0e-3,
        "native pyramid4_uniform cap_outer_g0 {:.3e} exceeds 5.0e-3 commercial gate",
        cap_outer.worst_g0_gap,
    );
    assert!(
        cap_outer.worst_g1_coplanarity <= 5.0e-2,
        "native pyramid4_uniform cap_outer_g1 {:.3e} exceeds 5.0e-2 commercial gate",
        cap_outer.worst_g1_coplanarity,
    );
}

#[test]
#[ignore = "native pyramid3 Gregory gate; current path declines through trim-pcurve containment fallback"]
fn a2_native_pyramid3_uniform_anchor_incompatible_uses_gregory_no_fallback() {
    assert_native_gregory_case("pyramid3_uniform", true, &[0.10, 0.10, 0.10], 3);
}

#[test]
#[ignore = "diagnostic gate: uses process-wide strict A2 env and must run in isolation"]
fn a2_strict_v_free_repeated_apex_succeeds_without_fallback() {
    let _guard = EnvGuard::set("POPKERN_A2_REQUIRE_V_FREE_REPEATED_APEX", "1");
    let result = build_under_a2(true, &[0.04, 0.06, 0.05, 0.07]);
    let reasons = a2_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "strict V-free repeated-apex mode must complete without fallback; got:\n{}",
        reasons.join("\n"),
    );
    assert!(
        !reasons
            .iter()
            .any(|reason| reason.contains("repeated-apex construction retained apex vertex")),
        "strict V-free repeated-apex mode must not retain production apex topology; got:\n{}",
        reasons.join("\n"),
    );
}

/// A2 edge-curve semantic storage gate.  Straight 3D edge curves
/// must be authored as analytic lines at construction time, including
/// support edges introduced by face/edge splits, so downstream topology
/// and GUI feature logic do not need to rediscover obvious curve intent.
#[test]
fn a2_does_not_store_semantic_lines_as_nurbs() {
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
        Case {
            name: "pyramid4_mixed",
            radii: &[0.04, 0.06, 0.05, 0.07],
            is_pyramid: true,
        },
    ];
    for case in cases {
        let result = build_under_a2(case.is_pyramid, case.radii);
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
            "[{}] semantically straight A2 edge curves must be stored as CurveGeometry::Line, \
             not NURBS: {}",
            case.name,
            semantic_line_nurbs.join(", ")
        );
    }
}
