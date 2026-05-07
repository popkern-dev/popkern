//! Phase 6c-A2.5 — `validate_all` matrix natively under A2.
//!
//! Locks in the commercial-grade Sphere-cap quality of A2 (per
//! ADR §6c-A2 and §6c-X.A.6's metric breakdown).  Runs the same
//! six-fixture matrix that `phase_6c_a6_validation` uses, but with
//! `FilletOptions::use_a2_apex_construction = true`.  Two
//! dispatch outcomes coexist in one matrix:
//!
//! * **Uniform fixtures** (`cube_uniform`, `pyramid3_uniform`,
//!   `pyramid4_uniform`) — A2 succeeds end-to-end, no
//!   `BlendWarning::A2FellBackToLegacy`.  The planner keeps
//!   Sphere only when collapsed Sphere targets are constructively
//!   anchored; otherwise it emits Gregory to preserve topology.
//! * **Gregory-cap fixtures** (`cube_mixed`, `cube_extreme`,
//!   `pyramid4_mixed`) — A2 should complete natively through the
//!   constructive repeated-apex normalization path.  Legacy fallback
//!   remains a rollback safety net, not an expected outcome.
//!
//! Either way: every fixture must pass the A.6 metric gate
//! (`structural == 0`, `path_a_relevant == 0`).  This is the
//! commercial-grade flip-gate's strictest form: under the new A2
//! architecture, the existing matrix doesn't regress.

#[path = "support/a2_validate_helpers.rs"]
mod a2_validate_helpers;
#[path = "kp_gate_helpers.rs"]
mod helpers;

use a2_validate_helpers::{
    a2_fallback_reasons, assert_native_cap_has_cylinder_seams, assert_native_gregory_case,
    assert_native_gregory_cube_case, breakdown, build_under_a2, count_a2_fallbacks,
};
use helpers::measure_cap_outer_boundary;
use popkern_core::geometry::CurveGeometry;

/// Whole-matrix gate: every fixture passes `validate_all` cleanly
/// under native A2 dispatch.  Strict gate: no legacy fallback,
/// `structural == 0`, and `path_a_relevant == 0` for all 6.
#[test]
fn a2_5_matrix_passes_under_a2() {
    struct Case {
        name: &'static str,
        radii: &'static [f64],
        is_pyramid: bool,
    }
    let cases: &[Case] = &[
        Case {
            name: "cube_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: false,
        },
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
            name: "pyramid4_uniform",
            radii: &[0.05, 0.05, 0.05, 0.05],
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
        let fallbacks = count_a2_fallbacks(&result);
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
            fallbacks, 0,
            "[{}] expected native A2 dispatch across the matrix; got {} A2FellBackToLegacy warning(s)",
            case.name, fallbacks,
        );
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
        assert_eq!(
            m.path_a_relevant,
            0,
            "[{}] path_a_relevant must be 0 under A2; errors:\n{}",
            case.name,
            errors
                .iter()
                .map(|e| format!("  - {e:?}"))
                .collect::<Vec<_>>()
                .join("\n"),
        );
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

/// Collapsed Sphere fixtures that should remain native Sphere under A2
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
fn a2_native_pyramid4_mixed_partition_no_fallback() {
    assert_native_gregory_case("pyramid4_mixed", true, &[0.04, 0.06, 0.05, 0.07], 4);
}

#[test]
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
}

#[test]
fn a2_native_pyramid3_uniform_anchor_incompatible_uses_gregory_no_fallback() {
    assert_native_gregory_case("pyramid3_uniform", true, &[0.10, 0.10, 0.10], 3);
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

#[test]
fn a2_native_cap_support_seam_ownership_inventory() {
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
            .collect::<std::collections::BTreeSet<_>>();
        let cap_methods = result
            .corner_patches
            .iter()
            .map(|p| (p.face, p.method.as_str()))
            .collect::<std::collections::BTreeMap<_, _>>();
        let mut findings = Vec::new();
        for (edge, edge_data) in result.brep.edges.iter() {
            let he = edge_data.half_edge;
            let Some(ha) = result.brep.half_edges.get(he) else {
                continue;
            };
            let Some(hb) = result.brep.half_edges.get(ha.twin) else {
                continue;
            };
            let Some(loop_a) = result.brep.loops.get(ha.loop_) else {
                continue;
            };
            let Some(loop_b) = result.brep.loops.get(hb.loop_) else {
                continue;
            };
            let face_a = loop_a.face;
            let face_b = loop_b.face;
            if face_a == face_b {
                continue;
            }
            let a_is_cap = cap_faces.contains(&face_a);
            let b_is_cap = cap_faces.contains(&face_b);
            if a_is_cap == b_is_cap {
                continue;
            }
            let cap_he = if a_is_cap { he } else { ha.twin };
            let support_he = if a_is_cap { ha.twin } else { he };
            let support_face = if a_is_cap { face_b } else { face_a };
            let support_kind = result
                .brep
                .faces
                .get(support_face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| match &s.geometry {
                    popkern_core::geometry::surface::SurfaceGeometry::Plane { .. } => "Plane",
                    popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. } => "Cylinder",
                    popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. } => "Sphere",
                    popkern_core::geometry::surface::SurfaceGeometry::Torus { .. } => "Torus",
                    popkern_core::geometry::surface::SurfaceGeometry::Cone { .. } => "Cone",
                    popkern_core::geometry::surface::SurfaceGeometry::Gregory(_) => "Gregory",
                    popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_) => "Nurbs",
                })
                .unwrap_or("None");
            let issue = match &edge_data.curve {
                None => Some("missing 3D curve"),
                Some(curve)
                    if matches!(curve.geometry, CurveGeometry::Nurbs(_))
                        && curve
                            .geometry
                            .as_line_segment_in_range((curve.t_start, curve.t_end), 1.0e-8)
                            .is_some() =>
                {
                    Some("straight NURBS")
                }
                Some(_) => None,
            };
            let missing_cap_pcurve = result.brep.half_edges[cap_he].pcurve.is_none();
            let missing_support_pcurve = result.brep.half_edges[support_he].pcurve.is_none();
            if issue.is_some() || missing_cap_pcurve || missing_support_pcurve {
                let start = result.brep.half_edges[cap_he].origin;
                let end = result.brep.half_edges[result.brep.half_edges[cap_he].twin].origin;
                let cap_face = if a_is_cap { face_a } else { face_b };
                let on_cap_outer = result
                    .brep
                    .faces
                    .get(cap_face)
                    .and_then(|f| result.brep.loops.get(f.outer_loop))
                    .is_some_and(|loop_data| {
                        let mut cursor = loop_data.half_edge;
                        for _ in 0..(result.brep.half_edges.len() + 1) {
                            if cursor == cap_he {
                                return true;
                            }
                            let Some(he_data) = result.brep.half_edges.get(cursor) else {
                                return false;
                            };
                            cursor = he_data.next;
                            if cursor == loop_data.half_edge {
                                return false;
                            }
                        }
                        false
                    });
                findings.push(format!(
                    "{edge:?} he_edge={:?} cap={cap_face:?}/{} support={support_face:?}/{support_kind} cap_he={cap_he:?} outer={on_cap_outer} support_he={support_he:?} verts={start:?}->{end:?} issue={} cap_pcurve={} support_pcurve={}",
                    result.brep.half_edges[he].edge,
                    cap_methods
                        .get(&cap_face)
                        .copied()
                        .unwrap_or("unknown"),
                    issue.unwrap_or("none"),
                    !missing_cap_pcurve,
                    !missing_support_pcurve
                ));
            }
        }
        eprintln!(
            "[A2 seam ownership] {}: {} finding(s){}",
            case.name,
            findings.len(),
            if findings.is_empty() {
                String::new()
            } else {
                format!(": {}", findings.join(", "))
            }
        );
        if !findings.is_empty() && std::env::var_os("POPKERN_A2_TRACE").is_some() {
            for cap_face in &cap_faces {
                let Some(face_data) = result.brep.faces.get(*cap_face) else {
                    continue;
                };
                let Some(loop_data) = result.brep.loops.get(face_data.outer_loop) else {
                    continue;
                };
                let mut cursor = loop_data.half_edge;
                let mut entries = Vec::new();
                for _ in 0..(result.brep.half_edges.len() + 1) {
                    let Some(he) = result.brep.half_edges.get(cursor) else {
                        break;
                    };
                    entries.push(format!("{:?}:{:?}", cursor, he.edge));
                    cursor = he.next;
                    if cursor == loop_data.half_edge {
                        break;
                    }
                }
                eprintln!(
                    "[A2 seam ownership loop] {} {cap_face:?}: {}",
                    case.name,
                    entries.join(" -> ")
                );
            }
        }
        assert!(
            findings.is_empty(),
            "[{}] A2 cap/support seam ownership findings: {}",
            case.name,
            findings.join(", ")
        );
    }
}

#[ignore]
#[test]
fn a2_cube_mixed_cap_outer_per_edge_breakdown() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

    let result = build_under_a2(false, &[0.08, 0.10, 0.12]);
    let cap_set = result
        .corner_patches
        .iter()
        .map(|p| p.face)
        .collect::<std::collections::BTreeSet<_>>();
    eprintln!(
        "\n=== A2 cube_mixed cap-outer per-edge breakdown: {} cap faces ===",
        cap_set.len()
    );
    let kind = |g: &SurfaceGeometry| -> &'static str {
        match g {
            SurfaceGeometry::Plane { .. } => "Plane",
            SurfaceGeometry::Cylinder { .. } => "Cyl",
            SurfaceGeometry::Sphere { .. } => "Sphere",
            SurfaceGeometry::Torus { .. } => "Torus",
            SurfaceGeometry::Cone { .. } => "Cone",
            SurfaceGeometry::Gregory(_) => "Gregory",
            SurfaceGeometry::Nurbs(_) => "Nurbs",
        }
    };
    for (eh, edge) in result.brep.edges.iter() {
        let Some(coe) = edge.curve.as_ref() else {
            continue;
        };
        let he_a = edge.half_edge;
        let Some(ha) = result.brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = result.brep.half_edges.get(ha.twin) else {
            continue;
        };
        let face_a = result.brep.loops[ha.loop_].face;
        let face_b = result.brep.loops[hb.loop_].face;
        let a_is_cap = cap_set.contains(&face_a);
        let b_is_cap = cap_set.contains(&face_b);
        if face_a == face_b || a_is_cap == b_is_cap {
            continue;
        }
        let (Some(pa), Some(pb)) = (ha.pcurve.as_ref(), hb.pcurve.as_ref()) else {
            continue;
        };
        let Some(sa) = result
            .brep
            .faces
            .get(face_a)
            .and_then(|f| f.surface.as_ref())
        else {
            continue;
        };
        let Some(sb) = result
            .brep
            .faces
            .get(face_b)
            .and_then(|f| f.surface.as_ref())
        else {
            continue;
        };
        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
        let t_a = 0.5 * (pa_lo + pa_hi);
        let t_b = 0.5 * (pb_lo + pb_hi);
        let uv_a = pa.point_at(t_a);
        let uv_b = pb.point_at(t_b);
        let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
        let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
        let (t_lo, t_hi) = (coe.t_start, coe.t_end);
        let p_mid = coe.geometry.point_at(0.5 * (t_lo + t_hi));
        eprintln!(
            "edge={eh:?} faces=({face_a:?}:{}, {face_b:?}:{}) cap_is_a={} g0={:.3e} dev_a={:.3e} dev_b={:.3e} curve_mid=({:.3},{:.3},{:.3}) a=({:.3},{:.3},{:.3}) b=({:.3},{:.3},{:.3})",
            kind(&sa.geometry),
            kind(&sb.geometry),
            a_is_cap,
            p_a.distance_to(&p_b),
            p_mid.distance_to(&p_a),
            p_mid.distance_to(&p_b),
            p_mid.x,
            p_mid.y,
            p_mid.z,
            p_a.x,
            p_a.y,
            p_a.z,
            p_b.x,
            p_b.y,
            p_b.z,
        );
    }
}

#[ignore]
#[test]
fn a2_native_cube_mixed_has_cap_to_cylinder_seams() {
    let result = build_under_a2(false, &[0.08, 0.10, 0.12]);
    let reasons = a2_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "diagnostic expects native A2 to be accepted; got fallback(s): {reasons:#?}"
    );

    assert_native_cap_has_cylinder_seams("cube_mixed", &result, 3);
}

/// Combined matrix dump (ignored, diagnostic).  Run with
/// `--ignored --nocapture` to inspect every fixture under A2.
#[ignore]
#[test]
fn a2_5_matrix_dump() {
    struct Case {
        name: &'static str,
        radii: &'static [f64],
        is_pyramid: bool,
    }
    let cases: &[Case] = &[
        Case {
            name: "cube_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: false,
        },
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
            name: "pyramid4_uniform",
            radii: &[0.05, 0.05, 0.05, 0.05],
            is_pyramid: true,
        },
        Case {
            name: "pyramid4_mixed",
            radii: &[0.04, 0.06, 0.05, 0.07],
            is_pyramid: true,
        },
    ];
    eprintln!("\n=== A2.5 validate_all matrix natively under A2 ===");
    eprintln!(
        "{:<18}  {:>8}  {:>5}  {:>5}  {:>5}  {:>10}",
        "fixture", "fallback", "struct", "pathA", "other", "worst gap"
    );
    eprintln!("{:-<60}", "");
    for case in cases {
        let result = build_under_a2(case.is_pyramid, case.radii);
        let fallbacks = count_a2_fallbacks(&result);
        let errors = result.brep.validate_all();
        let m = breakdown(&errors);
        eprintln!(
            "{:<18}  {:>8}  {:>5}  {:>5}  {:>5}  {:>10}",
            case.name,
            fallbacks,
            m.structural,
            m.path_a_relevant,
            m.other_geometric,
            if m.worst_path_a_gap > 0.0 {
                format!("{:.3e}", m.worst_path_a_gap)
            } else {
                "—".to_string()
            },
        );
        for error in errors {
            eprintln!("  - {error:?}");
        }
    }
}
