use super::*;

/// A.6 probe: measure each fixture under **rebuild only** —
/// no `POPKERN_SPHERE_FIT_BOUNDARY`, no Path A.  Confirms that
/// the existing topology rebuild (Phase 6c-C) closes
/// path_a_relevant gaps when run on the original 6-vertex
/// rolling-ball topology, not the post-substitution state.
/// This is the prerequisite check for the cap-kind dispatch
/// refactor: if rebuild-only on cube_uniform reaches zero
/// path_a_relevant, the dispatch design is sound.
#[ignore]
#[test]
fn phase_6c_a6_rebuild_only_probe() {
    fn build_with_rebuild_only(
        is_pyramid: bool,
        radii: &[f64],
    ) -> (BRep, Vec<popkern_core::topology::FaceHandle>) {
        // Crucially: NO POPKERN_SPHERE_FIT_BOUNDARY. The rebuild
        // operates on the original rolling-ball 6-vertex cap-loop
        // topology, which is the pre-state its plan recognises.
        let _g1 = EnvGuard::set("POPKERN_APEX_TOPOLOGY_REBUILD", "1");
        if is_pyramid {
            let pyr = BRep::n_gonal_pyramid(radii.len(), 1.0, 0.8, Tolerance::default()).unwrap();
            let vertex = pick_vertex_with_valence(&pyr, radii.len()).unwrap();
            let edges = incident_edges(&pyr, vertex);
            fillet_corner(&pyr, &edges, radii, true)
        } else {
            let cube = BRep::unit_cube(Tolerance::default()).unwrap();
            let vertex = cube.vertices.keys().next().unwrap();
            let edges = incident_edges(&cube, vertex);
            fillet_corner(&cube, &edges, radii, true)
        }
    }

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
            name: "pyramid3_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: true,
        },
        Case {
            name: "pyramid4_uniform",
            radii: &[0.05, 0.05, 0.05, 0.05],
            is_pyramid: true,
        },
    ];

    eprintln!("\n=== A.6 rebuild-only probe (no SPHERE_FIT_BOUNDARY, no Path A) ===");
    eprintln!(
        "{:<18}  {:>5}  {:>5}  {:>5}  {:>10}",
        "fixture", "struct", "pathA", "other", "worst gap"
    );
    eprintln!("{:-<48}", "");
    for case in cases {
        let (brep, _faces) = build_with_rebuild_only(case.is_pyramid, case.radii);
        let errors = brep.validate_all();
        let m = breakdown(&errors);
        eprintln!(
            "{:<18}  {:>5}  {:>5}  {:>5}  {:>10}",
            case.name,
            m.structural,
            m.path_a_relevant,
            m.other_geometric,
            if m.worst_path_a_gap > 0.0 {
                format!("{:.3e}", m.worst_path_a_gap)
            } else {
                "-".to_string()
            },
        );
    }
}

/// A.6 probe: measure each fixture's per-bucket counts under
/// both Path A and the Phase 6c-C apex topology rebuild.
#[ignore]
#[test]
fn phase_6c_a6_path_a_plus_rebuild_probe() {
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

    eprintln!(
        "\n=== A.6 Path-A + apex-topology-rebuild probe ===\n\
         (POPKERN_APEX_PATH_A=1, POPKERN_APEX_TOPOLOGY_REBUILD=1)"
    );
    eprintln!(
        "{:<18}  {:>5}  {:>5}  {:>5}  {:>10}",
        "fixture", "struct", "pathA", "other", "worst gap"
    );
    eprintln!("{:-<48}", "");
    for case in cases {
        let (brep, _faces) = build_with_path_a_and_rebuild(case.is_pyramid, case.radii);
        let errors = brep.validate_all();
        let m = breakdown(&errors);
        eprintln!(
            "{:<18}  {:>5}  {:>5}  {:>5}  {:>10}",
            case.name,
            m.structural,
            m.path_a_relevant,
            m.other_geometric,
            if m.worst_path_a_gap > 0.0 {
                format!("{:.3e}", m.worst_path_a_gap)
            } else {
                "-".to_string()
            },
        );
    }
}

/// Historical A.6 lock-in for variable-radius cube fixtures under the
/// pre-A2 sphere-fit substitution pipeline.
#[ignore]
#[test]
fn phase_6c_a6_variable_radius_cubes_path_a_clean() {
    for (name, radii) in &[
        ("cube_mixed", &[0.08, 0.10, 0.12][..]),
        ("cube_extreme", &[0.02, 0.10, 0.20][..]),
    ] {
        let (brep, _faces) = build_with_path_a(false, radii);
        let errors = brep.validate_all();
        let m = breakdown(&errors);
        eprintln!(
            "[A.6] {}: structural={} pathA={} other={} worst_pathA_gap={:.3e}",
            name, m.structural, m.path_a_relevant, m.other_geometric, m.worst_path_a_gap,
        );
        assert_eq!(
            m.structural,
            0,
            "{}: structural must be 0; errors:\n{}",
            name,
            errors
                .iter()
                .map(|e| format!("  - {e:?}"))
                .collect::<Vec<_>>()
                .join("\n"),
        );
        assert_eq!(
            m.path_a_relevant,
            0,
            "{}: path_a_relevant must be 0 under Path A's full pipeline; errors:\n{}",
            name,
            errors
                .iter()
                .map(|e| format!("  - {e:?}"))
                .collect::<Vec<_>>()
                .join("\n"),
        );
    }
}
