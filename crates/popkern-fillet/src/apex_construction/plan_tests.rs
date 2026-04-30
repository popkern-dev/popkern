//! Unit tests for [`plan_apex_resolution`].
//!
//! Per ADR §11 acceptance criterion #10, every degenerate case in
//! ADR §4.4 has a corresponding test here.  A2.1.b through A2.1.f
//! flesh out the test bodies as the planner gains capability.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::types::{FilletRadius, FilletSpec};

use super::plan::plan_apex_resolution;

#[test]
fn stub_returns_insufficient_fillets_for_n_lt_3() {
    // Build a minimal BRep with one vertex.
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep
        .mvfs(popkern_core::geometry::Point3::new(0.0, 0.0, 0.0))
        .unwrap();
    // Empty fillet_specs → InsufficientFillets.
    let err =
        plan_apex_resolution(&brep, mvfs.vertex, &[], &[], &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InsufficientFillets { n: 0, .. }
        ),
        "expected InsufficientFillets, got {err:?}"
    );
    // 2 fillets → still insufficient.
    let mev1 = brep
        .mev_at_vertex(
            mvfs.vertex,
            mvfs.loop_,
            popkern_core::geometry::Point3::new(1.0, 0.0, 0.0),
        )
        .unwrap();
    let dummy_specs = vec![
        FilletSpec {
            edge: mev1.edge,
            radius: FilletRadius::Constant(0.1),
        },
        FilletSpec {
            edge: mev1.edge,
            radius: FilletRadius::Constant(0.1),
        },
    ];
    let err = plan_apex_resolution(&brep, mvfs.vertex, &dummy_specs, &[], &Tolerance::default())
        .unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InsufficientFillets { n: 2, .. }
        ),
        "expected InsufficientFillets(n=2), got {err:?}"
    );
}

#[test]
fn stub_returns_invalid_input_for_n_ge_3() {
    // 3 fillets passes the InsufficientFillets check but the stub
    // bails at "planner not yet implemented".
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep
        .mvfs(popkern_core::geometry::Point3::new(0.0, 0.0, 0.0))
        .unwrap();
    let mev1 = brep
        .mev_at_vertex(
            mvfs.vertex,
            mvfs.loop_,
            popkern_core::geometry::Point3::new(1.0, 0.0, 0.0),
        )
        .unwrap();
    let dummy_edge = mev1.edge;
    let specs = vec![
        FilletSpec {
            edge: dummy_edge,
            radius: FilletRadius::Constant(0.1),
        };
        3
    ];
    // 3 fillets passes the InsufficientFillets check, validation
    // succeeds (all specs reuse the same valid edge), then sphere-
    // fit math runs into the empty fillet_surfaces and surfaces as
    // InvalidInput("boundary") (compute_blend_boundary errors).
    let err =
        plan_apex_resolution(&brep, mvfs.vertex, &specs, &[], &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InvalidInput {
                field: "boundary",
                ..
            }
        ),
        "expected InvalidInput(boundary), got {err:?}"
    );
}

/// Helper for building a 3-fillet test BRep with an apex vertex
/// and 3 incident edges.  Returns (brep, apex_vertex, edges[3]).
///
/// Uses a unit cube and picks any corner vertex (valence-3) plus
/// its 3 incident edges.  This matches what `fillet_corner` does
/// in the existing kp_gate_helpers test infrastructure.
fn build_apex_with_three_edges() -> (
    BRep,
    popkern_core::topology::VertexHandle,
    [popkern_core::topology::EdgeHandle; 3],
) {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let apex = cube.vertices.keys().next().unwrap();
    // Collect the 3 edges incident to this corner.
    let mut edge_set: std::collections::BTreeSet<popkern_core::topology::EdgeHandle> =
        std::collections::BTreeSet::new();
    for he in cube.half_edges.values() {
        if he.origin == apex {
            edge_set.insert(he.edge);
        }
    }
    let edges: Vec<_> = edge_set.into_iter().collect();
    assert_eq!(edges.len(), 3, "cube corner should have 3 incident edges");
    (cube, apex, [edges[0], edges[1], edges[2]])
}

#[test]
fn rejects_nonexistent_apex_vertex() {
    let (brep, _apex, edges) = build_apex_with_three_edges();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    // Use a vertex handle that doesn't exist (default).
    let bogus_apex = popkern_core::topology::VertexHandle::default();
    let err =
        plan_apex_resolution(&brep, bogus_apex, &specs, &[], &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InvalidInput { field: "apex", .. }
        ),
        "expected InvalidInput(apex), got {err:?}"
    );
}

#[test]
fn rejects_nonexistent_edge() {
    let (brep, apex, _edges) = build_apex_with_three_edges();
    let specs = vec![
        FilletSpec {
            edge: popkern_core::topology::EdgeHandle::default(),
            radius: FilletRadius::Constant(0.1),
        };
        3
    ];
    let err = plan_apex_resolution(&brep, apex, &specs, &[], &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InvalidInput {
                field: "fillet_specs[].edge",
                ..
            }
        ),
        "expected InvalidInput(edge), got {err:?}"
    );
}

#[test]
fn rejects_edge_not_incident_to_apex() {
    let (brep, apex, edges) = build_apex_with_three_edges();
    // Find an edge in the cube NOT incident to apex.
    let apex_edges: std::collections::BTreeSet<_> = edges.iter().copied().collect();
    let far_edge = brep
        .edges
        .keys()
        .find(|e| !apex_edges.contains(e))
        .expect("cube has more than 3 edges");
    let specs = vec![
        FilletSpec {
            edge: edges[0],
            radius: FilletRadius::Constant(0.1),
        },
        FilletSpec {
            edge: edges[1],
            radius: FilletRadius::Constant(0.1),
        },
        FilletSpec {
            edge: far_edge, // not incident to apex
            radius: FilletRadius::Constant(0.1),
        },
    ];
    let err = plan_apex_resolution(&brep, apex, &specs, &[], &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InvalidInput {
                field: "fillet_specs[].edge",
                ..
            }
        ),
        "expected InvalidInput(edge not incident), got {err:?}"
    );
}

#[test]
fn rejects_nan_radius() {
    let (brep, apex, edges) = build_apex_with_three_edges();
    let specs = vec![
        FilletSpec {
            edge: edges[0],
            radius: FilletRadius::Constant(f64::NAN),
        },
        FilletSpec {
            edge: edges[1],
            radius: FilletRadius::Constant(0.1),
        },
        FilletSpec {
            edge: edges[2],
            radius: FilletRadius::Constant(0.1),
        },
    ];
    let err = plan_apex_resolution(&brep, apex, &specs, &[], &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InvalidInput {
                field: "fillet_specs[].radius",
                ..
            }
        ),
        "expected InvalidInput(NaN radius), got {err:?}"
    );
}

#[test]
fn rejects_zero_radius() {
    let (brep, apex, edges) = build_apex_with_three_edges();
    let specs = vec![
        FilletSpec {
            edge: edges[0],
            radius: FilletRadius::Constant(0.0),
        },
        FilletSpec {
            edge: edges[1],
            radius: FilletRadius::Constant(0.1),
        },
        FilletSpec {
            edge: edges[2],
            radius: FilletRadius::Constant(0.1),
        },
    ];
    let err = plan_apex_resolution(&brep, apex, &specs, &[], &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InvalidInput {
                field: "fillet_specs[].radius",
                ..
            }
        ),
        "expected InvalidInput(zero radius), got {err:?}"
    );
}

#[test]
fn rejects_negative_radius() {
    let (brep, apex, edges) = build_apex_with_three_edges();
    let specs = vec![
        FilletSpec {
            edge: edges[0],
            radius: FilletRadius::Constant(-0.5),
        },
        FilletSpec {
            edge: edges[1],
            radius: FilletRadius::Constant(0.1),
        },
        FilletSpec {
            edge: edges[2],
            radius: FilletRadius::Constant(0.1),
        },
    ];
    let err = plan_apex_resolution(&brep, apex, &specs, &[], &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InvalidInput {
                field: "fillet_specs[].radius",
                ..
            }
        ),
        "expected InvalidInput(negative radius), got {err:?}"
    );
}

#[test]
fn rejects_inf_radius() {
    let (brep, apex, edges) = build_apex_with_three_edges();
    let specs = vec![
        FilletSpec {
            edge: edges[0],
            radius: FilletRadius::Constant(f64::INFINITY),
        },
        FilletSpec {
            edge: edges[1],
            radius: FilletRadius::Constant(0.1),
        },
        FilletSpec {
            edge: edges[2],
            radius: FilletRadius::Constant(0.1),
        },
    ];
    let err = plan_apex_resolution(&brep, apex, &specs, &[], &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InvalidInput {
                field: "fillet_specs[].radius",
                ..
            }
        ),
        "expected InvalidInput(infinite radius), got {err:?}"
    );
}

#[test]
fn rejects_variable_radius_with_nan_start() {
    let (brep, apex, edges) = build_apex_with_three_edges();
    let specs = vec![
        FilletSpec {
            edge: edges[0],
            radius: FilletRadius::Variable {
                start: f64::NAN,
                end: 0.1,
            },
        },
        FilletSpec {
            edge: edges[1],
            radius: FilletRadius::Constant(0.1),
        },
        FilletSpec {
            edge: edges[2],
            radius: FilletRadius::Constant(0.1),
        },
    ];
    let err = plan_apex_resolution(&brep, apex, &specs, &[], &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InvalidInput {
                field: "fillet_specs[].radius",
                ..
            }
        ),
        "expected InvalidInput(variable NaN), got {err:?}"
    );
}

#[test]
fn rejects_bypoints_radius_empty() {
    let (brep, apex, edges) = build_apex_with_three_edges();
    let specs = vec![
        FilletSpec {
            edge: edges[0],
            radius: FilletRadius::ByPoints(vec![]),
        },
        FilletSpec {
            edge: edges[1],
            radius: FilletRadius::Constant(0.1),
        },
        FilletSpec {
            edge: edges[2],
            radius: FilletRadius::Constant(0.1),
        },
    ];
    let err = plan_apex_resolution(&brep, apex, &specs, &[], &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InvalidInput {
                field: "fillet_specs[].radius",
                ..
            }
        ),
        "expected InvalidInput(empty ByPoints), got {err:?}"
    );
}

#[test]
fn rejects_bypoints_t_out_of_range() {
    let (brep, apex, edges) = build_apex_with_three_edges();
    let specs = vec![
        FilletSpec {
            edge: edges[0],
            radius: FilletRadius::ByPoints(vec![(2.0, 0.1)]),
        },
        FilletSpec {
            edge: edges[1],
            radius: FilletRadius::Constant(0.1),
        },
        FilletSpec {
            edge: edges[2],
            radius: FilletRadius::Constant(0.1),
        },
    ];
    let err = plan_apex_resolution(&brep, apex, &specs, &[], &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InvalidInput {
                field: "fillet_specs[].radius",
                ..
            }
        ),
        "expected InvalidInput(ByPoints t out of range), got {err:?}"
    );
}

/// Direct unit tests for `decide_cap_kind` — invokes the logic
/// without going through the full planner.  Constructs synthetic
/// `PerFilletTermination` arrays and verifies cap-kind dispatch.
mod decide_cap_kind_tests {
    use super::super::cap_kind::CapKind;
    use super::super::plan::{PerFilletTermination, decide_cap_kind};
    use super::super::telemetry::PlanTelemetry;
    use crate::types::{FilletRadius, FilletSpec};
    use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::{EdgeHandle, FaceHandle};

    fn dummy_curve() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        }
    }

    fn term(
        idx: usize,
        endpoint_a: Point3,
        endpoint_b: Point3,
        cap_centre: Point3,
        cap_radius: f64,
    ) -> PerFilletTermination {
        PerFilletTermination {
            fillet_idx: idx,
            endpoint_a,
            endpoint_b,
            cap_centre,
            cap_radius,
            curve_3d: dummy_curve(),
            face_a_normal: Vec3::new(1.0, 0.0, 0.0),
            face_b_normal: Vec3::new(0.0, 1.0, 0.0),
            face_a: FaceHandle::default(),
            face_b: FaceHandle::default(),
        }
    }

    fn dummy_specs(n: usize, r: f64) -> Vec<FilletSpec> {
        (0..n)
            .map(|_| FilletSpec {
                edge: EdgeHandle::default(),
                radius: FilletRadius::Constant(r),
            })
            .collect()
    }

    #[test]
    fn dispatches_sphere_for_uniform_radii_and_collapsed_bridges() {
        // Three fillets at uniform radius 0.1, bridges collapse
        // (T_b[k] == T_a[k+1] within tol).
        let p0 = Point3::new(0.1, 0.0, 0.1);
        let p1 = Point3::new(0.1, 0.1, 0.0);
        let p2 = Point3::new(0.0, 0.1, 0.1);
        // Fillet 0: T_a = p2, T_b = p0 (so T_b[0] = T_a[1] = p0).
        // Fillet 1: T_a = p0, T_b = p1.
        // Fillet 2: T_a = p1, T_b = p2.
        let centre = Point3::new(0.0, 0.0, 0.0);
        let terms = vec![
            term(0, p2, p0, centre, 0.1),
            term(1, p0, p1, centre, 0.1),
            term(2, p1, p2, centre, 0.1),
        ];
        let specs = dummy_specs(3, 0.1);
        let mut telemetry = PlanTelemetry::new();
        let kind = decide_cap_kind(&terms, &specs, &Tolerance::default(), &mut telemetry)
            .expect("decide_cap_kind should succeed");
        assert!(
            matches!(kind, CapKind::Sphere { .. }),
            "expected Sphere, got {kind:?}"
        );
        assert!(telemetry.cap_kind_reason.contains("Sphere"));
    }

    #[test]
    fn dispatches_gregory_for_variable_radii() {
        // Three fillets, radii 0.08 / 0.10 / 0.12 (not uniform).
        let p0 = Point3::new(0.1, 0.0, 0.1);
        let p1 = Point3::new(0.1, 0.1, 0.0);
        let p2 = Point3::new(0.0, 0.1, 0.1);
        let centre = Point3::new(0.0, 0.0, 0.0);
        let terms = vec![
            term(0, p2, p0, centre, 0.08),
            term(1, p0, p1, centre, 0.10),
            term(2, p1, p2, centre, 0.12),
        ];
        let specs = dummy_specs(3, 0.10); // ignored by the function
        let mut telemetry = PlanTelemetry::new();
        let kind = decide_cap_kind(&terms, &specs, &Tolerance::default(), &mut telemetry)
            .expect("decide_cap_kind should succeed");
        assert!(
            matches!(kind, CapKind::Gregory { .. }),
            "expected Gregory, got {kind:?}"
        );
        assert!(telemetry.cap_kind_reason.contains("variable radii"));
    }

    #[test]
    fn dispatches_gregory_for_uniform_radii_but_disjoint_bridges() {
        // Uniform radii but T_b[k] != T_a[k+1] (bridges don't collapse).
        let p0a = Point3::new(0.1, 0.0, 0.1);
        let p0b = Point3::new(0.1, 0.05, 0.1);
        let p1a = Point3::new(0.1, 0.1, 0.0);
        let p1b = Point3::new(0.05, 0.1, 0.05);
        let p2a = Point3::new(0.0, 0.1, 0.1);
        let p2b = Point3::new(0.05, 0.05, 0.1);
        let centre = Point3::new(0.0, 0.0, 0.0);
        let terms = vec![
            term(0, p0a, p0b, centre, 0.1),
            term(1, p1a, p1b, centre, 0.1),
            term(2, p2a, p2b, centre, 0.1),
        ];
        let specs = dummy_specs(3, 0.1);
        let mut telemetry = PlanTelemetry::new();
        let kind = decide_cap_kind(&terms, &specs, &Tolerance::default(), &mut telemetry)
            .expect("decide_cap_kind should succeed");
        assert!(
            matches!(kind, CapKind::Gregory { .. }),
            "expected Gregory (bridges open), got {kind:?}"
        );
    }

    #[test]
    fn telemetry_records_diagnostic_values() {
        let p0 = Point3::new(0.1, 0.0, 0.1);
        let p1 = Point3::new(0.1, 0.1, 0.0);
        let p2 = Point3::new(0.0, 0.1, 0.1);
        let centre = Point3::new(0.0, 0.0, 0.0);
        let terms = vec![
            term(0, p2, p0, centre, 0.08),
            term(1, p0, p1, centre, 0.10),
            term(2, p1, p2, centre, 0.12),
        ];
        let specs = dummy_specs(3, 0.10);
        let mut telemetry = PlanTelemetry::new();
        let _ = decide_cap_kind(&terms, &specs, &Tolerance::default(), &mut telemetry).unwrap();
        assert!(telemetry.diagnostic_values.contains_key("radii_min"));
        assert!(telemetry.diagnostic_values.contains_key("radii_max"));
        assert!(telemetry.diagnostic_values.contains_key("radii_spread"));
        assert!(
            telemetry
                .diagnostic_values
                .contains_key("max_bridge_collapse")
        );
        assert!(
            telemetry
                .diagnostic_values
                .contains_key("bridge_collapse_threshold")
        );
    }
}

#[test]
fn passes_validation_for_valid_three_edge_apex() {
    let (brep, apex, edges) = build_apex_with_three_edges();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    // No fillet_surfaces provided ⇒ compute_blend_boundary errors,
    // surfacing as InvalidInput("boundary").  Validation has passed
    // by this point, so the failure mode confirms we got past
    // input validation and entered the sphere-fit math (A2.1.c).
    let err = plan_apex_resolution(&brep, apex, &specs, &[], &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InvalidInput {
                field: "boundary",
                ..
            }
        ),
        "expected InvalidInput(boundary) — validation passed, sphere-fit reached, got {err:?}"
    );
}
