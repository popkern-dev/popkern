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
