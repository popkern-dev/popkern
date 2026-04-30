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
    let err = plan_apex_resolution(&brep, mvfs.vertex, &[], &Tolerance::default()).unwrap_err();
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
    let err =
        plan_apex_resolution(&brep, mvfs.vertex, &dummy_specs, &Tolerance::default()).unwrap_err();
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
    let err = plan_apex_resolution(&brep, mvfs.vertex, &specs, &Tolerance::default()).unwrap_err();
    assert!(
        matches!(
            err,
            super::error::ApexPlanError::InvalidInput {
                field: "planner",
                ..
            }
        ),
        "expected InvalidInput(planner stub), got {err:?}"
    );
}
