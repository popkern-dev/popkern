use crate::geometry::pcurve::Pcurve;
use crate::tolerance::Tolerance;
use crate::topology::{
    BRep, EdgeHandle, HalfEdgeHandle, SharedTrimEdgeFailureReason, SharedTrimEdgeInput,
    commit_shared_trim_edge_transactionally,
};

fn cube_edge_pair(brep: &BRep) -> (EdgeHandle, HalfEdgeHandle, HalfEdgeHandle) {
    let (edge, edge_data) = brep.edges.iter().next().expect("cube edge");
    let half_edge_a = edge_data.half_edge;
    let half_edge_b = brep.half_edges[half_edge_a].twin;
    (edge, half_edge_a, half_edge_b)
}

fn pcurve_midpoint_gap(a: &Pcurve, b: &Pcurve) -> f64 {
    let a_t = 0.5 * (a.t_start + a.t_end);
    let b_t = 0.5 * (b.t_start + b.t_end);
    a.point_at(a_t).distance_to(&b.point_at(b_t))
}

#[test]
fn shared_trim_edge_commits_valid_cube_edge_atomically() {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("unit cube");
    let (edge, half_edge_a, half_edge_b) = cube_edge_pair(&brep);
    let curve = brep.edges[edge].curve.clone().expect("edge curve");
    let pcurve_a = brep.half_edges[half_edge_a]
        .pcurve
        .clone()
        .expect("side a pcurve");
    let pcurve_b = brep.half_edges[half_edge_b]
        .pcurve
        .clone()
        .expect("side b pcurve");

    let commit = commit_shared_trim_edge_transactionally(
        &mut brep,
        SharedTrimEdgeInput {
            edge,
            half_edge_a,
            half_edge_b,
            curve,
            pcurve_a,
            pcurve_b,
            tolerance: Tolerance::default().linear(),
            provenance: "shared_trim_edge_test",
        },
    )
    .expect("valid cube edge commit");

    assert_eq!(commit.edge, edge);
    assert!(commit.max_endpoint_gap <= Tolerance::default().linear());
    assert!(commit.max_curve_surface_deviation <= Tolerance::default().linear());
    assert!(commit.max_side_agreement_deviation <= Tolerance::default().linear());
    brep.validate().expect("valid brep after commit");
}

#[test]
fn shared_trim_edge_rolls_back_invalid_side_pcurve() {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("unit cube");
    let before = brep.clone();
    let (edge, half_edge_a, half_edge_b) = cube_edge_pair(&brep);
    let curve = brep.edges[edge].curve.clone().expect("edge curve");
    let pcurve_a = brep.half_edges[half_edge_a]
        .pcurve
        .clone()
        .expect("side a pcurve");
    let existing_b = brep.half_edges[half_edge_b]
        .pcurve
        .clone()
        .expect("side b pcurve");
    let bad_pcurve_b = Pcurve::linear(
        (10.0, 10.0),
        (11.0, 10.0),
        existing_b.t_start,
        existing_b.t_end,
    );

    let err = commit_shared_trim_edge_transactionally(
        &mut brep,
        SharedTrimEdgeInput {
            edge,
            half_edge_a,
            half_edge_b,
            curve,
            pcurve_a,
            pcurve_b: bad_pcurve_b,
            tolerance: Tolerance::default().linear(),
            provenance: "shared_trim_edge_test",
        },
    )
    .expect_err("bad side pcurve must decline");

    assert_eq!(
        err.reason,
        SharedTrimEdgeFailureReason::CurveSurfaceMismatch
    );
    let after_b = brep.half_edges[half_edge_b]
        .pcurve
        .as_ref()
        .expect("rolled back side b");
    let before_b = before.half_edges[half_edge_b]
        .pcurve
        .as_ref()
        .expect("original side b");
    assert!(pcurve_midpoint_gap(after_b, before_b) <= Tolerance::default().numerical_floor());
    brep.validate().expect("valid brep after rollback");
}

#[test]
fn shared_trim_edge_rejects_non_twin_half_edges() {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("unit cube");
    let (edge, half_edge_a, _) = cube_edge_pair(&brep);
    let half_edge_b = brep
        .half_edges
        .iter()
        .map(|(handle, _)| handle)
        .find(|handle| *handle != half_edge_a && *handle != brep.half_edges[half_edge_a].twin)
        .expect("non twin half-edge");
    let curve = brep.edges[edge].curve.clone().expect("edge curve");
    let pcurve_a = brep.half_edges[half_edge_a]
        .pcurve
        .clone()
        .expect("side a pcurve");
    let pcurve_b = brep.half_edges[half_edge_b]
        .pcurve
        .clone()
        .expect("side b pcurve");

    let err = commit_shared_trim_edge_transactionally(
        &mut brep,
        SharedTrimEdgeInput {
            edge,
            half_edge_a,
            half_edge_b,
            curve,
            pcurve_a,
            pcurve_b,
            tolerance: Tolerance::default().linear(),
            provenance: "shared_trim_edge_test",
        },
    )
    .expect_err("non twin half-edge must decline");

    assert!(
        matches!(
            err.reason,
            SharedTrimEdgeFailureReason::HalfEdgeEdgeMismatch
                | SharedTrimEdgeFailureReason::HalfEdgesNotTwins
        ),
        "unexpected reason: {:?}",
        err.reason
    );
    brep.validate().expect("valid brep after decline");
}
