//! End-to-end integration tests for MEASURE expressions.
//!
//! Verifies that `Expression::Measure(MeasureQuery)` expressions evaluate
//! correctly through the full stack: expression evaluator → measure dispatch →
//! live BRep geometry.

use std::collections::HashMap;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_history::expression::evaluate;
use popkern_history::node::BodyId;
use popkern_history::{
    BinOp, BoundExtent, Expression, MeasureContext, MeasureQuery, ParameterTable,
};

// -- Helpers ------------------------------------------------------------------

fn unit_cube_ctx() -> (BRep, HashMap<BodyId, BRep>, Tolerance) {
    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).expect("unit cube should build");
    let mut bodies = HashMap::new();
    bodies.insert(BodyId(0), brep.clone());
    (brep, bodies, tol)
}

fn make_ctx<'a>(bodies: &'a HashMap<BodyId, BRep>, tol: &'a Tolerance) -> MeasureContext<'a> {
    MeasureContext {
        bodies,
        active_body: BodyId(0),
        tolerance: tol,
    }
}

// -- Tests --------------------------------------------------------------------

#[test]
fn measure_body_volume_of_unit_cube() {
    let (_brep, bodies, tol) = unit_cube_ctx();
    let ctx = make_ctx(&bodies, &tol);
    let table = ParameterTable::new();

    let expr = Expression::Measure(MeasureQuery::BodyVolume);
    let result = evaluate(&expr, &table, &ctx).expect("volume should evaluate");

    assert!(
        (result - 1.0).abs() < 0.01,
        "expected volume ≈ 1.0, got {result}"
    );
}

#[test]
fn measure_body_surface_area_of_unit_cube() {
    let (_brep, bodies, tol) = unit_cube_ctx();
    let ctx = make_ctx(&bodies, &tol);
    let table = ParameterTable::new();

    let expr = Expression::Measure(MeasureQuery::BodySurfaceArea);
    let result = evaluate(&expr, &table, &ctx).expect("surface area should evaluate");

    assert!(
        (result - 6.0).abs() < 0.01,
        "expected surface area ≈ 6.0, got {result}"
    );
}

#[test]
fn measure_body_bounding_box() {
    let (_brep, bodies, tol) = unit_cube_ctx();
    let ctx = make_ctx(&bodies, &tol);
    let table = ParameterTable::new();

    let expr = Expression::Measure(MeasureQuery::BodyBound(BoundExtent::SizeX));
    let result = evaluate(&expr, &table, &ctx).expect("bounding box size X should evaluate");

    assert!(
        (result - 1.0).abs() < 0.01,
        "expected bounding box SizeX ≈ 1.0, got {result}"
    );
}

#[test]
fn measure_in_arithmetic_expression() {
    let (_brep, bodies, tol) = unit_cube_ctx();
    let ctx = make_ctx(&bodies, &tol);
    let table = ParameterTable::new();

    // MEASURE(BodyVolume) * 2.0 — should equal 2.0
    let expr = Expression::BinOp(
        Box::new(Expression::Measure(MeasureQuery::BodyVolume)),
        BinOp::Mul,
        Box::new(Expression::Literal(2.0)),
    );
    let result = evaluate(&expr, &table, &ctx).expect("arithmetic expression should evaluate");

    assert!(
        (result - 2.0).abs() < 0.01,
        "expected volume * 2.0 ≈ 2.0, got {result}"
    );
}

#[test]
fn measure_without_context_returns_error() {
    let ctx = MeasureContext::empty();
    let table = ParameterTable::new();

    let expr = Expression::Measure(MeasureQuery::BodyVolume);
    let result = evaluate(&expr, &table, &ctx);

    assert!(
        result.is_err(),
        "expected error when no BRep is available, got Ok({:?})",
        result.ok()
    );
}

#[test]
fn measure_edge_length_of_unit_cube() {
    // Test edge_length via popkern_measure directly — the unit cube does not
    // have stable EntityNames assigned on its edges so we can't use the
    // MeasureQuery::EdgeLength path (which requires name resolution).
    // Body-level MEASURE queries cover the name-resolution path adequately.
    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).expect("unit cube should build");

    // Collect all edge handles.
    let edge_handles: Vec<_> = brep.edges.keys().collect();
    assert!(!edge_handles.is_empty(), "unit cube must have edges");

    // Every edge on a unit cube has length 1.0.
    for handle in edge_handles {
        let len = popkern_measure::edge_length(&brep, handle).expect("edge_length should succeed");
        assert!(
            (len - 1.0).abs() < 0.01,
            "expected unit cube edge length ≈ 1.0, got {len}"
        );
    }
}
