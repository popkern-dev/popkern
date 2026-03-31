//! Integration tests for the parametric expression pipeline in popkern-history
//! (tests 1–5).
//!
//! See `parametric_advanced.rs` for tests 6–10 (save/load, undo, delete,
//! bounds, built-ins).
//!
//! Convention: `add_global(name, value, unit)` stores `value` in SI base units
//! directly (no conversion). The `unit` field is display metadata only.
//! Op IDs are assigned only on `add_feature`, starting at 1.

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_history::*;

// ---------------------------------------------------------------------------
// Helpers (shared via pub fn — `parametric_advanced.rs` duplicates these)
// ---------------------------------------------------------------------------

pub fn rectangle_sketch_json() -> String {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0);
    let p3 = sk.add_point(0.0, 1.0);
    sk.add_line(p0, p1);
    sk.add_line(p1, p2);
    sk.add_line(p2, p3);
    sk.add_line(p3, p0);
    sk.add_constraint(Constraint::Fixed(p0));
    sk.add_constraint(Constraint::Fixed(p1));
    sk.add_constraint(Constraint::Fixed(p2));
    sk.add_constraint(Constraint::Fixed(p3));
    sk.solve().expect("rectangle should solve");
    serde_json::to_string(&sk).expect("rectangle should serialize")
}

pub fn sketch_face_name(op_id: u64) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: Vec::new(),
    }
}

/// Look up the op_id assigned to a feature.
pub fn feature_op_id(tree: &FeatureTree, id: FeatureId) -> u64 {
    tree.feature(id).expect("feature").operation_id.0
}

// ---------------------------------------------------------------------------
// Test 1: Parametric extrude with expression
// ---------------------------------------------------------------------------

#[test]
fn parametric_extrude_with_expression() {
    let mut tree = FeatureTree::new();

    // Global parameter: width = 0.05 SI (displayed as 50mm)
    let width = tree
        .add_global_parameter("width", 0.05, Unit::Millimeters)
        .expect("add width");

    // Sketch gets op_id = 1
    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch1",
    );
    let sketch_op = feature_op_id(&tree, sketch_id);

    // Extrude distance parameter, initially 0.0
    let dist = tree
        .add_global_parameter("_extrude1_dist", 0.0, Unit::Meters)
        .expect("add dist");

    // Set expression: distance = width * 2
    let expr = Expression::BinOp(
        Box::new(Expression::Param(width)),
        BinOp::Mul,
        Box::new(Expression::Literal(2.0)),
    );
    tree.set_parameter_expression(dist, expr)
        .expect("set expression");

    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: sketch_face_name(sketch_op),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist,
        },
        "Extrude1",
    );

    let result = tree.replay(Tolerance::default()).expect("replay");
    assert!(
        matches!(
            result.feature_results.get(&extrude_id),
            Some(FeatureResult::Success { .. })
        ),
        "extrude should succeed; got {:?}",
        result.feature_results.get(&extrude_id)
    );

    // After replay: dist should be 0.1m (0.05 * 2)
    let dist_val = tree.parameters().get(dist).expect("dist").value;
    assert!(
        (dist_val - 0.1).abs() < 1e-10,
        "dist should be 0.1m (0.05 * 2), got {dist_val}"
    );

    let brep1 = tree.current_brep().expect("brep").clone();
    let max_z_1: f64 = brep1
        .vertices
        .iter()
        .map(|(_, v)| v.position.z)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (max_z_1 - 0.1).abs() < 0.01,
        "first extrusion height should be 0.1m, got {max_z_1}"
    );

    // Change width to 0.1 SI (100mm) — extrude should now be 0.2m
    tree.set_parameter_value(width, 0.1, Unit::Meters)
        .expect("set width");
    tree.replay(Tolerance::default()).expect("second replay");

    let dist_val2 = tree.parameters().get(dist).expect("dist").value;
    assert!(
        (dist_val2 - 0.2).abs() < 1e-10,
        "dist should be 0.2m (0.1 * 2), got {dist_val2}"
    );

    let brep2 = tree.current_brep().expect("brep after change");
    let max_z_2: f64 = brep2
        .vertices
        .iter()
        .map(|(_, v)| v.position.z)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (max_z_2 - 0.2).abs() < 0.01,
        "second extrusion height should be 0.2m, got {max_z_2}"
    );
}

// ---------------------------------------------------------------------------
// Test 2: Expression chain across features
// ---------------------------------------------------------------------------

#[test]
fn expression_chain_across_features() {
    let mut tree = FeatureTree::new();

    // base = 0.01 SI (10mm)
    let base = tree
        .add_global_parameter("base", 0.01, Unit::Millimeters)
        .expect("add base");

    // fillet_ratio = 0.1 (dimensionless)
    let fillet_ratio = tree
        .add_global_parameter("fillet_ratio", 0.1, Unit::Ratio)
        .expect("add fillet_ratio");

    // Sketch gets op_id = 1
    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch1",
    );
    let sketch_op = feature_op_id(&tree, sketch_id);

    // Extrude distance = base (expression: Param(base))
    let extrude_dist = tree
        .add_global_parameter("_extrude_dist", 0.0, Unit::Meters)
        .expect("add extrude_dist");
    tree.set_parameter_expression(extrude_dist, Expression::Param(base))
        .expect("set extrude_dist expr");

    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: sketch_face_name(sketch_op),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: extrude_dist,
        },
        "Extrude1",
    );

    // Fillet radius = base * fillet_ratio
    let fillet_radius = tree
        .add_global_parameter("_fillet_radius", 0.0, Unit::Meters)
        .expect("add fillet_radius");
    let radius_expr = Expression::BinOp(
        Box::new(Expression::Param(base)),
        BinOp::Mul,
        Box::new(Expression::Param(fillet_ratio)),
    );
    tree.set_parameter_expression(fillet_radius, radius_expr)
        .expect("set fillet_radius expr");

    // Use a bogus edge name; fillet may fail but expression must propagate.
    let bogus_edge = EntityName {
        operation: OperationId(999),
        role: EntityRole {
            category: RoleCategory::Lateral,
            index: 0,
        },
        ancestors: Vec::new(),
    };
    let fillet_id = tree.add_feature(
        Feature::Fillet {
            edges: vec![bogus_edge],
            radius: fillet_radius,
        },
        "Fillet1",
    );

    tree.replay(Tolerance::default()).expect("replay");

    // Extrude should succeed
    assert!(
        matches!(
            tree.feature(extrude_id)
                .expect("extrude node")
                .last_result
                .as_ref(),
            Some(FeatureResult::Success { .. })
        ),
        "extrude should succeed; got {:?}",
        tree.feature(extrude_id).unwrap().last_result
    );

    // Fillet was attempted (Success or Failed — geometry may not support it)
    assert!(
        tree.feature(fillet_id)
            .expect("fillet node")
            .last_result
            .is_some(),
        "fillet should have been evaluated"
    );

    // Verify expression propagated: fillet_radius = 0.01 * 0.1 = 0.001m
    let fr_val = tree
        .parameters()
        .get(fillet_radius)
        .expect("fillet_radius")
        .value;
    assert!(
        (fr_val - 0.001).abs() < 1e-10,
        "fillet radius should be 0.001m (0.01 * 0.1), got {fr_val}"
    );
}

// ---------------------------------------------------------------------------
// Test 3: Conditional suppression end-to-end
// ---------------------------------------------------------------------------

#[test]
fn conditional_suppression_end_to_end() {
    let mut tree = FeatureTree::new();

    // Sketch gets op_id = 1
    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch1",
    );
    let sketch_op = feature_op_id(&tree, sketch_id);

    let dist = tree
        .add_global_parameter("_cse_dist", 0.1, Unit::Meters)
        .expect("add dist");

    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: sketch_face_name(sketch_op),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist,
        },
        "Extrude1",
    );
    tree.replay(Tolerance::default()).expect("initial replay");
    assert_eq!(tree.current_brep().expect("brep").faces.len(), 6);

    let fillet_radius = tree
        .add_global_parameter("_cse_fillet_r", 0.01, Unit::Meters)
        .expect("add fillet_r");

    // Find an actual edge from the extruded box to fillet.
    let edge_name = {
        let brep = tree.current_brep().expect("brep");
        let mut found = None;
        for (_h, edge) in &brep.edges {
            if let Some(nid) = edge.name
                && let Some(name) = brep.name_pool.get(nid)
            {
                found = Some(name.clone());
                break;
            }
        }
        found
    };

    let fillet_id = if let Some(en) = edge_name {
        tree.add_feature(
            Feature::Fillet {
                edges: vec![en],
                radius: fillet_radius,
            },
            "Fillet1",
        )
    } else {
        let bogus = EntityName {
            operation: OperationId(999),
            role: EntityRole {
                category: RoleCategory::Lateral,
                index: 0,
            },
            ancestors: Vec::new(),
        };
        tree.add_feature(
            Feature::Fillet {
                edges: vec![bogus],
                radius: fillet_radius,
            },
            "Fillet1",
        )
    };

    // suppress_if = Literal(1.0) → always suppressed
    tree.set_suppress_if(fillet_id, Some(Expression::Literal(1.0)))
        .expect("set suppress_if");
    let result = tree
        .replay(Tolerance::default())
        .expect("replay suppressed");

    assert!(
        matches!(
            result.feature_results.get(&fillet_id),
            Some(FeatureResult::Suppressed)
        ),
        "fillet should be suppressed when suppress_if = 1.0, got {:?}",
        result.feature_results.get(&fillet_id)
    );

    // Change suppress_if to Literal(0.0) → fillet should be attempted
    tree.set_suppress_if(fillet_id, Some(Expression::Literal(0.0)))
        .expect("set suppress_if");
    let result2 = tree
        .replay(Tolerance::default())
        .expect("replay unsuppressed");

    let fillet_result = result2.feature_results.get(&fillet_id);
    assert!(
        matches!(
            fillet_result,
            Some(FeatureResult::Success { .. }) | Some(FeatureResult::Failed { .. })
        ),
        "fillet should be attempted when suppress_if = 0.0, got {:?}",
        fillet_result
    );

    // Extrude should still succeed
    assert!(
        matches!(
            result2.feature_results.get(&extrude_id),
            Some(FeatureResult::Success { .. })
        ),
        "extrude should still succeed"
    );
}

// ---------------------------------------------------------------------------
// Test 4: Expression parser integration
// ---------------------------------------------------------------------------

#[test]
fn expression_parser_integration() {
    let mut tree = FeatureTree::new();

    // Store SI values: width = 0.05m, height = 0.1m
    tree.add_global_parameter("width", 0.05, Unit::Millimeters)
        .expect("add width");
    tree.add_global_parameter("height", 0.1, Unit::Millimeters)
        .expect("add height");

    // Parse "width * 2 + height"
    let expr =
        parse_expression("width * 2 + height", tree.parameters()).expect("parse should succeed");

    // Evaluate: 0.05 * 2 + 0.1 = 0.2
    let result = {
        use popkern_history::MeasureContext;
        use popkern_history::expression::evaluate;
        evaluate(&expr, tree.parameters(), &MeasureContext::empty()).expect("evaluate")
    };
    assert!(
        (result - 0.2).abs() < 1e-10,
        "width(0.05) * 2 + height(0.1) should be 0.2, got {result}"
    );
}

// ---------------------------------------------------------------------------
// Test 5: Expression formatter integration with round-trip
// ---------------------------------------------------------------------------

#[test]
fn expression_formatter_round_trip() {
    let mut tree = FeatureTree::new();

    let width = tree
        .add_global_parameter("width", 0.05, Unit::Millimeters)
        .expect("add width");

    // Build: width > 0.1 ? 0.005 : 0.003
    let cond = Expression::Compare(
        Box::new(Expression::Param(width)),
        CmpOp::Gt,
        Box::new(Expression::Literal(0.1)),
    );
    let then_expr = Expression::Literal(0.005);
    let else_expr = Expression::Literal(0.003);
    let original =
        Expression::Conditional(Box::new(cond), Box::new(then_expr), Box::new(else_expr));

    let formatted = format_expression(&original, tree.parameters());
    assert!(
        !formatted.is_empty(),
        "formatted expression should not be empty"
    );

    // Parse the formatted string back
    let reparsed = parse_expression(&formatted, tree.parameters())
        .expect("re-parse of formatted expression should succeed");

    // Evaluate both and compare
    let eval_orig = {
        use popkern_history::MeasureContext;
        use popkern_history::expression::evaluate;
        evaluate(&original, tree.parameters(), &MeasureContext::empty()).expect("eval original")
    };
    let eval_reparsed = {
        use popkern_history::MeasureContext;
        use popkern_history::expression::evaluate;
        evaluate(&reparsed, tree.parameters(), &MeasureContext::empty()).expect("eval reparsed")
    };
    assert!(
        (eval_orig - eval_reparsed).abs() < 1e-10,
        "round-tripped expression should evaluate to same value: orig={eval_orig}, reparsed={eval_reparsed}"
    );
}
