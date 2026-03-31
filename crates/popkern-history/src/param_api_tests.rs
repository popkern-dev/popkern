// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
#[cfg(test)]
use crate::FeatureTree;
use crate::expression::Expression;
use crate::feature::{Feature, SketchPlaneDesc};
use crate::unit::Unit;

fn sketch_xy() -> Feature {
    Feature::Sketch {
        plane: SketchPlaneDesc::XY,
        sketch_data: String::new(),
    }
}

// ---------------------------------------------------------------------------
// 1. set_parameter_expression → undo → expression cleared
// ---------------------------------------------------------------------------

#[test]
fn set_expression_undo_restores_none() {
    let mut tree = FeatureTree::new();
    let id = tree
        .add_global_parameter("width", 10.0, Unit::Millimeters)
        .expect("add param");
    // Clear undo log so we start fresh.
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    let expr = Expression::Literal(42.0);
    tree.set_parameter_expression(id, expr)
        .expect("set expression");

    assert!(
        tree.parameters().get(id).unwrap().expression.is_some(),
        "expression should be set"
    );

    tree.undo();

    assert!(
        tree.parameters().get(id).unwrap().expression.is_none(),
        "undo should clear the expression"
    );
}

#[test]
fn set_expression_undo_restores_old_expression() {
    let mut tree = FeatureTree::new();
    let id = tree
        .add_global_parameter("depth", 5.0, Unit::Millimeters)
        .expect("add param");
    // Set an initial expression.
    let first_expr = Expression::Literal(5.0);
    tree.set_parameter_expression(id, first_expr.clone())
        .expect("set first expression");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    let second_expr = Expression::Literal(99.0);
    tree.set_parameter_expression(id, second_expr)
        .expect("set second expression");
    tree.undo();

    // After undo the expression should be restored (Some, not None).
    assert!(
        tree.parameters().get(id).unwrap().expression.is_some(),
        "undo should restore the first expression"
    );
}

#[test]
fn set_expression_redo_reapplies() {
    let mut tree = FeatureTree::new();
    let id = tree
        .add_global_parameter("height", 20.0, Unit::Millimeters)
        .expect("add param");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    let expr = Expression::Literal(77.0);
    tree.set_parameter_expression(id, expr)
        .expect("set expression");
    tree.undo();
    assert!(tree.parameters().get(id).unwrap().expression.is_none());

    tree.redo();
    assert!(
        tree.parameters().get(id).unwrap().expression.is_some(),
        "redo should reapply the expression"
    );
}

// ---------------------------------------------------------------------------
// 2. set_parameter_value → undo → old value restored
// ---------------------------------------------------------------------------

#[test]
fn set_parameter_value_undo_restores_old() {
    let mut tree = FeatureTree::new();
    let id = tree
        .add_global_parameter("radius", 3.0, Unit::Meters)
        .expect("add param");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    tree.set_parameter_value(id, 9.0, Unit::Meters)
        .expect("set value");
    assert_eq!(tree.parameters().get(id).unwrap().value, 9.0);

    tree.undo();
    assert_eq!(
        tree.parameters().get(id).unwrap().value,
        3.0,
        "undo should restore old value"
    );
}

#[test]
fn set_parameter_value_redo_reapplies() {
    let mut tree = FeatureTree::new();
    let id = tree
        .add_global_parameter("len", 1.0, Unit::Meters)
        .expect("add param");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    tree.set_parameter_value(id, 5.0, Unit::Meters)
        .expect("set value");
    tree.undo();
    assert_eq!(tree.parameters().get(id).unwrap().value, 1.0);

    tree.redo();
    assert_eq!(
        tree.parameters().get(id).unwrap().value,
        5.0,
        "redo should reapply new value"
    );
}

// ---------------------------------------------------------------------------
// 3. set_parameter_bounds → undo → old bounds restored
// ---------------------------------------------------------------------------

#[test]
fn set_bounds_undo_restores_none() {
    let mut tree = FeatureTree::new();
    let id = tree
        .add_global_parameter("angle", 45.0, Unit::Degrees)
        .expect("add param");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    tree.set_parameter_bounds(id, 0.0, 90.0);
    assert_eq!(tree.parameters().get(id).unwrap().bounds, Some((0.0, 90.0)));

    tree.undo();
    assert_eq!(
        tree.parameters().get(id).unwrap().bounds,
        None,
        "undo should clear bounds"
    );
}

#[test]
fn clear_bounds_undo_restores_old_bounds() {
    let mut tree = FeatureTree::new();
    let id = tree
        .add_global_parameter("angle2", 30.0, Unit::Degrees)
        .expect("add param");
    // Set bounds first.
    tree.parameters_mut().set_bounds(id, 10.0, 80.0);
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    tree.clear_parameter_bounds(id);
    assert!(tree.parameters().get(id).unwrap().bounds.is_none());

    tree.undo();
    assert_eq!(
        tree.parameters().get(id).unwrap().bounds,
        Some((10.0, 80.0)),
        "undo should restore old bounds"
    );
}

// ---------------------------------------------------------------------------
// 4. add_global_parameter → undo → unreachable by name → redo → back
// ---------------------------------------------------------------------------

#[test]
fn add_global_parameter_undo_removes_from_name_index() {
    let mut tree = FeatureTree::new();
    let id = tree
        .add_global_parameter("thickness", 2.0, Unit::Millimeters)
        .expect("add param");

    assert!(
        tree.parameters().resolve_name("thickness").is_some(),
        "should be in name index after add"
    );

    tree.undo();

    assert!(
        tree.parameters().resolve_name("thickness").is_none(),
        "undo should remove from name index"
    );

    // Redo should re-register.
    tree.redo();
    assert_eq!(
        tree.parameters().resolve_name("thickness"),
        Some(id),
        "redo should restore name mapping"
    );
}

#[test]
fn add_global_parameter_multiple_undo_redo_cycles() {
    let mut tree = FeatureTree::new();
    tree.add_global_parameter("p1", 1.0, Unit::Meters)
        .expect("add p1");
    tree.add_global_parameter("p2", 2.0, Unit::Meters)
        .expect("add p2");

    // Undo p2.
    tree.undo();
    assert!(tree.parameters().resolve_name("p2").is_none());
    assert!(tree.parameters().resolve_name("p1").is_some());

    // Undo p1.
    tree.undo();
    assert!(tree.parameters().resolve_name("p1").is_none());

    // Redo p1.
    tree.redo();
    assert!(tree.parameters().resolve_name("p1").is_some());

    // Redo p2.
    tree.redo();
    assert!(tree.parameters().resolve_name("p2").is_some());
}

// ---------------------------------------------------------------------------
// 5. rename_parameter → undo → old name restored
// ---------------------------------------------------------------------------

#[test]
fn rename_parameter_undo_restores_old_name() {
    let mut tree = FeatureTree::new();
    let id = tree
        .add_global_parameter("old_name", 7.0, Unit::Meters)
        .expect("add param");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    tree.rename_parameter(id, "new_name").expect("rename");
    assert!(tree.parameters().resolve_name("new_name").is_some());
    assert!(tree.parameters().resolve_name("old_name").is_none());

    tree.undo();
    assert!(
        tree.parameters().resolve_name("old_name").is_some(),
        "undo should restore old name"
    );
    assert!(
        tree.parameters().resolve_name("new_name").is_none(),
        "undo should remove new name"
    );
}

#[test]
fn rename_parameter_redo_reapplies() {
    let mut tree = FeatureTree::new();
    let id = tree
        .add_global_parameter("aaa", 1.0, Unit::Meters)
        .expect("add param");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    tree.rename_parameter(id, "bbb").expect("rename");
    tree.undo();
    tree.redo();

    assert!(
        tree.parameters().resolve_name("bbb").is_some(),
        "redo should restore renamed name"
    );
}

// ---------------------------------------------------------------------------
// 6. set_suppress_if → undo → old condition restored
// ---------------------------------------------------------------------------

#[test]
fn set_suppress_if_undo_restores_none() {
    let mut tree = FeatureTree::new();
    let fid = tree.add_feature(sketch_xy(), "S1");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    let expr = Expression::Literal(1.0);
    tree.set_suppress_if(fid, Some(expr));
    assert!(tree.feature(fid).unwrap().suppress_if.is_some());

    tree.undo();
    assert!(
        tree.feature(fid).unwrap().suppress_if.is_none(),
        "undo should clear suppress_if"
    );
}

#[test]
fn set_suppress_if_redo_reapplies() {
    let mut tree = FeatureTree::new();
    let fid = tree.add_feature(sketch_xy(), "S1");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    let expr = Expression::Literal(0.5);
    tree.set_suppress_if(fid, Some(expr));
    tree.undo();
    tree.redo();

    assert!(
        tree.feature(fid).unwrap().suppress_if.is_some(),
        "redo should restore suppress_if"
    );
}

// ---------------------------------------------------------------------------
// 7. remove_global_parameter → undo → name restored
// ---------------------------------------------------------------------------

#[test]
fn remove_global_parameter_undo_restores_name() {
    let mut tree = FeatureTree::new();
    let id = tree
        .add_global_parameter("removable", 5.0, Unit::Meters)
        .expect("add param");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    tree.remove_global_parameter(id).expect("remove param");
    assert!(tree.parameters().resolve_name("removable").is_none());

    tree.undo();
    assert!(
        tree.parameters().resolve_name("removable").is_some(),
        "undo should restore parameter name"
    );
}

// ---------------------------------------------------------------------------
// 8. clear_parameter_expression → undo → expression restored
// ---------------------------------------------------------------------------

#[test]
fn clear_expression_undo_restores_expression() {
    let mut tree = FeatureTree::new();
    let id = tree
        .add_global_parameter("x", 1.0, Unit::Meters)
        .expect("add param");
    // Use a simple literal expression so we avoid PartialEq on Expression.
    let expr = Expression::Literal(42.0);
    tree.parameters_mut()
        .set_expression(id, expr)
        .expect("set expression");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    assert!(tree.parameters().get(id).unwrap().expression.is_some());
    tree.clear_parameter_expression(id);
    assert!(tree.parameters().get(id).unwrap().expression.is_none());

    tree.undo();
    assert!(
        tree.parameters().get(id).unwrap().expression.is_some(),
        "undo should restore expression"
    );
}
