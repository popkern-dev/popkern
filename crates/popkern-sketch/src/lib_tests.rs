//! Tests for SketchOnPlane: assign_element_names and sketch creation pipeline.

use popkern_constraint::ElementHandle;
use popkern_core::naming::{OperationId, RoleCategory};

use crate::SketchOnPlane;
use crate::plane::SketchPlane;

// ---------------------------------------------------------------------------
// assign_element_names with empty sketch
// ---------------------------------------------------------------------------

#[test]
fn assign_element_names_empty_sketch() {
    let mut sop = SketchOnPlane::new(SketchPlane::xy());
    let op = OperationId(42);

    sop.assign_element_names(op);

    assert!(
        sop.element_names.is_empty(),
        "empty sketch should produce no element names"
    );
}

// ---------------------------------------------------------------------------
// assign_element_names with a single line
// ---------------------------------------------------------------------------

#[test]
fn assign_element_names_single_line() {
    let mut sop = SketchOnPlane::new(SketchPlane::xy());
    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(1.0, 0.0);
    let _l = sop.sketch.add_line(p0, p1);

    let op = OperationId(7);
    sop.assign_element_names(op);

    assert_eq!(
        sop.element_names.len(),
        1,
        "one line should produce one element name"
    );

    // The single entry should have index 0 and Generated category.
    let name = sop.element_names.values().next().expect("has one entry");
    assert_eq!(name.operation, op);
    assert_eq!(name.role.category, RoleCategory::Generated);
    assert_eq!(name.role.index, 0);
    assert!(name.ancestors.is_empty());
}

// ---------------------------------------------------------------------------
// assign_element_names with multiple elements
// ---------------------------------------------------------------------------

#[test]
fn assign_element_names_multiple_elements() {
    let mut sop = SketchOnPlane::new(SketchPlane::xy());
    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(1.0, 0.0);
    let p2 = sop.sketch.add_point(1.0, 1.0);
    let p3 = sop.sketch.add_point(0.0, 1.0);
    let _l0 = sop.sketch.add_line(p0, p1);
    let _l1 = sop.sketch.add_line(p1, p2);
    let _l2 = sop.sketch.add_line(p2, p3);
    let _l3 = sop.sketch.add_line(p3, p0);

    let op = OperationId(100);
    sop.assign_element_names(op);

    assert_eq!(
        sop.element_names.len(),
        4,
        "four lines should produce four names"
    );

    // Indices should be 0..3 (in some order).
    let mut indices: Vec<u32> = sop.element_names.values().map(|n| n.role.index).collect();
    indices.sort();
    assert_eq!(indices, vec![0, 1, 2, 3]);

    // All should share the same operation and category.
    for name in sop.element_names.values() {
        assert_eq!(name.operation, op);
        assert_eq!(name.role.category, RoleCategory::Generated);
    }
}

// ---------------------------------------------------------------------------
// assign_element_names excludes construction geometry
// ---------------------------------------------------------------------------

#[test]
fn assign_element_names_excludes_construction() {
    let mut sop = SketchOnPlane::new(SketchPlane::xy());
    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(1.0, 0.0);
    let p2 = sop.sketch.add_point(1.0, 1.0);
    let l0 = sop.sketch.add_line(p0, p1);
    let _l1 = sop.sketch.add_line(p1, p2);

    // Mark l0 as construction.
    sop.sketch.set_construction(ElementHandle::Line(l0), true);

    let op = OperationId(55);
    sop.assign_element_names(op);

    assert_eq!(
        sop.element_names.len(),
        1,
        "construction line should be excluded"
    );
}

// ---------------------------------------------------------------------------
// assign_element_names clears previous names on re-call
// ---------------------------------------------------------------------------

#[test]
fn assign_element_names_clears_on_recall() {
    let mut sop = SketchOnPlane::new(SketchPlane::xy());
    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(1.0, 0.0);
    let _l = sop.sketch.add_line(p0, p1);

    sop.assign_element_names(OperationId(1));
    assert_eq!(sop.element_names.len(), 1);

    // Call again with a different operation; old names should be gone.
    sop.assign_element_names(OperationId(2));
    assert_eq!(sop.element_names.len(), 1);

    let name = sop.element_names.values().next().expect("has one entry");
    assert_eq!(
        name.operation,
        OperationId(2),
        "second call should overwrite with new operation id"
    );
}

// ---------------------------------------------------------------------------
// assign_element_names with arcs
// ---------------------------------------------------------------------------

#[test]
fn assign_element_names_includes_arcs() {
    use std::f64::consts::FRAC_PI_2;

    let mut sop = SketchOnPlane::new(SketchPlane::xy());
    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(1.0, 0.0);
    let _l = sop.sketch.add_line(p0, p1);
    let _arc = sop.sketch.add_arc(0.0, 0.0, 1.0, 0.0, FRAC_PI_2);

    let op = OperationId(99);
    sop.assign_element_names(op);

    // 1 line + 1 arc = 2 names.
    assert_eq!(sop.element_names.len(), 2, "line + arc should give 2 names");
}

// ---------------------------------------------------------------------------
// SketchOnPlane::new creates empty state
// ---------------------------------------------------------------------------

#[test]
fn new_sketch_on_plane_is_empty() {
    let sop = SketchOnPlane::new(SketchPlane::yz());

    assert!(sop.element_names.is_empty());
    assert_eq!(sop.sketch.non_construction_curves().len(), 0);
}
