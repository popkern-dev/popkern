//! Tests for the cellular decomposition.

use popkern_core::geometry::Point3;

use super::cell::{Cell, CellBoundaryFace, CellDecomposition};
use super::select::select_for_operation;
use crate::error::BoolOperation;

fn make_two_body_decomposition() -> CellDecomposition {
    // Simulate two overlapping cubes producing 3 cells:
    // Cell 0: in A only (membership = 0b01 = 1)
    // Cell 1: in A and B (membership = 0b11 = 3)
    // Cell 2: in B only (membership = 0b10 = 2)
    let cells = vec![
        Cell {
            id: 0,
            membership: 1,
            sample_point: Point3::new(0.25, 0.5, 0.5),
        },
        Cell {
            id: 1,
            membership: 3,
            sample_point: Point3::new(0.75, 0.5, 0.5),
        },
        Cell {
            id: 2,
            membership: 2,
            sample_point: Point3::new(1.25, 0.5, 0.5),
        },
    ];

    // 4 boundary faces:
    // Face 0: between cell 0 (A-only) and cell 1 (A∩B) — intersection boundary
    // Face 1: between cell 1 (A∩B) and cell 2 (B-only) — intersection boundary
    // Face 2: outer face of cell 0 (A-only) — between cell 0 and exterior (cell 3 implicit)
    // Face 3: outer face of cell 2 (B-only)
    // For simplicity, we model only the intersection boundaries + exterior as a large cell.
    let exterior_cell = 3; // not in cells vec — represents "outside both"
    let mut cells_with_exterior = cells;
    cells_with_exterior.push(Cell {
        id: 3,
        membership: 0, // outside both
        sample_point: Point3::new(5.0, 5.0, 5.0),
    });

    let boundaries = vec![
        CellBoundaryFace {
            fragment_index: 0,
            cell_plus: 0,
            cell_minus: 1,
        },
        CellBoundaryFace {
            fragment_index: 1,
            cell_plus: 1,
            cell_minus: 2,
        },
        CellBoundaryFace {
            fragment_index: 2,
            cell_plus: 3,
            cell_minus: 0,
        },
        CellBoundaryFace {
            fragment_index: 3,
            cell_plus: 2,
            cell_minus: 3,
        },
    ];

    CellDecomposition {
        cells: cells_with_exterior,
        boundaries,
    }
}

#[test]
fn cell_count() {
    let decomp = make_two_body_decomposition();
    assert_eq!(decomp.cell_count(), 4); // 3 interior + 1 exterior
}

#[test]
fn union_selects_all_body_cells() {
    let decomp = make_two_body_decomposition();
    let selected = decomp.select_union();
    assert_eq!(selected.len(), 3, "union should select cells 0, 1, 2");
    assert!(selected.contains(&0));
    assert!(selected.contains(&1));
    assert!(selected.contains(&2));
    assert!(!selected.contains(&3)); // exterior not selected
}

#[test]
fn intersection_selects_shared_cell() {
    let decomp = make_two_body_decomposition();
    let selected = decomp.select_intersection();
    assert_eq!(selected.len(), 1, "intersection should select cell 1 only");
    assert!(selected.contains(&1));
}

#[test]
fn subtraction_selects_a_minus_b() {
    let decomp = make_two_body_decomposition();
    let selected = decomp.select_subtraction();
    assert_eq!(selected.len(), 1, "subtraction should select cell 0 only");
    assert!(selected.contains(&0));
}

#[test]
fn symmetric_difference() {
    let decomp = make_two_body_decomposition();
    let selected = decomp.select_symmetric_difference();
    assert_eq!(selected.len(), 2, "xor should select cells 0 and 2");
    assert!(selected.contains(&0));
    assert!(selected.contains(&2));
}

#[test]
fn boundary_faces_of_union() {
    let decomp = make_two_body_decomposition();
    let selected = decomp.select_union();
    let faces = decomp.boundary_faces_of(&selected);
    // Union includes cells 0, 1, 2. Exterior is cell 3.
    // Faces between selected and exterior: faces 2, 3
    // Faces between two selected cells: faces 0, 1 (internal, excluded)
    assert_eq!(faces.len(), 2, "union boundary: faces 2 and 3");
}

#[test]
fn boundary_faces_of_intersection() {
    let decomp = make_two_body_decomposition();
    let selected = decomp.select_intersection();
    let faces = decomp.boundary_faces_of(&selected);
    // Intersection includes only cell 1.
    // Face 0: between cell 0 and 1 — cell 0 not selected, so this is a boundary
    // Face 1: between cell 1 and 2 — cell 2 not selected, so this is a boundary
    assert_eq!(faces.len(), 2, "intersection boundary: faces 0 and 1");
}

#[test]
fn select_for_operation_union() {
    let decomp = make_two_body_decomposition();
    let selected = select_for_operation(&decomp, BoolOperation::Union);
    assert_eq!(selected.len(), 3);
}

#[test]
fn custom_predicate() {
    let decomp = make_two_body_decomposition();
    // Select cells that are in NEITHER body (exterior)
    let selected = super::select::select_custom(&decomp, |m| m == 0);
    assert_eq!(selected.len(), 1);
    assert!(selected.contains(&3));
}

#[test]
fn cell_membership_queries() {
    let cell_a = Cell {
        id: 0,
        membership: 1,
        sample_point: Point3::new(0.0, 0.0, 0.0),
    };
    assert!(cell_a.in_a());
    assert!(!cell_a.in_b());

    let cell_ab = Cell {
        id: 1,
        membership: 3,
        sample_point: Point3::new(0.0, 0.0, 0.0),
    };
    assert!(cell_ab.in_a());
    assert!(cell_ab.in_b());
}
