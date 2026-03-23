//! CSG tree evaluator.

use popkern_core::{tolerance::Tolerance, topology::BRep};

use crate::error::{BoolError, BoolOperation};
use crate::operation::execute_boolean;

/// A node in a CSG (Constructive Solid Geometry) expression tree.
#[derive(Debug)]
pub enum CsgNode {
    /// A leaf node containing a solid.
    Leaf(Box<BRep>),
    /// Union of two subtrees.
    Union(Box<CsgNode>, Box<CsgNode>),
    /// Intersection of two subtrees.
    Intersection(Box<CsgNode>, Box<CsgNode>),
    /// Subtraction: left minus right.
    Subtraction(Box<CsgNode>, Box<CsgNode>),
}

/// Evaluate a CSG expression tree, consuming all leaf BReps.
///
/// Recursively evaluates the tree bottom-up, applying boolean operations at
/// each interior node. Returns the resulting solid on success, or a
/// [`BoolError`] if any operation fails.
///
/// # Arguments
///
/// * `tree` — the CSG expression tree to evaluate.
/// * `tolerance` — geometric tolerance for all distance and coincidence tests.
///
/// # Errors
///
/// Returns a [`BoolError`] if any boolean operation in the tree fails.
pub fn csg_eval(tree: CsgNode, tolerance: &Tolerance) -> Result<BRep, BoolError> {
    match tree {
        CsgNode::Leaf(brep) => Ok(*brep),
        CsgNode::Union(left, right) => {
            let a = csg_eval(*left, tolerance)?;
            let b = csg_eval(*right, tolerance)?;
            execute_boolean(BoolOperation::Union, a, b, tolerance)
        }
        CsgNode::Intersection(left, right) => {
            let a = csg_eval(*left, tolerance)?;
            let b = csg_eval(*right, tolerance)?;
            execute_boolean(BoolOperation::Intersection, a, b, tolerance)
        }
        CsgNode::Subtraction(left, right) => {
            let a = csg_eval(*left, tolerance)?;
            let b = csg_eval(*right, tolerance)?;
            execute_boolean(BoolOperation::Subtraction, a, b, tolerance)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::{geometry::Vec3, tolerance::Tolerance};

    #[test]
    fn leaf_returns_brep_unchanged() {
        let tol = Tolerance::default();
        let cube = BRep::unit_cube(tol.clone()).unwrap();
        let face_count = cube.faces.len();
        let result = csg_eval(CsgNode::Leaf(Box::new(cube)), &tol).unwrap();
        assert_eq!(result.faces.len(), face_count);
    }

    #[test]
    fn union_of_disjoint_cubes() {
        let tol = Tolerance::default();
        let a = BRep::unit_cube(tol.clone()).unwrap();
        let mut b = BRep::unit_cube(tol.clone()).unwrap();
        for (_, v) in &mut b.vertices {
            v.position = v.position + Vec3::new(5.0, 0.0, 0.0);
        }
        let tree = CsgNode::Union(
            Box::new(CsgNode::Leaf(Box::new(a))),
            Box::new(CsgNode::Leaf(Box::new(b))),
        );
        let result = csg_eval(tree, &tol).unwrap();
        assert_eq!(result.faces.len(), 12); // 6 + 6
    }

    #[test]
    fn subtraction_of_disjoint_returns_left() {
        // A - B where B is far away = A
        let tol = Tolerance::default();
        let a = BRep::unit_cube(tol.clone()).unwrap();
        let mut b = BRep::unit_cube(tol.clone()).unwrap();
        for (_, v) in &mut b.vertices {
            v.position = v.position + Vec3::new(5.0, 0.0, 0.0);
        }
        let tree = CsgNode::Subtraction(
            Box::new(CsgNode::Leaf(Box::new(a))),
            Box::new(CsgNode::Leaf(Box::new(b))),
        );
        let result = csg_eval(tree, &tol).unwrap();
        assert_eq!(result.faces.len(), 6);
    }

    #[test]
    fn intersection_of_contained_cubes() {
        // A is a unit cube, B is a large cube containing A.
        // Intersection = A (6 faces).
        let tol = Tolerance::default();
        let a = BRep::unit_cube(tol.clone()).unwrap();
        let mut b = BRep::unit_cube(tol.clone()).unwrap();
        for (_, v) in &mut b.vertices {
            v.position = popkern_core::geometry::Point3::new(
                v.position.x * 3.0 - 1.0,
                v.position.y * 3.0 - 1.0,
                v.position.z * 3.0 - 1.0,
            );
        }
        let tree = CsgNode::Intersection(
            Box::new(CsgNode::Leaf(Box::new(a))),
            Box::new(CsgNode::Leaf(Box::new(b))),
        );
        let result = csg_eval(tree, &tol).unwrap();
        // A is inside B so intersection = A = 6 faces.
        assert_eq!(result.faces.len(), 6);
    }

    #[test]
    fn nested_tree_union_then_subtract() {
        // (A ∪ B) - C where all are disjoint
        let tol = Tolerance::default();
        let a = BRep::unit_cube(tol.clone()).unwrap();
        let mut b = BRep::unit_cube(tol.clone()).unwrap();
        for (_, v) in &mut b.vertices {
            v.position = v.position + Vec3::new(5.0, 0.0, 0.0);
        }
        let mut c = BRep::unit_cube(tol.clone()).unwrap();
        for (_, v) in &mut c.vertices {
            v.position = v.position + Vec3::new(10.0, 0.0, 0.0);
        }

        let tree = CsgNode::Subtraction(
            Box::new(CsgNode::Union(
                Box::new(CsgNode::Leaf(Box::new(a))),
                Box::new(CsgNode::Leaf(Box::new(b))),
            )),
            Box::new(CsgNode::Leaf(Box::new(c))),
        );
        let result = csg_eval(tree, &tol).unwrap();
        assert_eq!(result.faces.len(), 12); // (6+6) - 0 = 12
    }
}
