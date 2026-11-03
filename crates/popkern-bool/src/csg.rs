//! CSG tree evaluator.

use popkern_core::{tolerance::Tolerance, topology::BRep};

use crate::error::BoolError;

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
    // One partition of every leaf, the tree read as a set expression
    // over them (design: operations on the cell complex, 5.1).
    let mut leaves = Vec::new();
    let expr = expression(tree, &mut leaves);
    if let [only] = leaves.as_slice() {
        return Ok(only.clone());
    }
    let refs: Vec<&BRep> = leaves.iter().collect();
    crate::partition::nary_boolean(&refs, &expr, *tolerance, &crate::BoolOptions::default())
}

/// `tree` as a set expression over its leaves, which it moves to
/// `leaves` in order.
fn expression(tree: CsgNode, leaves: &mut Vec<BRep>) -> crate::partition::SetExpr {
    use crate::partition::SetExpr;
    match tree {
        CsgNode::Leaf(brep) => {
            leaves.push(*brep);
            SetExpr::body(leaves.len() - 1)
        }
        CsgNode::Union(l, r) => expression(*l, leaves).or(expression(*r, leaves)),
        CsgNode::Intersection(l, r) => expression(*l, leaves).and(expression(*r, leaves)),
        CsgNode::Subtraction(l, r) => expression(*l, leaves).minus(expression(*r, leaves)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    #[test]
    fn leaf_returns_brep_unchanged() {
        let tol = Tolerance::default();
        let cube = BRep::unit_cube(tol).unwrap();
        let face_count = cube.faces.len();
        let result = csg_eval(CsgNode::Leaf(Box::new(cube)), &tol).unwrap();
        assert_eq!(result.faces.len(), face_count);
    }

    #[test]
    fn union_of_disjoint_cubes() {
        let tol = Tolerance::default();
        let a = BRep::unit_cube(tol).unwrap();
        // Moved whole (vertices and surfaces): moving only the
        // vertices left every face's plane behind.
        let b = BRep::unit_cube(tol)
            .unwrap()
            .transform(&popkern_core::geometry::Transform3::translation(
                5.0, 0.0, 0.0,
            ))
            .unwrap();
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
        let a = BRep::unit_cube(tol).unwrap();
        // Moved whole (vertices and surfaces): moving only the
        // vertices left every face's plane behind.
        let b = BRep::unit_cube(tol)
            .unwrap()
            .transform(&popkern_core::geometry::Transform3::translation(
                5.0, 0.0, 0.0,
            ))
            .unwrap();
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
        let a = BRep::unit_cube(tol).unwrap();
        let b = BRep::unit_cube(tol).unwrap();
        // Geometry-honest enlargement: vertex translation alone leaves the
        // stored surfaces stale (the classify hit-box gate rightly
        // rejects hits on a surface that no longer matches its face).
        let b = b
            .transform(&popkern_core::geometry::Transform3::scale_uniform(3.0))
            .unwrap()
            .transform(&popkern_core::geometry::Transform3::translation(
                -1.0, -1.0, -1.0,
            ))
            .unwrap();
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
        let a = BRep::unit_cube(tol).unwrap();
        // Moved whole (vertices and surfaces): moving only the
        // vertices left every face's plane behind.
        let b = BRep::unit_cube(tol)
            .unwrap()
            .transform(&popkern_core::geometry::Transform3::translation(
                5.0, 0.0, 0.0,
            ))
            .unwrap();
        let c = BRep::unit_cube(tol)
            .unwrap()
            .transform(&popkern_core::geometry::Transform3::translation(
                10.0, 0.0, 0.0,
            ))
            .unwrap();

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
