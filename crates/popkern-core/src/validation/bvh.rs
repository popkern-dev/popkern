//! Bounding Volume Hierarchy for spatial queries.

use crate::geometry::Point3;
use crate::geometry::intersection::surface_subdivision::Aabb;

/// A Bounding Volume Hierarchy built from a set of (id, AABB) items.
///
/// Supports efficient overlap queries between all pairs of items.
pub struct Bvh {
    nodes: Vec<BvhNode>,
    /// The original items stored by the BVH (id, bbox).
    items: Vec<(usize, Aabb)>,
}

struct BvhNode {
    bbox: Aabb,
    data: BvhNodeData,
}

enum BvhNodeData {
    Leaf(Vec<usize>), // indices into self.items
    Internal { left: usize, right: usize },
}

/// Compute the combined AABB of a set of items (by their indices into `items`).
fn combined_aabb(indices: &[usize], items: &[(usize, Aabb)]) -> Aabb {
    let mut min = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &idx in indices {
        let b = &items[idx].1;
        min = Point3::new(min.x.min(b.min.x), min.y.min(b.min.y), min.z.min(b.min.z));
        max = Point3::new(max.x.max(b.max.x), max.y.max(b.max.y), max.z.max(b.max.z));
    }
    Aabb { min, max }
}

/// Compute the centroid of an AABB along the specified axis (0=x, 1=y, 2=z).
fn centroid_along(aabb: &Aabb, axis: usize) -> f64 {
    match axis {
        0 => (aabb.min.x + aabb.max.x) * 0.5,
        1 => (aabb.min.y + aabb.max.y) * 0.5,
        _ => (aabb.min.z + aabb.max.z) * 0.5,
    }
}

impl Bvh {
    /// Build a BVH from a slice of `(item_id, bounding_box)` pairs.
    ///
    /// Uses recursive median splitting along the longest axis.
    /// Nodes with ≤ 4 items become leaves.
    pub fn build(items: &[(usize, Aabb)]) -> Self {
        if items.is_empty() {
            return Self {
                nodes: Vec::new(),
                items: Vec::new(),
            };
        }

        let items_vec: Vec<(usize, Aabb)> = items.to_vec();
        // Start with indices for all items.
        let all_indices: Vec<usize> = (0..items_vec.len()).collect();
        let mut nodes: Vec<BvhNode> = Vec::new();
        build_recursive(&all_indices, &items_vec, &mut nodes);

        Self {
            nodes,
            items: items_vec,
        }
    }

    /// Find all pairs of items whose bounding boxes overlap (with optional padding).
    ///
    /// Returns a `Vec` of `(id_a, id_b)` with `id_a < id_b`.
    pub fn find_overlapping_pairs(&self, padding: f64) -> Vec<(usize, usize)> {
        if self.nodes.is_empty() {
            return Vec::new();
        }
        let mut pairs = Vec::new();
        // Traverse the tree: check subtree pairs recursively.
        traverse_pairs(0, 0, true, &self.nodes, &self.items, padding, &mut pairs);
        pairs.sort_unstable();
        pairs.dedup();
        pairs
    }
}

/// Recursively build the BVH. Returns the index of the newly created node.
fn build_recursive(indices: &[usize], items: &[(usize, Aabb)], nodes: &mut Vec<BvhNode>) -> usize {
    let bbox = combined_aabb(indices, items);
    let node_idx = nodes.len();

    if indices.len() <= 4 {
        nodes.push(BvhNode {
            bbox,
            data: BvhNodeData::Leaf(indices.to_vec()),
        });
        return node_idx;
    }

    // Find longest axis.
    let dx = bbox.max.x - bbox.min.x;
    let dy = bbox.max.y - bbox.min.y;
    let dz = bbox.max.z - bbox.min.z;
    let axis = if dx >= dy && dx >= dz {
        0
    } else if dy >= dz {
        1
    } else {
        2
    };

    // Sort indices by centroid along the chosen axis.
    let mut sorted = indices.to_vec();
    sorted.sort_unstable_by(|&a, &b| {
        centroid_along(&items[a].1, axis)
            .partial_cmp(&centroid_along(&items[b].1, axis))
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mid = sorted.len() / 2;
    let left_indices = &sorted[..mid];
    let right_indices = &sorted[mid..];

    // Reserve a slot for this internal node first, then build children.
    nodes.push(BvhNode {
        bbox,
        data: BvhNodeData::Leaf(Vec::new()), // placeholder
    });

    let left_idx = build_recursive(left_indices, items, nodes);
    let right_idx = build_recursive(right_indices, items, nodes);

    // Patch the placeholder with the real Internal node data.
    nodes[node_idx].data = BvhNodeData::Internal {
        left: left_idx,
        right: right_idx,
    };

    node_idx
}

/// Traverse the BVH collecting overlapping item pairs.
///
/// `same_node` indicates whether `node_a` and `node_b` are the same node
/// (used to avoid re-testing a subtree against itself when it has no overlap).
fn traverse_pairs(
    node_a: usize,
    node_b: usize,
    same_node: bool,
    nodes: &[BvhNode],
    items: &[(usize, Aabb)],
    padding: f64,
    pairs: &mut Vec<(usize, usize)>,
) {
    let a = &nodes[node_a];
    let b = &nodes[node_b];

    // Early-out: if testing two different subtrees and their bboxes don't overlap, skip.
    if !same_node && !a.bbox.overlaps(&b.bbox, padding) {
        return;
    }

    match (&a.data, &b.data) {
        (BvhNodeData::Leaf(la), BvhNodeData::Leaf(lb)) => {
            // Test all pairs between the two (or within one) leaf.
            if same_node {
                // Within the same leaf: test all unique pairs.
                for i in 0..la.len() {
                    for j in (i + 1)..la.len() {
                        let ia = la[i];
                        let ib = la[j];
                        if items[ia].1.overlaps(&items[ib].1, padding) {
                            let id_a = items[ia].0;
                            let id_b = items[ib].0;
                            let pair = if id_a < id_b {
                                (id_a, id_b)
                            } else {
                                (id_b, id_a)
                            };
                            pairs.push(pair);
                        }
                    }
                }
            } else {
                // Between two different leaves.
                for &ia in la.iter() {
                    for &ib in lb.iter() {
                        if items[ia].1.overlaps(&items[ib].1, padding) {
                            let id_a = items[ia].0;
                            let id_b = items[ib].0;
                            if id_a != id_b {
                                let pair = if id_a < id_b {
                                    (id_a, id_b)
                                } else {
                                    (id_b, id_a)
                                };
                                pairs.push(pair);
                            }
                        }
                    }
                }
            }
        }
        (
            BvhNodeData::Internal {
                left: la,
                right: ra,
            },
            BvhNodeData::Leaf(_),
        ) => {
            // Descend into the internal node's children against the leaf.
            traverse_pairs(*la, node_b, false, nodes, items, padding, pairs);
            traverse_pairs(*ra, node_b, false, nodes, items, padding, pairs);
        }
        (
            BvhNodeData::Leaf(_),
            BvhNodeData::Internal {
                left: lb,
                right: rb,
            },
        ) => {
            traverse_pairs(node_a, *lb, false, nodes, items, padding, pairs);
            traverse_pairs(node_a, *rb, false, nodes, items, padding, pairs);
        }
        (
            BvhNodeData::Internal {
                left: la,
                right: ra,
            },
            BvhNodeData::Internal {
                left: lb,
                right: rb,
            },
        ) => {
            if same_node {
                // Same internal node: recurse into each child with itself, then cross-test.
                let (la, ra, lb, rb) = (*la, *ra, *lb, *rb);
                traverse_pairs(la, la, true, nodes, items, padding, pairs);
                traverse_pairs(ra, ra, true, nodes, items, padding, pairs);
                traverse_pairs(la, rb, false, nodes, items, padding, pairs);
                traverse_pairs(ra, lb, false, nodes, items, padding, pairs);
            } else {
                let (la, ra, lb, rb) = (*la, *ra, *lb, *rb);
                traverse_pairs(la, lb, false, nodes, items, padding, pairs);
                traverse_pairs(la, rb, false, nodes, items, padding, pairs);
                traverse_pairs(ra, lb, false, nodes, items, padding, pairs);
                traverse_pairs(ra, rb, false, nodes, items, padding, pairs);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bvh_finds_overlapping_pair() {
        let items = vec![
            (
                0,
                Aabb {
                    min: Point3::new(0.0, 0.0, 0.0),
                    max: Point3::new(1.0, 1.0, 1.0),
                },
            ),
            (
                1,
                Aabb {
                    min: Point3::new(0.5, 0.5, 0.5),
                    max: Point3::new(1.5, 1.5, 1.5),
                },
            ),
            (
                2,
                Aabb {
                    min: Point3::new(5.0, 5.0, 5.0),
                    max: Point3::new(6.0, 6.0, 6.0),
                },
            ),
        ];
        let bvh = Bvh::build(&items);
        let pairs = bvh.find_overlapping_pairs(0.0);
        assert!(pairs.contains(&(0, 1)));
        assert!(!pairs.iter().any(|p| p.0 == 2 || p.1 == 2));
    }

    #[test]
    fn bvh_empty_input() {
        let bvh = Bvh::build(&[]);
        assert!(bvh.find_overlapping_pairs(0.0).is_empty());
    }

    #[test]
    fn bvh_single_item() {
        let items = vec![(
            0,
            Aabb {
                min: Point3::new(0.0, 0.0, 0.0),
                max: Point3::new(1.0, 1.0, 1.0),
            },
        )];
        let bvh = Bvh::build(&items);
        assert!(bvh.find_overlapping_pairs(0.0).is_empty());
    }

    #[test]
    fn bvh_two_non_overlapping_items() {
        let items = vec![
            (
                0,
                Aabb {
                    min: Point3::new(0.0, 0.0, 0.0),
                    max: Point3::new(1.0, 1.0, 1.0),
                },
            ),
            (
                1,
                Aabb {
                    min: Point3::new(2.0, 2.0, 2.0),
                    max: Point3::new(3.0, 3.0, 3.0),
                },
            ),
        ];
        let bvh = Bvh::build(&items);
        assert!(bvh.find_overlapping_pairs(0.0).is_empty());
    }

    #[test]
    fn bvh_two_overlapping_items() {
        let items = vec![
            (
                0,
                Aabb {
                    min: Point3::new(0.0, 0.0, 0.0),
                    max: Point3::new(1.0, 1.0, 1.0),
                },
            ),
            (
                1,
                Aabb {
                    min: Point3::new(0.5, 0.5, 0.5),
                    max: Point3::new(1.5, 1.5, 1.5),
                },
            ),
        ];
        let bvh = Bvh::build(&items);
        let pairs = bvh.find_overlapping_pairs(0.0);
        assert_eq!(pairs, vec![(0, 1)]);
    }

    #[test]
    fn bvh_padding_bridges_gap() {
        let items = vec![
            (
                0,
                Aabb {
                    min: Point3::new(0.0, 0.0, 0.0),
                    max: Point3::new(1.0, 1.0, 1.0),
                },
            ),
            (
                1,
                Aabb {
                    min: Point3::new(1.5, 0.0, 0.0),
                    max: Point3::new(2.5, 1.0, 1.0),
                },
            ),
        ];
        let bvh = Bvh::build(&items);
        // Without padding: no overlap (gap of 0.5).
        assert!(bvh.find_overlapping_pairs(0.0).is_empty());
        // With padding of 0.6: should overlap.
        let pairs = bvh.find_overlapping_pairs(0.6);
        assert!(pairs.contains(&(0, 1)));
    }

    #[test]
    fn bvh_many_items_triggers_split() {
        // Build 10 items so we force an internal node split (> 4 threshold).
        let items: Vec<(usize, Aabb)> = (0..10)
            .map(|i| {
                let x = i as f64 * 2.0;
                (
                    i,
                    Aabb {
                        min: Point3::new(x, 0.0, 0.0),
                        max: Point3::new(x + 1.0, 1.0, 1.0),
                    },
                )
            })
            .collect();
        let bvh = Bvh::build(&items);
        // No items overlap (gap of 1.0 between each).
        assert!(bvh.find_overlapping_pairs(0.0).is_empty());
    }

    #[test]
    fn bvh_many_items_some_overlapping() {
        // Items at 0..1 and 0.5..1.5 overlap; item at 10..11 does not.
        let items = vec![
            (
                0,
                Aabb {
                    min: Point3::new(0.0, 0.0, 0.0),
                    max: Point3::new(1.0, 1.0, 1.0),
                },
            ),
            (
                1,
                Aabb {
                    min: Point3::new(0.5, 0.0, 0.0),
                    max: Point3::new(1.5, 1.0, 1.0),
                },
            ),
            (
                2,
                Aabb {
                    min: Point3::new(10.0, 0.0, 0.0),
                    max: Point3::new(11.0, 1.0, 1.0),
                },
            ),
            (
                3,
                Aabb {
                    min: Point3::new(20.0, 0.0, 0.0),
                    max: Point3::new(21.0, 1.0, 1.0),
                },
            ),
            (
                4,
                Aabb {
                    min: Point3::new(30.0, 0.0, 0.0),
                    max: Point3::new(31.0, 1.0, 1.0),
                },
            ),
        ];
        let bvh = Bvh::build(&items);
        let pairs = bvh.find_overlapping_pairs(0.0);
        assert!(pairs.contains(&(0, 1)));
        // items 2, 3, 4 should not overlap with anything
        assert!(!pairs.iter().any(|p| p.0 == 2 || p.1 == 2));
        assert!(!pairs.iter().any(|p| p.0 == 3 || p.1 == 3));
        assert!(!pairs.iter().any(|p| p.0 == 4 || p.1 == 4));
    }

    #[test]
    fn bvh_result_pairs_are_ordered() {
        // All returned pairs should have first element < second.
        let items = vec![
            (
                5,
                Aabb {
                    min: Point3::new(0.0, 0.0, 0.0),
                    max: Point3::new(1.0, 1.0, 1.0),
                },
            ),
            (
                3,
                Aabb {
                    min: Point3::new(0.5, 0.5, 0.5),
                    max: Point3::new(1.5, 1.5, 1.5),
                },
            ),
        ];
        let bvh = Bvh::build(&items);
        let pairs = bvh.find_overlapping_pairs(0.0);
        for (a, b) in &pairs {
            assert!(a < b, "Pair ({a}, {b}) not in canonical order");
        }
    }
}
