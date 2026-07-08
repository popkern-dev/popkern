//! Cell decomposition via flood-fill on fragment adjacency.
//!
//! Given a set of face fragments from the intersection graph, builds a
//! cellular decomposition by:
//! 1. Assigning each fragment to two cells (one on each side)
//! 2. Merging cells via union-find when fragments share edges
//! 3. Classifying each cell's membership via point-in-solid

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::classify::point_in_solid_final;
use crate::error::PointClassification;

use super::cell::{Cell, CellBoundaryFace, CellDecomposition, CellId};

/// Input parameters for cell decomposition.
pub struct DecomposeParams<'a> {
    /// Number of face fragments.
    pub fragment_count: usize,
    /// For each fragment, list of adjacent fragment indices.
    pub fragment_adjacency: &'a [Vec<usize>],
    /// Sample points for each fragment (used for cell classification).
    pub fragment_samples: &'a [Point3],
    /// Original body A.
    pub a: &'a BRep,
    /// Original body B.
    pub b: &'a BRep,
    /// Precomputed AABBs for body A's faces.
    pub aabbs_a: &'a [(FaceHandle, Aabb)],
    /// Precomputed AABBs for body B's faces.
    pub aabbs_b: &'a [(FaceHandle, Aabb)],
    /// Geometric tolerance.
    pub tolerance: &'a Tolerance,
}

/// Build a cellular decomposition from face fragments.
///
/// Each fragment separates two cells (positive-normal and negative-normal
/// sides). Adjacent fragments merge their corresponding cells via
/// union-find. Each resulting cell is classified by membership (which
/// input bodies contain it) via point-in-solid.
pub fn decompose_cells(params: &DecomposeParams<'_>) -> CellDecomposition {
    let n = params.fragment_count;
    let mut uf = UnionFind::new(n * 2);

    // For fragment i: cell_plus = 2*i, cell_minus = 2*i+1.
    // Adjacent fragments merge their same-side cells.
    for (i, neighbors) in params.fragment_adjacency.iter().enumerate() {
        for &j in neighbors {
            if j <= i {
                continue;
            }
            uf.union(2 * i, 2 * j);
            uf.union(2 * i + 1, 2 * j + 1);
        }
    }

    // Collect unique cells
    let mut cell_id_map: HashMap<usize, CellId> = HashMap::new();
    let mut cells: Vec<Cell> = Vec::new();

    for i in 0..(n * 2) {
        let root = uf.find(i);
        cell_id_map.entry(root).or_insert_with(|| {
            let id = cells.len();
            cells.push(Cell {
                id,
                membership: 0,
                sample_point: Point3::new(0.0, 0.0, 0.0),
            });
            id
        });
    }

    // Classify each cell by point-in-solid
    for cell in &mut cells {
        let mut sample = None;
        for (i, pt) in params.fragment_samples.iter().enumerate() {
            let plus_root = uf.find(2 * i);
            let minus_root = uf.find(2 * i + 1);
            if cell_id_map.get(&plus_root) == Some(&cell.id)
                || cell_id_map.get(&minus_root) == Some(&cell.id)
            {
                sample = Some(*pt);
                break;
            }
        }

        if let Some(pt) = sample {
            cell.sample_point = pt;
            if let Ok(PointClassification::Inside) =
                point_in_solid_final(&pt, params.a, params.aabbs_a, params.tolerance)
            {
                cell.membership |= 1;
            }
            if let Ok(PointClassification::Inside) =
                point_in_solid_final(&pt, params.b, params.aabbs_b, params.tolerance)
            {
                cell.membership |= 2;
            }
        }
    }

    // Build boundary face list
    let boundaries: Vec<_> = (0..n)
        .map(|i| {
            let cell_plus = cell_id_map[&uf.find(2 * i)];
            let cell_minus = cell_id_map[&uf.find(2 * i + 1)];
            CellBoundaryFace {
                fragment_index: i,
                cell_plus,
                cell_minus,
            }
        })
        .collect();

    CellDecomposition { cells, boundaries }
}

/// Simple union-find (disjoint set) data structure.
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        if self.rank[ra] < self.rank[rb] {
            self.parent[ra] = rb;
        } else if self.rank[ra] > self.rank[rb] {
            self.parent[rb] = ra;
        } else {
            self.parent[rb] = ra;
            self.rank[ra] += 1;
        }
    }
}
