//! Graph decomposition of the constraint system into independent sub-problems.
//!
//! Uses a union-find (disjoint set) data structure to partition the sketch's
//! parameter vector into connected components, where connectivity is defined
//! by shared constraint references. Independent components can be solved
//! separately, reducing the size of each Newton-Raphson system.

use crate::constraint::{ConstraintHandle, ConstraintKind};
use crate::elements::ElementHandle;
use crate::rigid_set;
use crate::Sketch;
use std::collections::HashMap;

/// Disjoint-set (union-find) data structure with path compression and
/// union-by-rank for near-constant amortised operations.
pub(crate) struct UnionFind {
    /// Parent pointers; `parent[x] == x` means `x` is a root.
    parent: Vec<usize>,
    /// Tree rank for union-by-rank heuristic.
    rank: Vec<usize>,
}

impl UnionFind {
    /// Create a new union-find with `n` singleton sets.
    pub(crate) fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    /// Find the representative (root) of the set containing `x`, with
    /// path compression.
    pub(crate) fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    /// Merge the sets containing `x` and `y` using union-by-rank.
    pub(crate) fn union(&mut self, x: usize, y: usize) {
        let rx = self.find(x);
        let ry = self.find(y);
        if rx == ry {
            return;
        }
        match self.rank[rx].cmp(&self.rank[ry]) {
            std::cmp::Ordering::Less => self.parent[rx] = ry,
            std::cmp::Ordering::Greater => self.parent[ry] = rx,
            std::cmp::Ordering::Equal => {
                self.parent[ry] = rx;
                self.rank[rx] += 1;
            }
        }
    }
}

/// A connected sub-problem: a subset of parameters and the constraints
/// that reference them.
pub(crate) struct SubProblem {
    /// Parameter indices in the global `sketch.params` vector, sorted.
    pub global_indices: Vec<usize>,
    /// Map from global parameter index to local (component) index.
    pub local_indices: HashMap<usize, usize>,
    /// Constraint handles belonging to this component.
    pub constraints: Vec<ConstraintHandle>,
}

/// Collect all parameter indices owned by an element.
///
/// For points: 2 indices. For arcs/ellipses/hyperbolas: 5. For parabolas: 4.
/// Lines have no own params (they reference point params). Splines reference
/// their control-point params.
fn element_param_indices(sketch: &Sketch, elem: &ElementHandle) -> Vec<usize> {
    match elem {
        ElementHandle::Point(ph) => {
            if let Some(pt) = sketch.points.get(*ph) {
                let [ix, iy] = pt.param_indices();
                vec![ix, iy]
            } else {
                Vec::new()
            }
        }
        ElementHandle::Line(lh) => {
            if let Some(line) = sketch.lines.get(*lh) {
                let mut indices = Vec::new();
                if let Some(pt) = sketch.points.get(line.start) {
                    let [ix, iy] = pt.param_indices();
                    indices.push(ix);
                    indices.push(iy);
                }
                if let Some(pt) = sketch.points.get(line.end) {
                    let [ix, iy] = pt.param_indices();
                    indices.push(ix);
                    indices.push(iy);
                }
                indices
            } else {
                Vec::new()
            }
        }
        ElementHandle::Arc(ah) => {
            if let Some(arc) = sketch.arcs.get(*ah) {
                arc.param_indices().to_vec()
            } else {
                Vec::new()
            }
        }
        ElementHandle::Ellipse(eh) => {
            if let Some(ellipse) = sketch.ellipses.get(*eh) {
                ellipse.param_indices().to_vec()
            } else {
                Vec::new()
            }
        }
        ElementHandle::Parabola(ph) => {
            if let Some(parabola) = sketch.parabolas.get(*ph) {
                parabola.param_indices().to_vec()
            } else {
                Vec::new()
            }
        }
        ElementHandle::Hyperbola(hh) => {
            if let Some(hyperbola) = sketch.hyperbolas.get(*hh) {
                hyperbola.param_indices().to_vec()
            } else {
                Vec::new()
            }
        }
        ElementHandle::Spline(sh) => {
            if let Some(spline) = sketch.splines.get(*sh) {
                let mut indices = Vec::new();
                for &cp in &spline.control_points {
                    if let Some(pt) = sketch.points.get(cp) {
                        let [ix, iy] = pt.param_indices();
                        indices.push(ix);
                        indices.push(iy);
                    }
                }
                indices
            } else {
                Vec::new()
            }
        }
    }
}

/// Collect all parameter indices touched by a constraint (element params +
/// auxiliary params).
fn constraint_param_indices(
    sketch: &Sketch,
    handle: ConstraintHandle,
    entry: &crate::constraint::ConstraintEntry,
) -> Vec<usize> {
    let elements = rigid_set::constraint_element_handles(&entry.constraint);
    let mut indices: Vec<usize> = elements
        .iter()
        .flat_map(|elem| element_param_indices(sketch, elem))
        .collect();
    // Include auxiliary params (e.g. spline contact parameter t).
    if let Some(&aux_start) = sketch.aux_param_map.get(&handle) {
        let count = entry.constraint.aux_param_count();
        for i in 0..count {
            indices.push(aux_start + i);
        }
    }
    indices
}

/// Decompose the sketch's constraint graph into connected components.
///
/// Each component is an independent sub-problem that can be solved
/// separately. Driven constraints, constraints in disabled groups, and
/// intra-rigid-set constraints are excluded.
///
/// Returns an empty vector when there are no driving constraints.
pub(crate) fn find_components(sketch: &Sketch) -> Vec<SubProblem> {
    let n = sketch.params.len();
    if n == 0 {
        return Vec::new();
    }

    let mut uf = UnionFind::new(n);

    // Phase 1: union params within each driving constraint.
    let mut active_constraints: Vec<(ConstraintHandle, Vec<usize>)> = Vec::new();

    for (handle, entry) in sketch.constraints.iter() {
        if entry.kind != ConstraintKind::Driving {
            continue;
        }
        if let Some(gh) = entry.group {
            if let Some(g) = sketch.groups.get(gh) {
                if !g.enabled {
                    continue;
                }
            }
        }
        if rigid_set::all_in_same_rigid_set(sketch, &entry.constraint).is_some() {
            continue;
        }

        let indices = constraint_param_indices(sketch, handle, entry);
        if indices.is_empty() {
            continue;
        }

        // Union all params in this constraint together.
        let first = indices[0];
        for &idx in &indices[1..] {
            uf.union(first, idx);
        }

        active_constraints.push((handle, indices));
    }

    // Phase 2: union rigid-set body params with their member point params.
    for rs in sketch.rigid_sets.values() {
        let rs_tx = rs.param_start;
        let rs_ty = rs.param_start + 1;
        let rs_theta = rs.param_start + 2;

        // Union the three rigid-body params together.
        uf.union(rs_tx, rs_ty);
        uf.union(rs_tx, rs_theta);

        // Union each member point's params with the rigid-body params.
        for &ph in &rs.points {
            if let Some(pt) = sketch.points.get(ph) {
                let [ix, iy] = pt.param_indices();
                uf.union(rs_tx, ix);
                uf.union(rs_tx, iy);
            }
        }

        // Union centre params of arc/conic elements.
        for (elem, _) in &rs.center_offsets {
            let ps = match elem {
                ElementHandle::Arc(ah) => sketch.arcs.get(*ah).map(|a| a.param_start),
                ElementHandle::Ellipse(eh) => sketch.ellipses.get(*eh).map(|e| e.param_start),
                ElementHandle::Parabola(ph) => sketch.parabolas.get(*ph).map(|p| p.param_start),
                ElementHandle::Hyperbola(hh) => sketch.hyperbolas.get(*hh).map(|h| h.param_start),
                _ => None,
            };
            if let Some(ps) = ps {
                uf.union(rs_tx, ps);
                uf.union(rs_tx, ps + 1);
            }
        }
    }

    // Phase 3: group params by root to form components.
    let mut root_to_params: HashMap<usize, Vec<usize>> = HashMap::new();
    // Only include params that are referenced by at least one active constraint.
    let mut referenced = vec![false; n];
    for (_, indices) in &active_constraints {
        for &idx in indices {
            referenced[idx] = true;
        }
    }
    // Also mark rigid-set params as referenced if any of their points are.
    for rs in sketch.rigid_sets.values() {
        let rs_tx = rs.param_start;
        let any_referenced = rs.points.iter().any(|&ph| {
            sketch.points.get(ph).is_some_and(|pt| {
                let [ix, _] = pt.param_indices();
                referenced[ix]
            })
        });
        if any_referenced {
            referenced[rs_tx] = true;
            referenced[rs_tx + 1] = true;
            referenced[rs_tx + 2] = true;
        }
    }
    for (i, &is_ref) in referenced.iter().enumerate() {
        if is_ref {
            let root = uf.find(i);
            root_to_params.entry(root).or_default().push(i);
        }
    }

    if root_to_params.is_empty() {
        return Vec::new();
    }

    // Phase 4: assign constraints to components.
    let mut root_to_constraints: HashMap<usize, Vec<ConstraintHandle>> = HashMap::new();
    for (handle, indices) in &active_constraints {
        if let Some(&first_idx) = indices.first() {
            let root = uf.find(first_idx);
            root_to_constraints.entry(root).or_default().push(*handle);
        }
    }

    // Phase 5: build SubProblem structs.
    let mut components: Vec<SubProblem> = Vec::new();
    for (root, mut params) in root_to_params {
        params.sort_unstable();
        params.dedup();
        let local_indices: HashMap<usize, usize> = params
            .iter()
            .enumerate()
            .map(|(local, &global)| (global, local))
            .collect();
        let constraints = root_to_constraints.remove(&root).unwrap_or_default();
        components.push(SubProblem {
            global_indices: params,
            local_indices,
            constraints,
        });
    }

    components
}

#[cfg(test)]
#[path = "decompose_tests.rs"]
mod tests;
