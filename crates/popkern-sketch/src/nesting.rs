// Items are pub(crate) for use by the BRep builder module.
// Suppress dead_code until those consumers are added.
#![allow(dead_code)]

//! Nesting classification for traced loops.
//!
//! Determines which loops contain which other loops, building a tree
//! structure where even-depth nodes are faces and odd-depth nodes are
//! holes.

use crate::containment::point_in_loop;
use crate::loop_tracer::TracedLoop;
use popkern_constraint::Sketch;

/// A node in the nesting tree of traced loops.
#[derive(Debug, Clone)]
pub(crate) struct NestingNode {
    /// Index into the original `loops` slice.
    pub loop_index: usize,
    /// Child loops directly contained by this one.
    pub children: Vec<NestingNode>,
    /// `true` at even depth (0, 2, 4...) = face region.
    /// `false` at odd depth (1, 3, 5...) = hole region.
    pub is_face: bool,
}

/// Build the nesting tree from a set of traced loops.
///
/// 1. Compute an interior test point for each loop.
/// 2. Build containment relationships via point-in-loop tests.
/// 3. Determine direct containment (no intermediate container).
/// 4. Assign face/hole status by tree depth.
pub(crate) fn classify_nesting(loops: &[TracedLoop], sketch: &Sketch) -> Vec<NestingNode> {
    if loops.is_empty() {
        return Vec::new();
    }

    let n = loops.len();

    // Step 1: Compute test points for each loop.
    let test_points: Vec<(f64, f64)> = loops.iter().map(compute_test_point).collect();

    // Step 2: Build containment matrix.
    // contains[i][j] = true means loop i geometrically contains loop j.
    // A loop can only contain another if it has strictly larger area.
    let areas: Vec<f64> = loops.iter().map(|l| l.signed_area.abs()).collect();
    let mut contains = vec![vec![false; n]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            // A smaller-area loop cannot contain a larger-area one.
            if areas[i] <= areas[j] {
                continue;
            }
            contains[i][j] = point_in_loop(test_points[j], &loops[i], sketch);
        }
    }

    // Step 3: Determine direct containment.
    // i directly contains j if i contains j and there is no k such that
    // i contains k and k contains j.
    let mut direct_parent: Vec<Option<usize>> = vec![None; n];
    for j in 0..n {
        // Find all containers of j, then pick the most nested one
        // (the one that is contained by all others).
        let containers: Vec<usize> = (0..n).filter(|&i| contains[i][j]).collect();
        if containers.is_empty() {
            continue;
        }
        // The direct parent is the container that is contained by all
        // other containers of j (the "innermost" one).
        let direct = containers.iter().copied().find(|&c| {
            containers
                .iter()
                .all(|&other| other == c || contains[other][c])
        });
        direct_parent[j] = direct;
    }

    // Step 4: Build tree. Roots = loops with no parent.
    build_tree(&direct_parent, 0)
}

/// Build the nesting tree from the direct-parent array.
///
/// Returns the roots (loops not contained by any other).
fn build_tree(direct_parent: &[Option<usize>], root_depth: usize) -> Vec<NestingNode> {
    let n = direct_parent.len();

    // Collect children for each potential parent.
    let mut children_of: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut roots: Vec<usize> = Vec::new();

    for (j, parent) in direct_parent.iter().enumerate() {
        match parent {
            Some(p) => children_of[*p].push(j),
            None => roots.push(j),
        }
    }

    // Recursively build nodes.
    fn build_node(idx: usize, depth: usize, children_of: &[Vec<usize>]) -> NestingNode {
        let children: Vec<NestingNode> = children_of[idx]
            .iter()
            .map(|&child| build_node(child, depth + 1, children_of))
            .collect();
        NestingNode {
            loop_index: idx,
            children,
            is_face: depth.is_multiple_of(2),
        }
    }

    roots
        .iter()
        .map(|&r| build_node(r, root_depth, &children_of))
        .collect()
}

/// Compute an interior test point for a loop.
///
/// Uses the centroid of all edge start positions.
fn compute_test_point(loop_: &TracedLoop) -> (f64, f64) {
    if loop_.edges.is_empty() {
        return (0.0, 0.0);
    }
    let n = loop_.edges.len() as f64;
    let mut sx = 0.0;
    let mut sy = 0.0;
    for edge in &loop_.edges {
        sx += edge.start.0;
        sy += edge.start.1;
    }
    (sx / n, sy / n)
}

#[cfg(test)]
#[path = "nesting_tests.rs"]
mod tests;
