//! The design body's surface: the iso-surface of the design thresholded
//! at ½, closed to face connectivity.
//!
//! The design body is the part ∩ {x̄ ≥ ½}. Marching the densities
//! themselves puts each crossing between a solid and a void cell where
//! the two densities straddle ½: a crisp design's solid cells sit at
//! 0.97 here and 1.0 there, so a flat member face comes out as slivers
//! tilted by a degree or two (measured on the bored-block fixture: the
//! kept plate's underside split into dozens of 3–4 vertex faces at
//! z = 9.000 and 9.007). Exact features cutting across such slivers
//! decline. Thresholded, every solid-void crossing lies on the shared
//! cell face and a member's flat sides are exactly planar; the part's
//! own boundary still follows its cut cells' fill.
//!
//! Two solid cells that touch only along an edge or at a corner carry
//! load in the solver through their shared nodes (the one-node hinges of
//! density designs, Sigmund & Petersson 1998), but their surface meets
//! in a pinch with no section: measured on a 12 × 4 × 4 mm cantilever,
//! the body bent 3.45 mm where the design's compliance says 0.2. So the
//! thresholded field is closed first: wherever two solid cells touch
//! without a face path between them, the intermediate cells of highest
//! density are made solid, so every connection the solver used is at
//! least a cell face wide in the body. Contacts whose intermediate cells
//! lie outside the part cannot be closed and are counted.

use std::collections::{BTreeMap, BTreeSet};

use popkern_analysis::structural::generative::Design;
use popkern_analysis::structural::generative::surface::{SurfaceMesh, design_surface};

/// The body's surface and what the closing did.
#[derive(Debug, Clone, PartialEq)]
pub struct BodySurface {
    /// The surface, mm.
    pub mesh: SurfaceMesh,
    /// Cells made solid to close edge and corner contacts.
    pub closed_cells: usize,
    /// Edge or corner contacts left open (their intermediate cells are
    /// outside the part).
    pub open_contacts: usize,
}

type Cell = [i64; 3];

fn add(a: Cell, b: Cell) -> Cell {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// The surface of `design`'s body: its densities thresholded at ½ (each
/// entry solid or void), closed to face connectivity, then the
/// iso-surface at ½.
pub fn body_surface(design: &Design) -> BodySurface {
    let h = design.cell_size;
    let origin = design.entries.first().map_or([0.0; 3], |e| e.corner);
    let key = |c: [f64; 3]| -> Cell { [0, 1, 2].map(|a| ((c[a] - origin[a]) / h).round() as i64) };
    // Each grid cell's entries (a split cell has several) and its
    // densest fragment.
    let mut cells: BTreeMap<Cell, (Vec<usize>, f64)> = BTreeMap::new();
    for (i, (e, &d)) in design.entries.iter().zip(&design.densities).enumerate() {
        let slot = cells.entry(key(e.corner)).or_insert((Vec::new(), 0.0));
        slot.0.push(i);
        slot.1 = slot.1.max(d);
    }
    let mut solid: BTreeSet<Cell> = cells
        .iter()
        .filter(|(_, (_, d))| *d >= 0.5)
        .map(|(c, _)| *c)
        .collect();
    let density = |c: &Cell| cells.get(c).map(|(_, d)| *d);
    let mut closed = Vec::new();
    let mut open;
    // Each contact once per pass, from its lower cell, in key order; a
    // closed cell can make new contacts, so passes repeat until none
    // closes (the solid set only grows: it ends).
    loop {
        let before = closed.len();
        open = 0;
        let order: Vec<Cell> = solid.iter().copied().collect();
        for c in order {
            for dz in -1..=1 {
                for dy in -1..=1 {
                    for dx in -1..=1 {
                        let off = [dx, dy, dz];
                        let axes: Vec<usize> = (0..3).filter(|&a| off[a] != 0).collect();
                        let n = add(c, off);
                        if axes.len() < 2 || n <= c || !solid.contains(&n) {
                            continue;
                        }
                        // Face paths from c to n: one step per axis, in
                        // every order; their intermediate cells.
                        let unit = |a: usize| {
                            let mut u = [0; 3];
                            u[a] = off[a];
                            u
                        };
                        let orders: Vec<Vec<usize>> = if axes.len() == 2 {
                            vec![vec![axes[0], axes[1]], vec![axes[1], axes[0]]]
                        } else {
                            vec![
                                vec![0, 1, 2],
                                vec![0, 2, 1],
                                vec![1, 0, 2],
                                vec![1, 2, 0],
                                vec![2, 0, 1],
                                vec![2, 1, 0],
                            ]
                        };
                        let paths: Vec<Vec<Cell>> = orders
                            .iter()
                            .map(|o| {
                                let mut at = c;
                                o[..o.len() - 1]
                                    .iter()
                                    .map(|&a| {
                                        at = add(at, unit(a));
                                        at
                                    })
                                    .collect()
                            })
                            .collect();
                        if paths.iter().any(|p| p.iter().all(|q| solid.contains(q))) {
                            continue;
                        }
                        // The path inside the part with the densest weakest
                        // cell (ties: the first in order).
                        let best = paths
                            .iter()
                            .filter_map(|p| {
                                let weakest = p
                                    .iter()
                                    .map(&density)
                                    .try_fold(f64::INFINITY, |m, d| d.map(|d| m.min(d)))?;
                                Some((weakest, p))
                            })
                            .fold(None::<(f64, &Vec<Cell>)>, |best, (w, p)| match best {
                                Some((bw, _)) if bw >= w => best,
                                _ => Some((w, p)),
                            });
                        match best {
                            Some((_, path)) => {
                                for q in path {
                                    if solid.insert(*q) {
                                        closed.push(*q);
                                    }
                                }
                            }
                            None => open += 1,
                        }
                    }
                }
            }
        }
        if closed.len() == before {
            break;
        }
    }
    let mut densities: Vec<f64> = design
        .densities
        .iter()
        .map(|&d| if d >= 0.5 { 1.0 } else { 0.0 })
        .collect();
    for q in &closed {
        if let Some((entries, _)) = cells.get(q) {
            for &i in entries {
                densities[i] = 1.0;
            }
        }
    }
    let crisp = Design {
        densities,
        ..design.clone()
    };
    BodySurface {
        mesh: design_surface(&crisp, 0.5),
        closed_cells: closed.len(),
        open_contacts: open,
    }
}

#[cfg(test)]
#[path = "surface_tests.rs"]
mod tests;
