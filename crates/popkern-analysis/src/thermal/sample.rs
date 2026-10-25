//! Field values and gradients at mesh vertices by trilinear
//! interpolation of nodal values in an active cell containing the
//! vertex.

use crate::aggregate::Extension;
use crate::element::poly::{eval, shape, shape_derivative};
use crate::error::AnalysisError;
use crate::grid::Grid;
use crate::model::SolidMesh;
use crate::moments::SolidMoments;

/// Relative slack (of the cell size) within which a point on or just
/// past a cell plane counts as inside the cell on either side. Surface
/// vertices can sit an ulp outside every active cell holding them.
pub const PLANE_SLACK: f64 = 1.0e-9;

/// Candidate cells along one axis: the cell containing `x`, then any
/// neighbour whose closed extent holds `x` within [`PLANE_SLACK`].
fn axis_candidates(grid: &Grid, axis: usize, x: f64) -> Vec<usize> {
    let t = (x - grid.origin()[axis]) / grid.cell_size();
    let n = grid.dims()[axis] as f64;
    let mut out: Vec<usize> = Vec::new();
    for v in [t, t - PLANE_SLACK, t + PLANE_SLACK] {
        let f = v.floor();
        if f >= 0.0 && f < n && !out.contains(&(f as usize)) {
            out.push(f as usize);
        }
    }
    out
}

/// The value of nodal field `values` (parallel to `extension.nodes`) at
/// point `p`, from an active cell of `moments` containing it; `None`
/// when no active cell contains it.
pub fn point_value(
    grid: &Grid,
    moments: &SolidMoments,
    extension: &Extension,
    values: &[f64],
    p: [f64; 3],
) -> Result<Option<f64>, AnalysisError> {
    Ok(point_eval(grid, moments, extension, values, p)?.map(|(v, _)| v))
}

/// The value and gradient (per mm) of nodal field `values` at `p`, from
/// an active cell containing it; `None` when no active cell does.
pub fn point_eval(
    grid: &Grid,
    moments: &SolidMoments,
    extension: &Extension,
    values: &[f64],
    p: [f64; 3],
) -> Result<Option<(f64, [f64; 3])>, AnalysisError> {
    let h = grid.cell_size();
    let [cx, cy, cz] = [0, 1, 2].map(|a| axis_candidates(grid, a, p[a]));
    for &k in &cz {
        for &j in &cy {
            for &i in &cx {
                let Some(cell) = grid.cell_index([i, j, k]) else {
                    continue;
                };
                let Some(entry) = moments.entry_at(grid, cell, p, None) else {
                    continue;
                };
                let lo = grid.cell_min([i, j, k]);
                let xi = [0, 1, 2].map(|a| 2.0 * (p[a] - lo[a]) / h - 1.0);
                let nodes = moments.node_keys(grid, entry).unwrap_or_default();
                let mut sum = 0.0;
                let mut gradient = [0.0; 3];
                for (a, node) in nodes.iter().enumerate() {
                    let q = extension
                        .position(*node)
                        .ok_or_else(|| AnalysisError::Internal {
                            stage: "sample",
                            detail: format!("node {node} of an active cell has no value"),
                        })?;
                    sum += eval(&shape(a), xi) * values[q];
                    for (k, g) in gradient.iter_mut().enumerate() {
                        *g += 2.0 / h * eval(&shape_derivative(a, k), xi) * values[q];
                    }
                }
                return Ok(Some((sum, gradient)));
            }
        }
    }
    Ok(None)
}

/// The value and gradient (per mm) of nodal field `values` at every
/// vertex of `mesh`.
pub fn vertex_values(
    grid: &Grid,
    mesh: &SolidMesh,
    moments: &SolidMoments,
    extension: &Extension,
    values: &[f64],
) -> Result<Vec<(f64, [f64; 3])>, AnalysisError> {
    mesh.vertices
        .iter()
        .enumerate()
        .map(|(vi, v)| {
            point_eval(grid, moments, extension, values, *v)?.ok_or_else(|| {
                AnalysisError::Internal {
                    stage: "sample",
                    detail: format!("vertex {vi} lies in no active cell"),
                }
            })
        })
        .collect()
}

/// The temperature gradient (per mm) at every vertex of `mesh`,
/// recovered like stresses (map §1.6): each entry's Q1 gradient at its
/// fragment centroid feeds the nodal SPR fit, and a vertex interpolates
/// its cell's nodal values; where a node has no value the vertex takes
/// its cell's own gradient.
pub fn recovered_gradients(
    grid: &Grid,
    mesh: &SolidMesh,
    moments: &SolidMoments,
    extension: &Extension,
    values: &[f64],
) -> Result<Vec<[f64; 3]>, AnalysisError> {
    use crate::recover::cell::centroid;
    use crate::recover::spr::Sample;
    let h = grid.cell_size();
    let missing = |node: usize| AnalysisError::Internal {
        stage: "sample",
        detail: format!("node {node} of an active cell has no value"),
    };
    // The entry's nodal values and its gradient at local point xi.
    let cell_gradient = |entry: usize, xi: [f64; 3]| -> Result<[f64; 3], AnalysisError> {
        let nodes = moments.node_keys(grid, entry).unwrap_or_default();
        let mut g = [0.0; 3];
        for (a, node) in nodes.iter().enumerate() {
            let q = extension.position(*node).ok_or_else(|| missing(*node))?;
            for (k, gk) in g.iter_mut().enumerate() {
                *gk += 2.0 / h * eval(&shape_derivative(a, k), xi) * values[q];
            }
        }
        Ok(g)
    };
    let mut samples = Vec::with_capacity(moments.cells.len());
    for (entry, c) in moments.cells.iter().enumerate() {
        let xc = centroid(&c.m);
        let lo = grid.cell_min(grid.cell_ijk(c.cell).ok_or_else(|| missing(c.cell))?);
        let g = cell_gradient(entry, xc)?;
        samples.push(Sample {
            x: [0, 1, 2].map(|a| lo[a] + 0.5 * h * (xc[a] + 1.0)),
            w: c.m[0].max(0.0),
            s: [g[0], g[1], g[2], 0.0, 0.0, 0.0],
        });
    }
    let nodal = crate::recover::nodal_values(grid, moments, extension, &samples);
    mesh.vertices
        .iter()
        .enumerate()
        .map(|(vi, p)| {
            let (entry, xi) = locate(grid, moments, *p).ok_or_else(|| AnalysisError::Internal {
                stage: "sample",
                detail: format!("vertex {vi} lies in no active cell"),
            })?;
            let nodes = moments.node_keys(grid, entry).unwrap_or_default();
            let star: Option<Vec<[f64; 6]>> = nodes
                .iter()
                .map(|n| {
                    extension
                        .position(*n)
                        .and_then(|q| nodal.get(q).copied().flatten())
                })
                .collect();
            match star {
                Some(vals) => Ok(std::array::from_fn(|k| {
                    (0..8).map(|a| eval(&shape(a), xi) * vals[a][k]).sum()
                })),
                None => cell_gradient(entry, xi),
            }
        })
        .collect()
}

/// The active entry containing `p` and `p` in its local coordinates.
fn locate(grid: &Grid, moments: &SolidMoments, p: [f64; 3]) -> Option<(usize, [f64; 3])> {
    let h = grid.cell_size();
    let [cx, cy, cz] = [0, 1, 2].map(|a| axis_candidates(grid, a, p[a]));
    for &k in &cz {
        for &j in &cy {
            for &i in &cx {
                let cell = grid.cell_index([i, j, k])?;
                if let Some(entry) = moments.entry_at(grid, cell, p, None) {
                    let lo = grid.cell_min([i, j, k]);
                    return Some((entry, [0, 1, 2].map(|a| 2.0 * (p[a] - lo[a]) / h - 1.0)));
                }
            }
        }
    }
    None
}
