//! Exact midline subdivision of polyhedral-spline Bézier patches.
//!
//! The TOMS mask tables tile per quarter-cell for most families but
//! per half-quarter for stars and n-gons of valence/sides >= 6. These
//! helpers subdivide patches at their parametric midlines (de
//! Casteljau, exact) so mixed complexes can be brought to one
//! boundary granularity and welded 1:1. See
//! `docs/architecture/smooth-mesh-cage-normalization.md`.

use super::types::{PolyhedralPatch, PolyhedralSpline, PolyhedralSplineError};
use crate::geometry::Point3;

/// Exactly subdivide every patch of a spline at its parametric
/// midlines (each patch becomes four). Used to bring a coarse-tiled
/// component to the half-quarter boundary granularity of the
/// valence/sides ≥ 6 star and n-gon tables, so mixed complexes weld
/// 1:1 (see `complex.rs`).
pub(super) fn subdivide_spline_patch_grid(
    spline: &PolyhedralSpline,
) -> Result<PolyhedralSpline, PolyhedralSplineError> {
    let mut patches = Vec::with_capacity(spline.patches.len() * 4);
    for patch in &spline.patches {
        patches.extend(quarter_subdivided(patch)?);
    }
    Ok(PolyhedralSpline {
        family: spline.family,
        patches,
    })
}

/// Split a tensor-product Bézier patch at (1/2, 1/2) into four
/// quarter patches covering [0,.5]x[0,.5] … [.5,1]x[.5,1]. De
/// Casteljau subdivision: the four patches reproduce the input
/// surface exactly.
pub(super) fn quarter_subdivided(
    patch: &PolyhedralPatch,
) -> Result<Vec<PolyhedralPatch>, PolyhedralSplineError> {
    let (lo_u, hi_u) = split_u(patch)?;
    let mut out = Vec::with_capacity(4);
    for half in [lo_u, hi_u] {
        let (lo_v, hi_v) = split_v(&half)?;
        out.push(lo_v);
        out.push(hi_v);
    }
    Ok(out)
}

fn split_u(
    patch: &PolyhedralPatch,
) -> Result<(PolyhedralPatch, PolyhedralPatch), PolyhedralSplineError> {
    let nv = patch.degree_v + 1;
    let mut lo = Vec::with_capacity((patch.degree_u + 1) * nv);
    let mut hi = Vec::with_capacity((patch.degree_u + 1) * nv);
    let column: Vec<Vec<Point3>> = (0..nv)
        .map(|j| {
            (0..=patch.degree_u)
                .map(|i| patch.cp(i, j))
                .collect::<Vec<_>>()
        })
        .collect();
    let split: Vec<(Vec<Point3>, Vec<Point3>)> = column
        .iter()
        .map(|points| de_casteljau_halves(points))
        .collect();
    for i in 0..=patch.degree_u {
        for half_column in &split {
            lo.push(half_column.0[i]);
            hi.push(half_column.1[i]);
        }
    }
    Ok((
        PolyhedralPatch::new(patch.family, patch.degree_u, patch.degree_v, lo)?,
        PolyhedralPatch::new(patch.family, patch.degree_u, patch.degree_v, hi)?,
    ))
}

fn split_v(
    patch: &PolyhedralPatch,
) -> Result<(PolyhedralPatch, PolyhedralPatch), PolyhedralSplineError> {
    let mut lo = Vec::with_capacity((patch.degree_u + 1) * (patch.degree_v + 1));
    let mut hi = Vec::with_capacity((patch.degree_u + 1) * (patch.degree_v + 1));
    for i in 0..=patch.degree_u {
        let row: Vec<Point3> = (0..=patch.degree_v).map(|j| patch.cp(i, j)).collect();
        let (row_lo, row_hi) = de_casteljau_halves(&row);
        lo.extend(row_lo);
        hi.extend(row_hi);
    }
    Ok((
        PolyhedralPatch::new(patch.family, patch.degree_u, patch.degree_v, lo)?,
        PolyhedralPatch::new(patch.family, patch.degree_u, patch.degree_v, hi)?,
    ))
}

/// De Casteljau split of a Bézier control polygon at t = 1/2: the
/// left half's points are the first points of each elimination
/// round, the right half's are the last points in reverse.
fn de_casteljau_halves(points: &[Point3]) -> (Vec<Point3>, Vec<Point3>) {
    let mut work = points.to_vec();
    let mut left = Vec::with_capacity(points.len());
    let mut right = Vec::with_capacity(points.len());
    left.push(work[0]);
    right.push(work[work.len() - 1]);
    while work.len() > 1 {
        let mut next = Vec::with_capacity(work.len() - 1);
        for pair in work.windows(2) {
            next.push(Point3::new(
                (pair[0].x + pair[1].x) * 0.5,
                (pair[0].y + pair[1].y) * 0.5,
                (pair[0].z + pair[1].z) * 0.5,
            ));
        }
        left.push(next[0]);
        right.push(next[next.len() - 1]);
        work = next;
    }
    right.reverse();
    (left, right)
}
