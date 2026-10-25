//! Evaluating the recovered field at a point of a cell: trilinear
//! interpolation of nodal σ*, or the cell's own EAS field where a node
//! has no σ* (step 3 of the module docs). Used for display vertices and
//! for the cell-resolution surface samples that peaks are searched on.

use serde::Serialize;

use crate::clip::ClippedSolid;
use crate::element::poly::{Poly27, eval};
use crate::model::FaceId;

use super::cell::{StressPoly, stress_at};
use super::normal::FlatCells;
use super::spr::Sample;
use super::tensor::{Voigt, principal, von_mises};
use super::{Fallback, Recovery, VertexStress};

/// Finds (entry, local coordinates) of a point inside the solid, `None`
/// outside it.
pub(super) type Inside<'a> = dyn Fn([f64; 3]) -> Option<(usize, [f64; 3])> + 'a;

/// Normal extrapolation at flat surfaces (`super::normal`).
pub(super) struct Normal<'a> {
    /// Flat cells.
    pub flat: &'a FlatCells,
    /// Locates interior points.
    pub inside: &'a Inside<'a>,
    /// Cell size.
    pub h: f64,
}

/// Everything needed to evaluate σ* in any active cell.
pub(super) struct Sampler<'a> {
    pub shapes: &'a [Poly27; 8],
    pub cell_nodes: &'a [[usize; 8]],
    pub nodal: &'a [Option<Voigt>],
    pub polys: &'a [StressPoly],
    pub samples: &'a [Sample],
    pub cell_error: &'a [Option<f64>],
    pub fallback: Fallback,
    pub normal: Option<Normal<'a>>,
}

impl Sampler<'_> {
    /// Nodal σ* interpolated at local point `xi` of cell `k`, where all
    /// its nodes have one.
    fn interpolated(&self, k: usize, xi: [f64; 3]) -> Option<Voigt> {
        let vals = (0..8)
            .map(|a| self.nodal[self.cell_nodes[k][a]])
            .collect::<Option<Vec<Voigt>>>()?;
        Some(std::array::from_fn(|q| {
            (0..8).map(|a| eval(&self.shapes[a], xi) * vals[a][q]).sum()
        }))
    }

    /// The recovered stress at surface point `x` (local `xi` of cell
    /// `k`, in grid cell `cell`, unit inward normal `inward`): normal
    /// extrapolation where the surface is flat and both inward points
    /// are inside the solid with σ*, else [`Sampler::at`].
    pub fn at_surface(
        &self,
        (k, xi): (usize, [f64; 3]),
        (x, cell, inward): ([f64; 3], usize, [f64; 3]),
        thin: bool,
    ) -> VertexStress {
        let extrapolated = self
            .normal
            .as_ref()
            .filter(|_| !thin && self.fallback == Fallback::InCell)
            .filter(|nm| nm.flat.is_flat(cell, x))
            .and_then(|nm| {
                let value = |t: f64| {
                    let q = [0, 1, 2].map(|a| x[a] + t * nm.h * inward[a]);
                    let (kq, xq) = (nm.inside)(q)?;
                    self.interpolated(kq, xq)
                };
                let (s1, s2) = (value(1.0)?, value(2.0)?);
                Some(std::array::from_fn::<f64, 6, _>(|q| 2.0 * s1[q] - s2[q]))
            });
        match extrapolated {
            Some(stress) => VertexStress {
                stress,
                von_mises: von_mises(&stress),
                max_principal: principal(&stress)[0],
                recovery: Recovery::Patch,
                error: self.cell_error[k],
            },
            None => self.at(k, xi, thin),
        }
    }

    /// The recovered stress at local point `xi` of active cell `k`; no
    /// error estimate when `thin` (a wall under the threshold) or when
    /// the in-cell field is used.
    pub fn at(&self, k: usize, xi: [f64; 3], thin: bool) -> VertexStress {
        let interpolated = self.interpolated(k, xi);
        let (stress, recovery) = match (interpolated, self.fallback) {
            (_, Fallback::ForceConstant) => (self.samples[k].s, Recovery::InCell),
            (_, Fallback::ForceInCell) | (None, Fallback::InCell) => {
                (stress_at(&self.polys[k], xi), Recovery::InCell)
            }
            (Some(s), Fallback::InCell) => (s, Recovery::Patch),
        };
        let error = match recovery {
            Recovery::Patch if !thin => self.cell_error[k],
            _ => None,
        };
        VertexStress {
            stress,
            von_mises: von_mises(&stress),
            max_principal: principal(&stress)[0],
            recovery,
            error,
        }
    }
}

/// The recovered stress at one point of the surface, at cell
/// resolution: a facet piece's centroid (weighted by its area) or one
/// of its corners (weight 0).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SurfaceSample {
    /// The face the piece lies on.
    pub face: FaceId,
    /// Position, mm.
    pub position: [f64; 3],
    /// Surface area the sample stands for, mm² (0 for corners).
    pub area: f64,
    /// The stress there.
    pub stress: VertexStress,
}

/// Finds (entry, local coordinates) of a point of clipped piece
/// `index` in cell `cell`.
pub(super) type Locate<'a> = dyn Fn(usize, usize, [f64; 3]) -> Option<(usize, [f64; 3])> + 'a;

/// Samples at every facet piece's centroid and corners. `thin[t]` marks
/// display triangles on a wall under the threshold.
pub(super) fn surface_samples(
    sampler: &Sampler<'_>,
    clipped: &ClippedSolid,
    locate: &Locate<'_>,
    thin: &[bool],
) -> Vec<SurfaceSample> {
    let mut out = Vec::new();
    for (index, piece) in clipped.pieces.iter().enumerate() {
        let v = piece.polygon.vertices();
        if v.is_empty() {
            continue;
        }
        let n = v.len() as f64;
        let c = [0, 1, 2].map(|a| v.iter().map(|p| p[a]).sum::<f64>() / n);
        let thin = thin.get(piece.triangle as usize).copied().unwrap_or(false);
        let points = std::iter::once((c, piece.area)).chain(v.iter().map(|&p| (p, 0.0)));
        for (p, area) in points {
            if let Some(at) = locate(index, piece.cell, p) {
                let inward = piece.normal.map(|n| -n);
                out.push(SurfaceSample {
                    face: piece.face,
                    position: p,
                    area,
                    stress: sampler.at_surface(at, (p, piece.cell, inward), thin),
                });
            }
        }
    }
    out
}
