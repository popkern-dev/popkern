//! Wall stresses inside an infilled part (PH-113): the walls meet the
//! soft core on an internal material boundary, and their stress can
//! peak there, away from every surface sample. Measured on a pulled
//! infilled bar: the inner corner of the end wall reads σ₁ 18% above
//! the surface peak.
//!
//! Sampled are full wall cells (no core, or mostly wall) at the
//! interface: holding some core themselves or next to a cell that does.
//! Each is read at its centroid and corners from its own field σ_h with
//! the walls' material; SPR's patches mix both materials across the
//! interface, so the recovered σ* is not used here.

use crate::element::operators::Core;
use crate::grid::Grid;
use crate::moments::{CellKind, SolidMoments};

use super::cell::{StressPoly, centroid, stress_at};
use super::tensor::{principal, von_mises};
use super::{Recovery, VertexStress};

/// A wall sample inside the part.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InteriorSample {
    /// Position, mm.
    pub position: [f64; 3],
    /// The walls' stress there (`Recovery::InCell`).
    pub stress: VertexStress,
}

fn has_core(core: &Core) -> bool {
    !matches!(core, Core::None)
}

/// The interface wall samples of one solid.
pub fn interior_wall_samples(
    grid: &Grid,
    moments: &SolidMoments,
    cores: &[Core],
    polys: &[StressPoly],
    cell_error: &[Option<f64>],
) -> Vec<InteriorSample> {
    let h = grid.cell_size();
    let dims = grid.dims();
    let mut out = Vec::new();
    for (entry, (cell, core)) in moments.cells.iter().zip(cores).enumerate() {
        let wall = match core {
            Core::None => true,
            Core::Full => false,
            Core::Partial(m) => m[0] <= 0.5 * cell.m[0],
        };
        if !wall || cell.kind != CellKind::Full || moments.is_fragment(entry) {
            continue;
        }
        let Some(ijk) = grid.cell_ijk(cell.cell) else {
            continue;
        };
        let at_interface = has_core(core)
            || (-1i64..=1).any(|dz| {
                (-1i64..=1).any(|dy| {
                    (-1i64..=1).any(|dx| {
                        let q = [ijk[0] as i64 + dx, ijk[1] as i64 + dy, ijk[2] as i64 + dz];
                        if (0..3).any(|a| q[a] < 0 || q[a] >= dims[a] as i64) {
                            return false;
                        }
                        grid.cell_index(q.map(|v| v as usize)).is_some_and(|c| {
                            moments
                                .entries(c)
                                .any(|e| cores.get(e).is_some_and(has_core))
                        })
                    })
                })
            });
        if !at_interface {
            continue;
        }
        let lo = grid.cell_min(ijk);
        let corners =
            (0..8usize).map(|c| [c & 1, (c >> 1) & 1, (c >> 2) & 1].map(|b| 2.0 * b as f64 - 1.0));
        for xi in std::iter::once(centroid(&cell.m)).chain(corners) {
            let s = stress_at(&polys[entry], xi);
            out.push(InteriorSample {
                position: [0, 1, 2].map(|a| lo[a] + 0.5 * h * (xi[a] + 1.0)),
                stress: VertexStress {
                    stress: s,
                    von_mises: von_mises(&s),
                    max_principal: principal(&s)[0],
                    recovery: Recovery::InCell,
                    error: cell_error.get(entry).copied().flatten(),
                },
            });
        }
    }
    out
}
