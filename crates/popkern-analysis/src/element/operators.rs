//! Cell operators for one solve: every full cell has the same operator,
//! built once from the analytic moments of a full cube (exact, where a
//! full cell's rasterised moments carry the column sums' round-off);
//! cut cells and fragments get their own.

use crate::linalg::CholeskyError;
use crate::moments::{CellKind, CellMoments, moment_index};

use super::vector::{VectorCell, VectorTables};

/// The moments of a full cell of size `h`: (h/2)³ ∏ ∫ξ^k over [−1, 1].
pub fn full_cell_moments(h: f64) -> [f64; 27] {
    let one = |k: usize| {
        if k.is_multiple_of(2) {
            2.0 / (k + 1) as f64
        } else {
            0.0
        }
    };
    let s = 0.125 * h * h * h;
    let mut m = [0.0; 27];
    for c in 0..3 {
        for b in 0..3 {
            for a in 0..3 {
                m[moment_index(a, b, c)] = s * one(a) * one(b) * one(c);
            }
        }
    }
    m
}

/// A cell's operator: shared for full cells, owned otherwise.
pub enum Operator<'a> {
    /// The full-cell operator.
    Full(&'a (VectorCell, bool)),
    /// A cut cell's or fragment's own.
    Own(Box<(VectorCell, bool)>),
}

impl Operator<'_> {
    /// The operator and whether EAS was used.
    pub fn get(&self) -> &(VectorCell, bool) {
        match self {
            Operator::Full(op) => op,
            Operator::Own(op) => op,
        }
    }
}

/// How much of an entry is an infill core (a second material,
/// `structural::infill`).
#[derive(Debug, Clone, PartialEq)]
pub enum Core {
    /// None of it.
    None,
    /// All of it.
    Full,
    /// The core's part: its moments (the rest is the walls' material).
    Partial(Box<[f64; 27]>),
}

impl Core {
    /// The materials of a cell with moments `m`: (moments, λ, μ) of the
    /// walls (Lamé `walls`) and of the core (`core`, when given),
    /// summing to `m`.
    pub fn parts(
        &self,
        m: &[f64; 27],
        walls: (f64, f64),
        core: Option<(f64, f64)>,
    ) -> Vec<([f64; 27], f64, f64)> {
        match (self, core) {
            (Core::None, _) | (_, None) => vec![(*m, walls.0, walls.1)],
            (Core::Full, Some(c)) => vec![(*m, c.0, c.1)],
            (Core::Partial(mc), Some(c)) => {
                let w: [f64; 27] = std::array::from_fn(|i| m[i] - mc[i]);
                vec![(w, walls.0, walls.1), (**mc, c.0, c.1)]
            }
        }
    }
}

/// The operators of one solve's cells.
pub struct CellOperators<'a> {
    tables: &'a VectorTables,
    h: f64,
    lambda: f64,
    mu: f64,
    full: (VectorCell, bool),
    /// The core's Lamé parameters and its full-cell operator.
    core: Option<((f64, f64), (VectorCell, bool))>,
}

impl<'a> CellOperators<'a> {
    /// Builds the shared full-cell operator for cells of size `h` and
    /// Lamé parameters (λ, μ).
    pub fn new(
        tables: &'a VectorTables,
        h: f64,
        lambda: f64,
        mu: f64,
    ) -> Result<CellOperators<'a>, CholeskyError> {
        let full = tables.entry_operator(&full_cell_moments(h), h, lambda, mu, false)?;
        Ok(CellOperators {
            tables,
            h,
            lambda,
            mu,
            full,
            core: None,
        })
    }

    /// The same, with an infill core of Lamé parameters (λ_c, μ_c).
    pub fn with_core(mut self, (lambda, mu): (f64, f64)) -> Result<Self, CholeskyError> {
        let full =
            self.tables
                .entry_operator(&full_cell_moments(self.h), self.h, lambda, mu, false)?;
        self.core = Some(((lambda, mu), full));
        Ok(self)
    }

    /// [`Self::of`] for an entry with infill core `core`.
    pub fn of_core(
        &self,
        cell: &CellMoments,
        fragment: bool,
        core: &Core,
    ) -> Result<Operator<'_>, CholeskyError> {
        let Some(((lc, mc), full_core)) = &self.core else {
            return self.of(cell, fragment);
        };
        match core {
            Core::None => self.of(cell, fragment),
            Core::Full if cell.kind == CellKind::Full && !fragment => Ok(Operator::Full(full_core)),
            Core::Full => {
                let own = self
                    .tables
                    .entry_operator(&cell.m, self.h, *lc, *mc, fragment)?;
                Ok(Operator::Own(Box::new(own)))
            }
            Core::Partial(m_core) => {
                let walls: [f64; 27] = std::array::from_fn(|i| cell.m[i] - m_core[i]);
                let parts = [(&walls, self.lambda, self.mu), (&**m_core, *lc, *mc)];
                let use_eas = cell.m[0] / (self.h * self.h * self.h) >= super::vector::ETA_EAS;
                let enhancement = if use_eas {
                    super::vector::Enhancement::FragmentMean
                } else {
                    super::vector::Enhancement::None
                };
                let op = self
                    .tables
                    .cell_parts(&cell.m, &parts, self.h, enhancement)?;
                Ok(Operator::Own(Box::new((op, use_eas))))
            }
        }
    }

    /// The operator of `cell` (`fragment` when it is one of a split
    /// cell's entries).
    pub fn of(&self, cell: &CellMoments, fragment: bool) -> Result<Operator<'_>, CholeskyError> {
        if cell.kind == CellKind::Full && !fragment {
            return Ok(Operator::Full(&self.full));
        }
        let own = self
            .tables
            .entry_operator(&cell.m, self.h, self.lambda, self.mu, fragment)?;
        Ok(Operator::Own(Box::new(own)))
    }
}
