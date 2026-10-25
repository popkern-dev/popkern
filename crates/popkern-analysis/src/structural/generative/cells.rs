//! A design's entries binned to their grid cells: a dense block of
//! cells over the design's bounding box with a zero layer around it,
//! each cell's fill (its fragments' volumes over the cell's) and
//! density (the fragments' densities weighted by volume).

use super::DesignEntry;

/// The design's cells.
pub(super) struct CellField {
    /// Cell size, mm.
    pub h: f64,
    /// Cells along each axis, the zero layer included.
    pub dims: [usize; 3],
    /// The lower corner of cell (0, 0, 0), mm.
    pub corner: [f64; 3],
    /// Per cell: its fill, at most 1.
    pub fill: Vec<f64>,
    /// Per cell: its density (0 where empty).
    pub density: Vec<f64>,
}

impl CellField {
    /// The cells of `entries` at `densities` on cells of size `h`;
    /// `None` without entries or with `h` not positive.
    pub fn new(entries: &[DesignEntry], densities: &[f64], h: f64) -> Option<CellField> {
        let origin = entries.first()?.corner;
        if h <= 0.0 {
            return None;
        }
        let cell_volume = h * h * h;
        // Cell keys relative to the first entry's cell.
        let key = |c: [f64; 3]| [0, 1, 2].map(|a| ((c[a] - origin[a]) / h).round() as i64);
        let mut lo = [i64::MAX; 3];
        let mut hi = [i64::MIN; 3];
        for e in entries {
            let k = key(e.corner);
            for a in 0..3 {
                lo[a] = lo[a].min(k[a]);
                hi[a] = hi[a].max(k[a]);
            }
        }
        // Cells over [lo − 1, hi + 1].
        let dims = [0, 1, 2].map(|a| (hi[a] - lo[a] + 3) as usize);
        let at = |k: [i64; 3]| -> usize {
            let q = [0, 1, 2].map(|a| (k[a] - lo[a] + 1) as usize);
            q[0] + dims[0] * (q[1] + dims[1] * q[2])
        };
        let n = dims[0] * dims[1] * dims[2];
        let mut fill = vec![0.0f64; n];
        let mut mass = vec![0.0f64; n];
        for (e, rho) in entries.iter().zip(densities) {
            let f = (e.volume / cell_volume).max(0.0);
            let c = at(key(e.corner));
            fill[c] += f;
            mass[c] += rho * f;
        }
        let density = (0..n)
            .map(|c| {
                if fill[c] > 0.0 {
                    mass[c] / fill[c]
                } else {
                    0.0
                }
            })
            .collect();
        for f in fill.iter_mut() {
            *f = f.min(1.0);
        }
        Some(CellField {
            h,
            dims,
            corner: [0, 1, 2].map(|a| origin[a] + h * (lo[a] - 1) as f64),
            fill,
            density,
        })
    }

    /// The index of cell `q`.
    pub fn index(&self, q: [usize; 3]) -> usize {
        q[0] + self.dims[0] * (q[1] + self.dims[1] * q[2])
    }

    /// The number of cells.
    pub fn len(&self) -> usize {
        self.fill.len()
    }

    /// The centre of cell `q`, mm.
    pub fn centre(&self, q: [usize; 3]) -> [f64; 3] {
        [0, 1, 2].map(|a| self.corner[a] + self.h * (q[a] as f64 + 0.5))
    }

    /// Cell `c`'s material: the smaller of its density and its fill (a
    /// cell the part only grazes holds little whatever its density).
    pub fn material(&self, c: usize) -> f64 {
        if self.fill[c] > 0.0 {
            self.density[c].min(self.fill[c])
        } else {
            0.0
        }
    }
}
