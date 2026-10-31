//! Axis-aligned boxes for the table's spatial filters: a curve's reach
//! for the interior-incidence scan, and the interference sites that
//! bound which operand vertices can interfere at all.

use popkern_core::geometry::Point3;

/// An axis-aligned box.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Bounds {
    /// Lower corner.
    pub lo: [f64; 3],
    /// Upper corner.
    pub hi: [f64; 3],
}

impl Bounds {
    /// The box of `points` grown by `grow` on every side; `None` when
    /// there are no points.
    pub(crate) fn of(points: &[Point3], grow: f64) -> Option<Self> {
        let first = points.first()?;
        let mut b = Bounds {
            lo: [first.x, first.y, first.z],
            hi: [first.x, first.y, first.z],
        };
        for p in points {
            for (i, c) in [p.x, p.y, p.z].into_iter().enumerate() {
                b.lo[i] = b.lo[i].min(c);
                b.hi[i] = b.hi[i].max(c);
            }
        }
        Some(b.grown(grow))
    }

    /// This box grown by `grow` on every side.
    pub(crate) fn grown(mut self, grow: f64) -> Self {
        for i in 0..3 {
            self.lo[i] -= grow;
            self.hi[i] += grow;
        }
        self
    }

    /// Whether `p` lies in the box.
    pub(crate) fn holds(&self, p: &Point3) -> bool {
        [p.x, p.y, p.z]
            .into_iter()
            .enumerate()
            .all(|(i, c)| c >= self.lo[i] && c <= self.hi[i])
    }
}
