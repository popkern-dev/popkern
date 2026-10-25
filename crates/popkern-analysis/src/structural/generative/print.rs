//! The layer-by-layer print filter (Langelaar 2017, "An additive
//! manufacturing filter for topology optimization of print-ready
//! designs"): a design is printed as built, each layer's material at
//! most what the layer below holds up, so an unsupported overhang never
//! appears in the printed densities the analysis sees.
//!
//! Layers are grid layers along the build direction, the first (on the
//! bed) the lowest the part occupies. An entry's support region is the
//! cell below it and that cell's four lateral face neighbours (the
//! paper's 5-element region in 3D), so a wall may lean 45° along a grid
//! axis and about 35° along a diagonal: conservative for FDM's usual
//! 45° limit. The printed density is ξ = smin(x, Ξ) with
//! Ξ = smax(ξ over the support region) (eqs 7–8), smoothed by eq 3,
//! smin = ½(x + Ξ − √((x − Ξ)² + ε) + √ε), and eq 5, the P-norm with
//! exponent 1/Q, Q = P + ln n_S / ln ξ₀ (eq 6), which removes the
//! overshoot of uniform layers of density ξ₀. Material the design must
//! keep is printed regardless (the part is printed with whatever
//! supports it needs there) and holds up what is above it.
//!
//! Sensitivities go back by the paper's adjoint (eqs 13–14): one sweep
//! from the top layer down.

use popkern_math::PortableMath;

use super::filter::DesignEntry;

/// Paper defaults (§3): ε = 1e-4, P = 40, ξ₀ = 0.5.
const EPSILON: f64 = 1e-4;
const P: f64 = 40.0;
const XI0: f64 = 0.5;
/// Cells in the support region.
const REGION: f64 = 5.0;
/// An entry below this fraction of its cell is empty.
const EMPTY: f64 = 1e-12;

/// The build direction: the bed's normal, along a grid axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub struct BuildDirection {
    /// Axis: 0 x, 1 y, 2 z.
    pub axis: usize,
    /// Whether the part grows towards −axis.
    pub negative: bool,
}

/// The print filter over a design's entries.
#[derive(Debug, Clone, PartialEq)]
pub struct PrintFilter {
    /// Entries in sweep order: by layer, then entry index.
    order: Vec<usize>,
    /// Each entry's supporting entries (the layer below).
    support: Vec<Vec<usize>>,
    /// Whether each entry lies on the bed.
    base: Vec<bool>,
    q: f64,
}

/// A printed design with what the adjoint needs.
#[derive(Debug, Clone, PartialEq)]
pub struct Printed {
    /// Printed densities ξ (at least the floors; an empty entry's is its
    /// floor).
    pub densities: Vec<f64>,
    /// Ξ per entry (the support's smax; 0 on the bed).
    support: Vec<f64>,
    /// The support's largest density m and R = Σ (ξ/m)ᴾ ≥ 1, so that
    /// Ξ = m^(P/Q) R^(1/Q): Σ ξᴾ itself underflows when the support is
    /// nearly empty, and its power 1/Q − 1 in the adjoint then overflowed
    /// to infinite sensitivities.
    peaks: Vec<(f64, f64)>,
    /// Whether the floor holds the entry.
    held: Vec<bool>,
}

impl Printed {
    /// Whether the floor holds each entry.
    pub fn held(&self) -> &[bool] {
        &self.held
    }
}

fn smin(x: f64, s: f64) -> f64 {
    0.5 * (x + s - ((x - s) * (x - s) + EPSILON).sqrt() + EPSILON.sqrt())
}

/// (∂smin/∂x, ∂smin/∂Ξ).
fn smin_slopes(x: f64, s: f64) -> (f64, f64) {
    let r = (x - s) / ((x - s) * (x - s) + EPSILON).sqrt();
    (0.5 * (1.0 - r), 0.5 * (1.0 + r))
}

impl PrintFilter {
    /// The filter of `entries` on a grid of cell size `h`, built along
    /// `direction`.
    pub fn new(entries: &[DesignEntry], h: f64, direction: BuildDirection) -> PrintFilter {
        let n = entries.len();
        let origin = entries.first().map_or([0.0; 3], |e| e.corner);
        let key = |e: &DesignEntry| -> [i64; 3] {
            let mut k = [0, 1, 2].map(|a| ((e.corner[a] - origin[a]) / h).round() as i64);
            if direction.negative {
                k[direction.axis] = -k[direction.axis];
            }
            k
        };
        // (cell key, entry), sorted: lookups by binary search. Entries
        // without volume (cells that only touch the part) are not
        // material: they neither print nor hold anything up.
        let empty = |e: &DesignEntry| e.volume <= EMPTY * h * h * h;
        let mut cells: Vec<([i64; 3], usize)> = entries
            .iter()
            .enumerate()
            .filter(|(_, e)| !empty(e))
            .map(|(i, e)| (key(e), i))
            .collect();
        cells.sort_unstable();
        let lowest = cells.iter().map(|c| c.0[direction.axis]).min().unwrap_or(0);
        let at = |k: [i64; 3]| -> &[([i64; 3], usize)] {
            let from = cells.partition_point(|c| c.0 < k);
            let to = cells.partition_point(|c| c.0 <= k);
            &cells[from..to]
        };
        let (a, b) = ((direction.axis + 1) % 3, (direction.axis + 2) % 3);
        let mut support = vec![Vec::new(); n];
        let mut base = vec![false; n];
        for &(k, i) in &cells {
            if k[direction.axis] == lowest {
                base[i] = true;
                continue;
            }
            let mut below = k;
            below[direction.axis] -= 1;
            let mut region = vec![below];
            for (axis, step) in [(a, -1), (a, 1), (b, -1), (b, 1)] {
                let mut c = below;
                c[axis] += step;
                region.push(c);
            }
            for c in region {
                support[i].extend(at(c).iter().map(|&(_, j)| j));
            }
        }
        let mut order: Vec<usize> = cells.iter().map(|&(_, i)| i).collect();
        order.sort_by_key(|&i| (key(&entries[i])[direction.axis], i));
        PrintFilter {
            order,
            support,
            base,
            q: P + REGION.pln() / XI0.pln(),
        }
    }

    /// The printed design of blueprint `x` with density floors `floors`.
    pub fn apply(&self, x: &[f64], floors: &[f64]) -> Printed {
        let n = x.len();
        let mut out = Printed {
            densities: floors.to_vec(),
            support: vec![0.0; n],
            peaks: vec![(0.0, 0.0); n],
            held: vec![false; n],
        };
        for &i in &self.order {
            let printed = if self.base[i] {
                x[i]
            } else {
                let m = self.support[i]
                    .iter()
                    .map(|&j| out.densities[j].max(0.0))
                    .fold(0.0, f64::max);
                let (s, r) = if m > 0.0 {
                    let r: f64 = self.support[i]
                        .iter()
                        .map(|&j| (out.densities[j].max(0.0) / m).ppowf(P))
                        .sum();
                    (m.ppowf(P / self.q) * r.ppowf(1.0 / self.q), r)
                } else {
                    (0.0, 0.0)
                };
                out.peaks[i] = (m, r);
                out.support[i] = s;
                smin(x[i], s)
            };
            out.held[i] = printed < floors[i];
            out.densities[i] = printed.max(floors[i]);
        }
        out
    }

    /// The fraction of a design's solid volume that would print
    /// unsupported: entries above `iso` (with `volumes`), printed by the
    /// exact rule (eqs 1–2: on the bed, or over a printed entry of the
    /// support region), kept entries (`floors` at 1) printing
    /// regardless.
    pub fn unsupported(&self, densities: &[f64], floors: &[f64], volumes: &[f64], iso: f64) -> f64 {
        let n = densities.len();
        let solid: Vec<bool> = densities.iter().map(|&d| d > iso).collect();
        let mut printed = vec![false; n];
        let (mut total, mut missing) = (0.0, 0.0);
        for &i in &self.order {
            let kept = floors[i] >= 1.0;
            printed[i] =
                kept || (solid[i] && (self.base[i] || self.support[i].iter().any(|&j| printed[j])));
            if solid[i] && !kept {
                total += volumes[i];
                if !printed[i] {
                    missing += volumes[i];
                }
            }
        }
        if total > 0.0 { missing / total } else { 0.0 }
    }

    /// ∂f/∂x from ∂f/∂ξ (`g`), for the printing of `x` into `printed`.
    pub fn back(&self, x: &[f64], printed: &Printed, g: &[f64]) -> Vec<f64> {
        let n = x.len();
        // λ_i: ∂f/∂ξ_i, the entries above it included, built top down.
        let mut lambda: Vec<f64> = g.to_vec();
        let mut out = vec![0.0; n];
        for &i in self.order.iter().rev() {
            if printed.held[i] {
                continue;
            }
            if self.base[i] {
                out[i] = lambda[i];
                continue;
            }
            let (dx, ds) = smin_slopes(x[i], printed.support[i]);
            out[i] = lambda[i] * dx;
            let (m, r) = printed.peaks[i];
            if m <= 0.0 {
                continue;
            }
            // ∂Ξ/∂ξ_j = (P/Q) ξ_jᴾ⁻¹ (Σ ξᴾ)^(1/Q − 1) (eq 14)
            // = (P/Q) m^(P/Q − 1) (ξ_j/m)ᴾ⁻¹ R^(1/Q − 1), each factor
            // bounded.
            let common = lambda[i]
                * ds
                * (P / self.q)
                * m.ppowf(P / self.q - 1.0)
                * r.ppowf(1.0 / self.q - 1.0);
            for &j in &self.support[i] {
                let xi = printed.densities[j].max(0.0) / m;
                lambda[j] += common * xi.ppowf(P - 1.0);
            }
        }
        out
    }
}

#[cfg(test)]
#[path = "print_tests.rs"]
mod tests;
