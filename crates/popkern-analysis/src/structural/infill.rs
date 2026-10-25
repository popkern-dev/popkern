//! Gyroid infill as solid walls around a homogenised core (map Stage 2,
//! item 8; notes 11, Bean, Lopez-Anido & Vel 2022).
//!
//! The walls are the material within the wall thickness of the surface
//! (the perimeters); the rest is the core, a homogenised gyroid of
//! relative density D with Bean's fits for E and G (`gyroid`; nearly
//! isotropic, Zener 0.95–1.25). The core is the set at least the wall thickness
//! from the surface, the inward offset popkern-shell would build, found
//! here by distance: an entry farther than the wall thickness (plus its
//! half-diagonal) is all core, a nearer one has none, and one between
//! takes its core moments from an 8³ lattice. Cut cells and fragments
//! are all walls (walls at least a cell thick put them there).

use popkern_math::PortableMath;

use crate::clip::ClippedSolid;
use crate::element::operators::Core;
use crate::grid::Grid;
use crate::moments::{CellKind, SolidMoments, moment_index};

/// An infill: relative density and wall thickness.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Infill {
    /// Relative density D of the gyroid core, in (0, 1].
    pub density: f64,
    /// Solid wall thickness, mm (perimeters × line width).
    pub walls_mm: f64,
    /// The gyroid's wall (bead) width h, mm.
    pub line_width_mm: f64,
}

/// The gyroid's surface area per unit-cell volume times λ: a thin
/// gyroid wall of width h fills D ≈ 3.09 h/λ of its cell (Bean 2025
/// §2.1: D ∝ h/λ). Bean 2022's double-bead 25% blocks, 6 cells across
/// 50.8 mm, give 3.0 with Cura's 0.35 mm line width.
pub const GYROID_AREA_PER_CELL: f64 = 3.09;

impl Infill {
    /// The gyroid's unit-cell size λ = 3.09 h/D, mm: the scale its
    /// homogenised properties, and so its strength, average over.
    pub fn cell_mm(&self) -> f64 {
        GYROID_AREA_PER_CELL * self.line_width_mm / self.density
    }
}

/// Bean 2022's homogenised gyroid of relative density `d` in a material
/// (E_m, ν_m), as an isotropic (E, ν). The gyroid is cubic; E and G
/// are Bean's (eqns 12–13, errors 0.44% and 0.35%):
/// E/E_m = D − (1 − D)^0.77 D^1.2, G/G_m = E/E_m + 0.5 ν_m D² sin(πD),
/// and ν = E/2G − 1 keeps both. Bean's own ν (eqn 14) is the cubic
/// one, which passes 0.5 at low density (0.50 at D = 0.25, ν_m = 0.35)
/// and so has no isotropic counterpart.
pub fn gyroid(d: f64, youngs: f64, poisson: f64) -> (f64, f64) {
    let e_ratio = d - (1.0 - d).ppowf(0.77) * d.ppowf(1.2);
    let g_ratio = e_ratio + 0.5 * poisson * d * d * (std::f64::consts::PI * d).psin();
    let g = g_ratio * youngs / (2.0 * (1.0 + poisson));
    let e = e_ratio * youngs;
    (e, e / (2.0 * g) - 1.0)
}

/// Bean, Vel & Lopez-Anido 2025 (eqn 8): the gyroid's effective yield
/// strength over its parent's, S_g/S_m = D − (1 − D)^C₁ D^C₂, with
/// C₁ = 0.67, C₂ = 1.2 in compression (RMS error 2.2%) and C₁ = C₂ = 1.3
/// in shear, normalised by the parent's shear yield (1.4%). Returns
/// (compressive, shear) ratios; both fall to 0 at D = 0 and rise to 1
/// at D = 1.
pub fn gyroid_yield(d: f64) -> (f64, f64) {
    let ratio = |c1: f64, c2: f64| d - (1.0 - d).ppowf(c1) * d.ppowf(c2);
    (ratio(0.67, 1.2), ratio(1.3, 1.3))
}

/// Bean 2025's scatter of printed gyroid yield strengths, as a fraction
/// of the strength: the average standard deviation of its compression
/// tests (3.6%), larger than either fit's error.
pub const GYROID_YIELD_SCATTER: f64 = 0.036;

/// The core's yield strength for a parent of yield strength `parent`
/// (MPa), as one isotropic von Mises allowable: the lower of the
/// compressive ratio and the shear ratio (the shear fit is normalised
/// by the parent's von Mises shear yield, S_m/√3, so S_m times its ratio
/// is the matching uniaxial value). The homogenised core is isotropic,
/// so the gyroid's axes are unknown here and the weaker mode governs.
pub fn core_yield(d: f64, parent: f64) -> f64 {
    let (compressive, shear) = gyroid_yield(d);
    parent * compressive.min(shear)
}

/// The core's volume (model units cubed) of one solid's entries.
pub fn core_volume(moments: &SolidMoments, cores: Option<&[Core]>) -> f64 {
    let Some(cores) = cores else {
        return 0.0;
    };
    moments
        .cells
        .iter()
        .zip(cores)
        .map(|(cell, core)| match core {
            Core::None => 0.0,
            Core::Full => cell.m[0],
            Core::Partial(m) => m[0],
        })
        .sum()
}

/// Lattice points per cell side for a partial core.
const LATTICE: usize = 8;

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Distance from `p` to the convex planar polygon `v` (normal `n`).
pub(crate) fn polygon_distance(p: [f64; 3], v: &[[f64; 3]], n: [f64; 3]) -> f64 {
    let Some(&a) = v.first() else {
        return f64::INFINITY;
    };
    let h = dot(sub(p, a), n);
    let q = [0, 1, 2].map(|k| p[k] - h * n[k]);
    let cross = |x: [f64; 3], y: [f64; 3]| {
        [
            x[1] * y[2] - x[2] * y[1],
            x[2] * y[0] - x[0] * y[2],
            x[0] * y[1] - x[1] * y[0],
        ]
    };
    let side = |i: usize| {
        let (s, e) = (v[i], v[(i + 1) % v.len()]);
        dot(cross(sub(e, s), sub(q, s)), n)
    };
    if (0..v.len()).all(|i| side(i) >= 0.0) || (0..v.len()).all(|i| side(i) <= 0.0) {
        return h.abs();
    }
    // Nearest point on the boundary.
    (0..v.len())
        .map(|i| {
            let (s, e) = (v[i], v[(i + 1) % v.len()]);
            let d = sub(e, s);
            let t = (dot(sub(p, s), d) / dot(d, d).max(f64::MIN_POSITIVE)).clamp(0.0, 1.0);
            let x = [0, 1, 2].map(|k| s[k] + t * d[k]);
            dot(sub(p, x), sub(p, x)).sqrt()
        })
        .fold(f64::INFINITY, f64::min)
}

/// Each entry's core for walls `walls` (mm).
pub fn cores(grid: &Grid, clipped: &ClippedSolid, moments: &SolidMoments, walls: f64) -> Vec<Core> {
    let h = grid.cell_size();
    let reach = (walls / h).ceil() as i64 + 1;
    let half_diagonal = 0.5 * h * 3.0f64.sqrt();
    let pieces = &clipped.pieces;
    let dims = grid.dims();
    let mut out = Vec::with_capacity(moments.cells.len());
    for (entry, cell) in moments.cells.iter().enumerate() {
        if cell.kind != CellKind::Full || moments.is_fragment(entry) {
            out.push(Core::None);
            continue;
        }
        let Some(ijk) = grid.cell_ijk(cell.cell) else {
            out.push(Core::None);
            continue;
        };
        // The surface pieces within reach.
        let mut near: Vec<&crate::clip::FacetPiece> = Vec::new();
        for dz in -reach..=reach {
            for dy in -reach..=reach {
                for dx in -reach..=reach {
                    let q = [ijk[0] as i64 + dx, ijk[1] as i64 + dy, ijk[2] as i64 + dz];
                    if (0..3).any(|a| q[a] < 0 || q[a] >= dims[a] as i64) {
                        continue;
                    }
                    let Some(c) = grid.cell_index(q.map(|v| v as usize)) else {
                        continue;
                    };
                    let start = pieces.partition_point(|p| p.cell < c);
                    let end = pieces.partition_point(|p| p.cell <= c);
                    near.extend(&pieces[start..end]);
                }
            }
        }
        let distance = |p: [f64; 3]| {
            near.iter()
                .map(|piece| polygon_distance(p, piece.polygon.vertices(), piece.normal))
                .fold(f64::INFINITY, f64::min)
        };
        let lo = grid.cell_min(ijk);
        let centre = lo.map(|x| x + 0.5 * h);
        let dc = distance(centre);
        if dc >= walls + half_diagonal {
            out.push(Core::Full);
            continue;
        }
        if dc + half_diagonal < walls {
            out.push(Core::None);
            continue;
        }
        let s = h / LATTICE as f64;
        let w = s * s * s;
        let mut m = [0.0; 27];
        for k in 0..LATTICE {
            for j in 0..LATTICE {
                for i in 0..LATTICE {
                    let xi = [i, j, k].map(|q| -1.0 + (2.0 * q as f64 + 1.0) / LATTICE as f64);
                    let p = [0, 1, 2].map(|a| lo[a] + 0.5 * h * (xi[a] + 1.0));
                    if distance(p) < walls {
                        continue;
                    }
                    for c in 0..3 {
                        for b in 0..3 {
                            for a in 0..3 {
                                m[moment_index(a, b, c)] += w
                                    * xi[0].powi(a as i32)
                                    * xi[1].powi(b as i32)
                                    * xi[2].powi(c as i32);
                            }
                        }
                    }
                }
            }
        }
        out.push(if m[0] > 0.0 {
            Core::Partial(Box::new(m))
        } else {
            Core::None
        });
    }
    out
}

#[cfg(test)]
#[path = "infill_tests.rs"]
mod tests;
