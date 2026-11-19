//! Rigorous curvature enclosures over a parameter region (PK-510).
//!
//! Gaussian and mean curvature are enclosed cell by cell from certified
//! enclosures of the surface's partial derivatives over each parameter
//! box (`geometry::certify`: outward-rounded intervals, Bernstein hulls
//! for NURBS and Gregory exact forms, closed forms for the analytic
//! kinds). With `N = S_u × S_v`, `ℓ = S_uu·N`, `μ = S_uv·N`,
//! `ν = S_vv·N`:
//!
//! ```text
//! K = (ℓν − μ²) / |N|⁴
//! H = (|S_v|² ℓ − 2 (S_u·S_v) μ + |S_u|² ν) / (2 |N|³)
//! ```
//!
//! Each is evaluated in interval arithmetic over the cell. The result
//! contains every value the true curvature takes there, wherever the
//! enclosure of `|N|²` excludes 0.
//!
//! This replaces a sampled bound that was documented as proven but was
//! not. It took five curvature samples per cell and estimated a
//! Lipschitz constant from their differences, with a 2× safety factor,
//! so a curvature feature between the samples escaped it.
//!
//! A cell whose box cannot be enclosed stays unbounded after the last
//! subdivision: an unsupported kind, a chart singularity (`|N|²`
//! reaching 0), or a Gregory corner. It is counted in
//! [`GuaranteedCurvatureBounds::degenerate_cells`], and the bounds are
//! then the whole line, never a guess.

use crate::geometry::certify::{Iv, SurfaceJetSource};
use crate::geometry::surface::SurfaceGeometry;
use crate::healing::interval::Interval;

/// Guaranteed curvature bounds over a parameter region.
///
/// When `degenerate_cells` is 0, `gaussian` and `mean` contain every
/// value the surface's Gaussian and mean curvature take over the
/// region. Otherwise they are the whole line.
#[derive(Debug, Clone, Copy)]
pub struct GuaranteedCurvatureBounds {
    /// Interval containing all Gaussian curvature values.
    pub gaussian: Interval,
    /// Interval containing all mean curvature values.
    pub mean: Interval,
    /// Width of the Gaussian interval.
    pub gaussian_width: f64,
    /// Deepest subdivision level used.
    pub subdivisions: usize,
    /// Leaf cells that could not be enclosed (see the module docs).
    pub degenerate_cells: usize,
}

/// A rectangular cell in parameter space.
#[derive(Debug, Clone, Copy)]
struct Cell {
    u: (f64, f64),
    v: (f64, f64),
}

impl Cell {
    fn subdivide(&self) -> [Cell; 4] {
        let um = 0.5 * (self.u.0 + self.u.1);
        let vm = 0.5 * (self.v.0 + self.v.1);
        [
            Cell {
                u: (self.u.0, um),
                v: (self.v.0, vm),
            },
            Cell {
                u: (um, self.u.1),
                v: (self.v.0, vm),
            },
            Cell {
                u: (self.u.0, um),
                v: (vm, self.v.1),
            },
            Cell {
                u: (um, self.u.1),
                v: (vm, self.v.1),
            },
        ]
    }
}

/// Enclosures of `(K, H)` over one cell, or `None` when the cell cannot
/// be enclosed.
fn cell_bounds(source: &mut SurfaceJetSource, cell: &Cell) -> Option<(Iv, Iv)> {
    let (jet, _) = source
        .jet_piecewise(Iv::new(cell.u.0, cell.u.1), Iv::new(cell.v.0, cell.v.1))
        .ok()?;
    let (su, sv) = (jet.partial(1, 0), jet.partial(0, 1));
    let (suu, suv, svv) = (jet.partial(2, 0), jet.partial(1, 1), jet.partial(0, 2));
    let n = su.cross(sv);
    let nn = n.norm_sqr();
    if nn.lo.is_nan() || nn.lo <= 0.0 {
        return None;
    }
    let (l, m, nu) = (suu.dot(n), suv.dot(n), svv.dot(n));
    let gaussian = (l * nu - m.sqr()) / nn.sqr();
    let (e, f, g) = (su.norm_sqr(), su.dot(sv), sv.norm_sqr());
    let mean = (g * l - f.scale(2.0) * m + e * nu) / (nn * nn.sqrt()).scale(2.0);
    (gaussian.is_finite() && mean.is_finite()).then_some((gaussian, mean))
}

/// Guaranteed curvature bounds over `u_range × v_range` by adaptive
/// subdivision. A cell is subdivided while its Gaussian enclosure is
/// wider than `target_width`, or while it cannot be enclosed, up to
/// `max_subdivisions` levels.
pub fn guaranteed_curvature_bounds(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    target_width: f64,
    max_subdivisions: usize,
) -> GuaranteedCurvatureBounds {
    let unbounded = |subdivisions: usize, degenerate_cells: usize| GuaranteedCurvatureBounds {
        gaussian: Interval::new(f64::NEG_INFINITY, f64::INFINITY),
        mean: Interval::new(f64::NEG_INFINITY, f64::INFINITY),
        gaussian_width: f64::INFINITY,
        subdivisions,
        degenerate_cells,
    };
    let Ok(mut source) = SurfaceJetSource::new(surface) else {
        return unbounded(0, 1);
    };
    let mut cells = vec![(
        Cell {
            u: u_range,
            v: v_range,
        },
        0usize,
    )];
    let mut hull: Option<(Iv, Iv)> = None;
    let (mut max_depth, mut degenerate_cells) = (0usize, 0usize);
    while let Some((cell, depth)) = cells.pop() {
        max_depth = max_depth.max(depth);
        let bounds = cell_bounds(&mut source, &cell);
        let refine = depth < max_subdivisions
            && bounds.is_none_or(|(gaussian, _)| gaussian.width() > target_width);
        if refine {
            cells.extend(cell.subdivide().into_iter().map(|c| (c, depth + 1)));
            continue;
        }
        match bounds {
            Some((gaussian, mean)) => {
                hull = Some(match hull {
                    Some((g, h)) => (g.hull(gaussian), h.hull(mean)),
                    None => (gaussian, mean),
                });
            }
            None => degenerate_cells += 1,
        }
    }
    match hull {
        Some((gaussian, mean)) if degenerate_cells == 0 => GuaranteedCurvatureBounds {
            gaussian: Interval::new(gaussian.lo, gaussian.hi),
            mean: Interval::new(mean.lo, mean.hi),
            gaussian_width: gaussian.width(),
            subdivisions: max_depth,
            degenerate_cells: 0,
        },
        _ => unbounded(max_depth, degenerate_cells.max(1)),
    }
}

#[cfg(test)]
#[path = "interval_curvature_tests.rs"]
mod tests;
