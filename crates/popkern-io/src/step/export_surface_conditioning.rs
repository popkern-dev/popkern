//! Interchange conditioning for Gregory surfaces (PK-524).
//!
//! STEP has no Gregory entity, so a Gregory patch is written as the
//! B-spline its `to_nurbs()` conversion returns. That conversion is
//! EXACT, but the dual cross-jet corner cells (PK-519) blend their
//! per-boundary claims with squared quadrant weights that vanish at
//! the cell corner; the kernel ε-regularizes those slots, so the exact
//! form is a rational bi-11 surface whose weights span 1e-14 .. 2 — a
//! near-base-point at the corner. popkern evaluates it exactly, but
//! consumers do not: OpenCASCADE's STEP import inflated the
//! tolerances of the triangular-pyramid cap's rim-quarter edges to
//! 5.5e-2 (square: 1.1e-2) while every curve, pcurve and vertex was
//! exact on the popkern side, and its boolean pre-check then read
//! neighbouring cap corners as self-interfering.
//!
//! An ill-conditioned exact form is therefore replaced by a
//! non-rational bi-quintic interpolation of the patch over the SAME
//! `(u, v)` chart, so every pcurve the exporter writes stays valid.
//! The fit is measured, not assumed: the patch is sampled on a
//! doubling grid and the fit is checked against the exact surface at
//! every cell midpoint until the worst deviation is within the
//! export's conditioning tolerance (a tenth of the representation's
//! length uncertainty). A fit that cannot get there keeps the exact
//! form, and well-conditioned patches never enter this path.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsSurface, Point3, SurfaceGeometry};

/// Weight spread below which a rational form counts as
/// ill-conditioned for interchange: `min / max` of its weights.
const MIN_INTERCHANGE_WEIGHT_RATIO: f64 = 1.0e-8;

/// Interpolation grid sizes tried, coarse to fine (`2^k + 1`).
const GRID_COUNTS: [usize; 4] = [17, 33, 65, 129];

/// Degree of the conditioned fit in both directions.
const FIT_DEGREE: usize = 5;

/// The B-spline to write for a Gregory `surface`: its exact rational
/// conversion, or — when that conversion's weights are too spread for
/// consumers — the first non-rational fit on the doubling grid whose
/// measured midpoint deviation is within `tolerance`.
pub(super) fn conditioned_gregory_nurbs(surface: &SurfaceGeometry, tolerance: f64) -> NurbsSurface {
    let exact = surface.to_nurbs();
    let (min_w, max_w) = exact
        .weights
        .iter()
        .fold((f64::INFINITY, 0.0_f64), |(lo, hi), &w| {
            (lo.min(w), hi.max(w))
        });
    if max_w <= 0.0 || max_w.is_nan() || min_w / max_w >= MIN_INTERCHANGE_WEIGHT_RATIO {
        return exact;
    }
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    if !(u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite()) {
        return exact;
    }
    for &count in &GRID_COUNTS {
        let spread = |a: f64, b: f64| -> Vec<f64> {
            (0..count)
                .map(|k| a + (b - a) * k as f64 / (count - 1) as f64)
                .collect()
        };
        let (us, vs) = (spread(u0, u1), spread(v0, v1));
        let grid: Vec<Vec<Point3>> = vs
            .iter()
            .map(|&v| us.iter().map(|&u| surface.point_at(u, v)).collect())
            .collect();
        let Ok(fit) = popkern_core::geometry::fitting::interpolate_surface_grid(
            &grid, &us, &vs, FIT_DEGREE, FIT_DEGREE,
        ) else {
            continue;
        };
        let deviation = us
            .windows(2)
            .flat_map(|du| vs.windows(2).map(move |dv| (du, dv)))
            .map(|(du, dv)| {
                let (u, v) = (0.5 * (du[0] + du[1]), 0.5 * (dv[0] + dv[1]));
                fit.point_at(u, v).distance_to(&surface.point_at(u, v))
            })
            .fold(0.0_f64, f64::max);
        if deviation <= tolerance {
            return fit;
        }
    }
    exact
}
