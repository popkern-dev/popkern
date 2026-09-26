//! Interchange conditioning for Gregory surfaces (PK-524, PK-551).
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
//! The fit is measured, not assumed, on a doubling grid:
//! - position: the worst deviation at every cell midpoint must be
//!   within the export's conditioning tolerance (a tenth of the
//!   representation's length uncertainty);
//! - curvature (PK-551): the 3D normal curvatures along and across
//!   each boundary, sampled densely and clustered toward the corners,
//!   must match the patch within [`CURVATURE_FIDELITY`] of its
//!   boundary curvature scale. Cap seams certify G2 across these
//!   boundaries; a position-only fit kept positions to 1e-9 but moved
//!   the cube-corner cap's transverse curvature by 0.115 at a corner
//!   (2.4 %), a round-trip G2 -> G1 drop on every corner-cell seam.
//!
//! The nodes are Chebyshev-Lobatto (clustered at the boundaries, where
//! the corner cells vary fastest): measured on that cap, 33 clustered
//! nodes hold the boundary curvature to 2.6e-4 where 33 uniform nodes
//! missed by 0.115. A patch whose fits never reach the curvature bar
//! keeps its most faithful fit up to 65 x 65: the pyramid and saddle
//! caps' ε-regularized corner cells (a near-base-point 2x2 weight block
//! of 1e-14 at one corner) pollute interpolation near that corner, and
//! no tensor grid resolves them at a writable size (PK-556). One whose
//! fits never reach the position bar keeps the exact form.
//! Well-conditioned patches never enter this path.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsSurface, Point3, SurfaceGeometry, Vec3};

/// Weight spread below which a rational form counts as
/// ill-conditioned for interchange: `min / max` of its weights.
const MIN_INTERCHANGE_WEIGHT_RATIO: f64 = 1.0e-8;

/// Interpolation grid sizes tried, coarse to fine. The cap bounds the
/// written patch at 65 x 65 control points: 129 x 129 held the pyramid
/// cap's corner layer to 2.4e-2 but wrote 9 MB per cap and slowed
/// OpenCASCADE's self-interference check more than fourfold.
const GRID_COUNTS: [usize; 3] = [17, 33, 65];

/// Degree of the conditioned fit in both directions.
const FIT_DEGREE: usize = 5;

/// Allowed boundary normal-curvature error, relative to the largest
/// boundary curvature of the patch (with a floor of the inverse patch
/// size for nearly flat patches). At the measured cap scales (|k| up to
/// ~13) this holds the error to ~1e-3, a few percent of the commercial
/// G2 gate.
const CURVATURE_FIDELITY: f64 = 1.0e-4;

/// Uniform boundary samples, plus log-spaced samples toward each end
/// (from 1e-6 of the boundary) where the corner cells live.
const BOUNDARY_SAMPLES: usize = 200;
const CORNER_SAMPLES: usize = 40;

/// The B-spline to write for a Gregory `surface`: its exact rational
/// conversion, or — when that conversion's weights are too spread for
/// consumers — the first fit on the doubling grid whose measured
/// midpoint deviation is within `tolerance` and whose boundary
/// curvature matches the patch (see the module docs).
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
    conditioned_fit(surface, (u0, u1), (v0, v1), tolerance).unwrap_or(exact)
}

/// The conditioned non-rational fit of `surface` over its chart
/// `[u0, u1] x [v0, v1]`: the first clustered grid meeting both the
/// position and the curvature bar, else the most curvature-faithful
/// position-valid fit, else `None`.
fn conditioned_fit(
    surface: &SurfaceGeometry,
    (u0, u1): (f64, f64),
    (v0, v1): (f64, f64),
    tolerance: f64,
) -> Option<NurbsSurface> {
    let boundary = boundary_samples(surface, (u0, u1), (v0, v1));
    let curvature_tolerance = CURVATURE_FIDELITY * boundary.curvature_scale;
    let mut best: Option<(NurbsSurface, f64)> = None;
    for &count in &GRID_COUNTS {
        let (us, vs) = (clustered(u0, u1, count), clustered(v0, v1, count));
        let grid: Vec<Vec<Point3>> = vs
            .iter()
            .map(|&v| us.iter().map(|&u| surface.point_at(u, v)).collect())
            .collect();
        let Ok(fit) = popkern_core::geometry::fitting::interpolate_surface_grid(
            &grid, &us, &vs, FIT_DEGREE, FIT_DEGREE,
        ) else {
            continue;
        };
        if midpoint_deviation(&fit, surface, &us, &vs) > tolerance {
            continue;
        }
        let curvature_error = boundary.curvature_error(&fit);
        if curvature_error <= curvature_tolerance {
            return Some(fit);
        }
        if best.as_ref().is_none_or(|(_, e)| curvature_error < *e) {
            best = Some((fit, curvature_error));
        }
    }
    best.map(|(fit, _)| fit)
}

/// `count` Chebyshev-Lobatto nodes on `[a, b]` (endpoints included).
fn clustered(a: f64, b: f64, count: usize) -> Vec<f64> {
    (0..count)
        .map(|k| {
            let x = k as f64 / (count - 1) as f64;
            a + (b - a) * (0.5 - 0.5 * (std::f64::consts::PI * x).cos())
        })
        .collect()
}

/// Worst fit deviation at the midpoints of the interpolation cells.
fn midpoint_deviation(
    fit: &NurbsSurface,
    surface: &SurfaceGeometry,
    us: &[f64],
    vs: &[f64],
) -> f64 {
    us.windows(2)
        .flat_map(|du| vs.windows(2).map(move |dv| (du, dv)))
        .map(|(du, dv)| {
            let (u, v) = (0.5 * (du[0] + du[1]), 0.5 * (dv[0] + dv[1]));
            fit.point_at(u, v).distance_to(&surface.point_at(u, v))
        })
        .fold(0.0_f64, f64::max)
}

/// The patch's boundary: sample `(u, v)` with the boundary direction,
/// its normal curvatures there, and their scale.
struct BoundarySamples {
    points: Vec<((f64, f64), bool, [f64; 2])>,
    curvature_scale: f64,
}

impl BoundarySamples {
    /// Largest normal-curvature difference between `fit` and the patch
    /// over the samples (non-finite readings count as infinite).
    fn curvature_error(&self, fit: &NurbsSurface) -> f64 {
        self.points
            .iter()
            .map(
                |&((u, v), along_u, patch)| match boundary_curvatures(fit, u, v, along_u) {
                    Some(k) => (k[0] - patch[0]).abs().max((k[1] - patch[1]).abs()),
                    None => f64::INFINITY,
                },
            )
            .fold(0.0_f64, f64::max)
    }
}

fn boundary_samples(
    surface: &SurfaceGeometry,
    (u0, u1): (f64, f64),
    (v0, v1): (f64, f64),
) -> BoundarySamples {
    let fractions: Vec<f64> = (0..=BOUNDARY_SAMPLES)
        .map(|k| k as f64 / BOUNDARY_SAMPLES as f64)
        .chain((0..CORNER_SAMPLES).flat_map(|k| {
            let d = 10f64.powf(-6.0 + 5.0 * k as f64 / CORNER_SAMPLES as f64);
            [d, 1.0 - d]
        }))
        .collect();
    let mut points = Vec::with_capacity(4 * fractions.len());
    let mut scale = 0.0_f64;
    for &f in &fractions {
        let (u, v) = (u0 + f * (u1 - u0), v0 + f * (v1 - v0));
        for (at, along_u) in [
            ((u, v0), true),
            ((u, v1), true),
            ((u0, v), false),
            ((u1, v), false),
        ] {
            if let Some(k) = boundary_curvatures(surface, at.0, at.1, along_u) {
                scale = scale.max(k[0].abs()).max(k[1].abs());
                points.push((at, along_u, k));
            }
        }
    }
    let size = surface
        .point_at(u0, v0)
        .distance_to(&surface.point_at(u1, v1))
        .max(
            surface
                .point_at(u1, v0)
                .distance_to(&surface.point_at(u0, v1)),
        );
    let floor = if size > 0.0 { 1.0 / size } else { 1.0 };
    BoundarySamples {
        points,
        curvature_scale: scale.max(floor),
    }
}

/// Normal curvatures at `(u, v)` along the boundary direction (`u` when
/// `along_u`, else `v`) and across it, in 3D: `II(x, x) / I(x, x)`.
fn boundary_curvatures(s: &dyn SurfaceEval, u: f64, v: f64, along_u: bool) -> Option<[f64; 2]> {
    let (su, sv) = (s.partial_u(u, v), s.partial_v(u, v));
    let n = su.cross(&sv).normalized()?;
    let along = if along_u { su } else { sv }.normalized()?;
    let across = n.cross(&along);
    let (e, f, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
    let det = e * g - f * f;
    if det.is_nan() || det <= 0.0 {
        return None;
    }
    let (l, m, nn) = (
        s.second_partial_uu(u, v).dot(&n),
        s.second_partial_uv(u, v).dot(&n),
        s.second_partial_vv(u, v).dot(&n),
    );
    let curvature = |w: Vec3| {
        let (r1, r2) = (w.dot(&su), w.dot(&sv));
        let (a, b) = ((r1 * g - r2 * f) / det, (r2 * e - r1 * f) / det);
        (l * a * a + 2.0 * m * a * b + nn * b * b) / (e * a * a + 2.0 * f * a * b + g * b * b)
    };
    let k = [curvature(along), curvature(across)];
    k.iter().all(|x| x.is_finite()).then_some(k)
}

#[cfg(test)]
#[path = "export_surface_conditioning_tests.rs"]
mod tests;
