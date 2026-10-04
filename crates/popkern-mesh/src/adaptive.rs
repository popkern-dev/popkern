//! Chord-height adaptive UV-grid subdivision for parametric surfaces.

use popkern_core::geometry::{Point3, SurfaceEval, SurfaceGeometry};

/// A rectangular cell in UV parameter space.
#[derive(Debug, Clone)]
pub(crate) struct UvCell {
    /// Minimum u parameter of this cell.
    pub u_min: f64,
    /// Maximum u parameter of this cell.
    pub u_max: f64,
    /// Minimum v parameter of this cell.
    pub v_min: f64,
    /// Maximum v parameter of this cell.
    pub v_max: f64,
}

/// The per-face subdivision cap, scaled to the surface's parametric
/// extent: a multi-turn swept NURBS (a 6.5-turn spring is 26 knot
/// spans in v) needs ~17 samples per coil turn to hold the chord
/// tolerance, far past a fixed cap of 64 — the capped mesh rendered
/// visibly faceted coils (PK-440). Eight cells per knot span keeps
/// single-span faces at the configured cap while long periodic sweeps
/// refine to their true chord tolerance; 512 bounds the quad-tree.
pub(crate) fn effective_max_subdivisions(
    surface: &SurfaceGeometry,
    options: &crate::TessellationOptions,
) -> usize {
    let span_count = |knots: &[f64]| knots.windows(2).filter(|w| w[1] > w[0] + 1e-12).count();
    match surface {
        SurfaceGeometry::Nurbs(ns) => options
            .max_subdivisions
            .max(8 * span_count(&ns.knots_u).max(span_count(&ns.knots_v)))
            .min(512),
        _ => options.max_subdivisions,
    }
}

/// Stations per chart direction at which the worst iso-line chord
/// sagitta (three cross-stations per direction) meets `chord_height`:
/// the larger direction's count, doubled from `start` and capped at
/// `cap`. Sizes the quadtree cap of whole periodic charts (a sphere or
/// torus face bounded only by its seam), whose fixed cap of 64 left an
/// r 254 sphere's facets 0.34 off at a 0.005 chord (PK-654).
pub(crate) fn chord_station_count(
    surface: &SurfaceGeometry,
    (u0, u1): (f64, f64),
    (v0, v1): (f64, f64),
    chord_height: f64,
    (start, cap): (usize, usize),
) -> usize {
    let sag_at = |n: usize, along_u: bool| -> f64 {
        let mut worst = 0.0_f64;
        for cross in [0.25, 0.5, 0.75] {
            let at = |f: f64| {
                if along_u {
                    surface.point_at(u0 + (u1 - u0) * f, v0 + (v1 - v0) * cross)
                } else {
                    surface.point_at(u0 + (u1 - u0) * cross, v0 + (v1 - v0) * f)
                }
            };
            for k in 0..n {
                let (a, b) = (at(k as f64 / n as f64), at((k + 1) as f64 / n as f64));
                let mid = at((k as f64 + 0.5) / n as f64);
                let chord_mid =
                    Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
                worst = worst.max(mid.distance_to(&chord_mid));
            }
        }
        worst
    };
    let mut n = start.max(1);
    while n < cap && (sag_at(n, true) > chord_height || sag_at(n, false) > chord_height) {
        n *= 2;
    }
    n.min(cap)
}

/// Compute the chord-height error for a UV cell.
///
/// This measures the distance from the bilinear interpolation of the four
/// corner points to the actual surface point at five interior samples,
/// and from the emitted facets (the shorter-diagonal split) at the cell's
/// UV midpoint.
/// A flat (zero-curvature) surface will return approximately zero; a curved
/// surface will return a positive value proportional to its local curvature.
pub(crate) fn chord_height(surface: &SurfaceGeometry, cell: &UvCell) -> f64 {
    let p00 = surface.point_at(cell.u_min, cell.v_min);
    let p10 = surface.point_at(cell.u_max, cell.v_min);
    let p01 = surface.point_at(cell.u_min, cell.v_max);
    let p11 = surface.point_at(cell.u_max, cell.v_max);
    let samples = [
        (0.5, 0.5),
        (0.25, 0.25),
        (0.75, 0.25),
        (0.25, 0.75),
        (0.75, 0.75),
    ];
    let bilinear = samples
        .iter()
        .map(|&(su, sv)| {
            let u = cell.u_min + su * (cell.u_max - cell.u_min);
            let v = cell.v_min + sv * (cell.v_max - cell.v_min);
            let interp = bilinear_point(p00, p10, p01, p11, su, sv);
            interp.distance_to(&surface.point_at(u, v))
        })
        .fold(0.0, f64::max);
    // The cell is emitted as two triangles split along its shorter
    // diagonal, whose midpoint is the facet point over the cell centre.
    // A bilinear patch reproduces a saddle exactly, so the bilinear
    // samples read a torus's inner equator as flat while the facets sat
    // 0.017 off at a 0.005 chord (PK-654).
    let centre = surface.point_at(
        0.5 * (cell.u_min + cell.u_max),
        0.5 * (cell.v_min + cell.v_max),
    );
    let mid =
        |a: Point3, b: Point3| Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
    let diagonal = if p00.distance_to(&p11) <= p10.distance_to(&p01) {
        mid(p00, p11)
    } else {
        mid(p10, p01)
    };
    // Each cell side is an emitted edge: its own chord sag (a saddle's
    // opposite-signed side sags cancel in every interior sample).
    let (um, vm) = (
        0.5 * (cell.u_min + cell.u_max),
        0.5 * (cell.v_min + cell.v_max),
    );
    let sides = [
        (p00, p10, (um, cell.v_min)),
        (p01, p11, (um, cell.v_max)),
        (p00, p01, (cell.u_min, vm)),
        (p10, p11, (cell.u_max, vm)),
    ]
    .iter()
    .map(|&(a, b, (u, v))| surface.point_at(u, v).distance_to(&mid(a, b)))
    .fold(0.0, f64::max);
    bilinear.max(centre.distance_to(&diagonal)).max(sides)
}

fn bilinear_point(p00: Point3, p10: Point3, p01: Point3, p11: Point3, u: f64, v: f64) -> Point3 {
    let w00 = (1.0 - u) * (1.0 - v);
    let w10 = u * (1.0 - v);
    let w01 = (1.0 - u) * v;
    let w11 = u * v;
    Point3::new(
        p00.x * w00 + p10.x * w10 + p01.x * w01 + p11.x * w11,
        p00.y * w00 + p10.y * w10 + p01.y * w01 + p11.y * w11,
        p00.z * w00 + p10.z * w10 + p01.z * w01 + p11.z * w11,
    )
}

/// Adaptively subdivide a surface's UV parameter space into a grid of cells.
///
/// Starting from an initial `min_subdivisions × min_subdivisions` grid, each
/// cell whose chord-height exceeds `chord_tolerance` is refined into four
/// sub-cells (quad-tree subdivision). Refinement continues until all leaf
/// cells satisfy the tolerance or until the subdivision depth would exceed the
/// level implied by `max_subdivisions`.
///
/// # Parameters
/// - `surface` — the parametric surface to subdivide.
/// - `u_range` — `(u_min, u_max)` parameter extent in the u direction.
/// - `v_range` — `(v_min, v_max)` parameter extent in the v direction.
/// - `min_subdivisions` — initial grid resolution along each axis (e.g. `2`
///   creates a 2×2 = 4-cell starting grid).
/// - `max_subdivisions` — maximum resolution along each axis; no cell will be
///   subdivided further once the axis resolution would exceed this value.
/// - `chord_tolerance` — maximum allowed chord-height deviation.
///
/// # Returns
/// A `Vec<UvCell>` of all leaf cells covering the parameter domain.
#[cfg(test)]
pub(crate) fn adaptive_subdivide(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    min_subdivisions: usize,
    max_subdivisions: usize,
    chord_tolerance: f64,
) -> Vec<UvCell> {
    adaptive_subdivide_bounded(
        surface,
        u_range,
        v_range,
        (min_subdivisions, max_subdivisions),
        chord_tolerance,
        None,
    )
}

/// [`adaptive_subdivide`] with an optional angular facet bound: a
/// cell also refines while the surface normal spreads more than
/// `max_normal_spread` radians across it (see
/// [`crate::facet_angle::normal_spread`]), within the same depth cap.
pub(crate) fn adaptive_subdivide_bounded(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    (min_subdivisions, max_subdivisions): (usize, usize),
    chord_tolerance: f64,
    max_normal_spread: Option<f64>,
) -> Vec<UvCell> {
    // Compute max refinement depth: how many times we can double min_subdivisions
    // before reaching max_subdivisions.
    //   depth 0: min_subdivisions cells/axis
    //   depth 1: 2*min_subdivisions cells/axis
    //   depth d: 2^d * min_subdivisions cells/axis
    // We want 2^d * min_subdivisions <= max_subdivisions
    //   => 2^d <= max_subdivisions / min_subdivisions
    //   => d <= log2(max_subdivisions / min_subdivisions)
    let max_depth = if min_subdivisions == 0 || max_subdivisions <= min_subdivisions {
        0_u32
    } else {
        let ratio = max_subdivisions / min_subdivisions;
        // ilog2 gives floor(log2(ratio))
        (ratio as u32).ilog2()
    };

    let (u_min, u_max) = u_range;
    let (v_min, v_max) = v_range;

    // Build initial grid.
    let mut queue: Vec<(UvCell, u32)> = Vec::new();
    let n = min_subdivisions.max(1);
    let du = (u_max - u_min) / n as f64;
    let dv = (v_max - v_min) / n as f64;
    for i in 0..n {
        for j in 0..n {
            queue.push((
                UvCell {
                    u_min: u_min + i as f64 * du,
                    u_max: u_min + (i + 1) as f64 * du,
                    v_min: v_min + j as f64 * dv,
                    v_max: v_min + (j + 1) as f64 * dv,
                },
                0,
            ));
        }
    }

    let mut leaves: Vec<UvCell> = Vec::new();

    while let Some((cell, depth)) = queue.pop() {
        let ch = chord_height(surface, &cell);
        let flat_enough = max_normal_spread.is_none_or(|m| {
            crate::facet_angle::normal_spread(
                surface,
                (cell.u_min, cell.u_max),
                (cell.v_min, cell.v_max),
            ) <= m
        });
        if (ch <= chord_tolerance && flat_enough) || depth >= max_depth {
            leaves.push(cell);
        } else {
            // Subdivide into 4 children.
            let u_mid = (cell.u_min + cell.u_max) * 0.5;
            let v_mid = (cell.v_min + cell.v_max) * 0.5;
            let next_depth = depth + 1;
            queue.push((
                UvCell {
                    u_min: cell.u_min,
                    u_max: u_mid,
                    v_min: cell.v_min,
                    v_max: v_mid,
                },
                next_depth,
            ));
            queue.push((
                UvCell {
                    u_min: u_mid,
                    u_max: cell.u_max,
                    v_min: cell.v_min,
                    v_max: v_mid,
                },
                next_depth,
            ));
            queue.push((
                UvCell {
                    u_min: cell.u_min,
                    u_max: u_mid,
                    v_min: v_mid,
                    v_max: cell.v_max,
                },
                next_depth,
            ));
            queue.push((
                UvCell {
                    u_min: u_mid,
                    u_max: cell.u_max,
                    v_min: v_mid,
                    v_max: cell.v_max,
                },
                next_depth,
            ));
        }
    }

    leaves
}

#[cfg(test)]
#[path = "adaptive_tests.rs"]
mod tests;
