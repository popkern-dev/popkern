//! Visual quality analysis: isophotes, reflection lines, highlight lines,
//! and zebra stripes computed directly on parametric surfaces.
//!
//! These tools reveal continuity defects that are invisible in standard
//! rendering but visible in reflective/glossy surface display.

use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};

// ---------------------------------------------------------------------------
// Isophotes
// ---------------------------------------------------------------------------

/// A sampled isophote curve on a surface.
///
/// An isophote is a contour of constant light intensity `N . L`,
/// where `N` is the surface normal and `L` is the light direction.
/// Kinks in isophotes across patch boundaries reveal G1 discontinuities.
#[derive(Debug, Clone)]
pub struct Isophote {
    /// Light intensity level for this isophote.
    pub level: f64,
    /// Sampled points along the isophote (u, v, point3d).
    pub samples: Vec<(f64, f64, Point3)>,
}

/// Compute isophotes on a surface for a given light direction.
///
/// Samples the surface on a grid and traces contours of constant
/// `N . light_dir` at the specified intensity levels. Each contour
/// is extracted using marching-squares style linear interpolation
/// along grid cell edges.
///
/// # Parameters
///
/// - `surface`: parametric surface to analyse
/// - `u_range`, `v_range`: parameter domain to sample
/// - `light_dir`: directional light (need not be unit length)
/// - `levels`: intensity values at which to trace isophotes
/// - `grid_density`: number of sample points in each parameter direction
pub fn compute_isophotes(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
    light_dir: Vec3,
    levels: &[f64],
    grid_density: usize,
) -> Vec<Isophote> {
    let grid_density = grid_density.max(2);
    let field = sample_scalar_field(surface, u_range, v_range, grid_density, |n| {
        n.dot(&light_dir)
    });
    extract_contours(surface, u_range, v_range, &field, grid_density, levels)
}

// ---------------------------------------------------------------------------
// Reflection lines
// ---------------------------------------------------------------------------

/// A sampled reflection line on a surface.
///
/// Reflection lines simulate the reflection of parallel light strips
/// on a glossy surface. Kinks in these lines reveal G2 breaks between
/// surface patches.
#[derive(Debug, Clone)]
pub struct ReflectionLine {
    /// Stripe index.
    pub index: i32,
    /// Sampled points along the reflection line.
    pub samples: Vec<(f64, f64, Point3)>,
}

/// Compute reflection lines on a surface.
///
/// Simulates the reflection of parallel light strips viewed from
/// `view_dir`. For each grid point, the reflected view direction is
/// computed as `R = 2*(N.V)*N - V` and then projected onto
/// `stripe_dir`. Contours are traced at integer multiples of
/// `stripe_spacing`, each representing one reflected stripe.
///
/// # Parameters
///
/// - `surface`: parametric surface to analyse
/// - `u_range`, `v_range`: parameter domain to sample
/// - `view_dir`: viewer direction (need not be unit length)
/// - `stripe_dir`: direction of the light strips
/// - `stripe_spacing`: distance between adjacent stripes
/// - `grid_density`: number of sample points in each parameter direction
pub fn compute_reflection_lines(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
    view_dir: Vec3,
    stripe_dir: Vec3,
    stripe_spacing: f64,
    grid_density: usize,
) -> Vec<ReflectionLine> {
    let grid_density = grid_density.max(2);
    let spacing = if stripe_spacing.abs() < 1e-30 {
        1.0
    } else {
        stripe_spacing
    };

    // Sample the stripe coordinate field.
    let field = sample_scalar_field(surface, u_range, v_range, grid_density, |n| {
        let n_dot_v = n.dot(&view_dir);
        let reflected = n * (2.0 * n_dot_v) - view_dir;
        reflected.dot(&stripe_dir) / spacing
    });

    // Determine the range of integer stripe indices in the field.
    let (min_val, max_val) = field_range(&field);
    let i_min = min_val.floor() as i32;
    let i_max = max_val.ceil() as i32;

    let mut lines = Vec::new();
    for idx in i_min..=i_max {
        let level = idx as f64;
        let contours = extract_contours(surface, u_range, v_range, &field, grid_density, &[level]);
        for iso in contours {
            lines.push(ReflectionLine {
                index: idx,
                samples: iso.samples,
            });
        }
    }
    lines
}

// ---------------------------------------------------------------------------
// Zebra stripes
// ---------------------------------------------------------------------------

/// Zebra stripe analysis result.
///
/// A grid of binary stripe values over the parameter domain.
/// Discontinuities in the stripe pattern across patch boundaries
/// reveal G1 breaks.
#[derive(Debug, Clone)]
pub struct ZebraStripeField {
    /// Grid of stripe values: 0.0 (black) or 1.0 (white).
    pub grid: Vec<Vec<f64>>,
    /// Parameter grid: (u, v) for each cell.
    pub params: Vec<Vec<(f64, f64)>>,
    /// Grid dimensions (u_count, v_count).
    pub dimensions: (usize, usize),
}

/// Compute zebra stripe field on a surface.
///
/// Alternating black/white reflected strips. The stripe coordinate
/// is computed identically to [`compute_reflection_lines`], then
/// quantised: `floor(s) % 2 == 0` yields white (1.0), otherwise
/// black (0.0). Discontinuities in the stripe pattern reveal G1
/// breaks.
///
/// # Parameters
///
/// - `surface`: parametric surface to analyse
/// - `u_range`, `v_range`: parameter domain to sample
/// - `view_dir`: viewer direction
/// - `stripe_dir`: direction of the light strips
/// - `stripe_spacing`: distance between adjacent stripes
/// - `grid_density`: number of sample points in each parameter direction
pub fn compute_zebra_stripes(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
    view_dir: Vec3,
    stripe_dir: Vec3,
    stripe_spacing: f64,
    grid_density: usize,
) -> ZebraStripeField {
    let grid_density = grid_density.max(2);
    let spacing = if stripe_spacing.abs() < 1e-30 {
        1.0
    } else {
        stripe_spacing
    };

    let du = (u_range.1 - u_range.0) / (grid_density as f64 - 1.0);
    let dv = (v_range.1 - v_range.0) / (grid_density as f64 - 1.0);

    let mut grid = Vec::with_capacity(grid_density);
    let mut params = Vec::with_capacity(grid_density);

    for i in 0..grid_density {
        let u = u_range.0 + du * i as f64;
        let mut row = Vec::with_capacity(grid_density);
        let mut param_row = Vec::with_capacity(grid_density);

        for j in 0..grid_density {
            let v = v_range.0 + dv * j as f64;
            let normal = surface.normal_at(u, v);
            let n = normal.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
            let n_dot_v = n.dot(&view_dir);
            let reflected = n * (2.0 * n_dot_v) - view_dir;
            let s = reflected.dot(&stripe_dir) / spacing;
            let stripe_val = if (s.floor() as i64) % 2 == 0 {
                1.0
            } else {
                0.0
            };
            row.push(stripe_val);
            param_row.push((u, v));
        }
        grid.push(row);
        params.push(param_row);
    }

    ZebraStripeField {
        grid,
        params,
        dimensions: (grid_density, grid_density),
    }
}

// ---------------------------------------------------------------------------
// Highlight lines
// ---------------------------------------------------------------------------

/// Compute highlight lines from a line light source.
///
/// Highlight lines are loci of specular highlights: the contours
/// where `N . H` is near 1.0, with `H = normalize(L + V)` being
/// the half-vector between the light and view directions. These
/// contours are extracted at a fixed set of levels near the specular
/// peak. They are returned as [`Isophote`] values because the
/// structure is identical to isophotes.
///
/// # Parameters
///
/// - `surface`: parametric surface to analyse
/// - `u_range`, `v_range`: parameter domain to sample
/// - `light_dir`: light direction
/// - `view_dir`: view direction
/// - `grid_density`: number of sample points in each parameter direction
pub fn compute_highlight_lines(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
    light_dir: Vec3,
    view_dir: Vec3,
    grid_density: usize,
) -> Vec<Isophote> {
    let grid_density = grid_density.max(2);
    let half_vec = (light_dir + view_dir)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 0.0, 1.0));

    let field = sample_scalar_field(surface, u_range, v_range, grid_density, |n| {
        n.dot(&half_vec)
    });

    // Highlight lines are near-specular: trace contours at high N.H values.
    let levels: Vec<f64> = (0..10).map(|i| 0.91 + 0.01 * i as f64).collect();
    extract_contours(surface, u_range, v_range, &field, grid_density, &levels)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Sample a scalar field over the parameter domain.
///
/// Evaluates the surface normal at each grid point, normalises it,
/// and feeds it to `scalar_fn` which returns the field value.
fn sample_scalar_field(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
    grid_density: usize,
    scalar_fn: impl Fn(Vec3) -> f64,
) -> Vec<Vec<f64>> {
    let du = (u_range.1 - u_range.0) / (grid_density as f64 - 1.0);
    let dv = (v_range.1 - v_range.0) / (grid_density as f64 - 1.0);

    let mut field = Vec::with_capacity(grid_density);
    for i in 0..grid_density {
        let u = u_range.0 + du * i as f64;
        let mut row = Vec::with_capacity(grid_density);
        for j in 0..grid_density {
            let v = v_range.0 + dv * j as f64;
            let normal = surface.normal_at(u, v);
            let n = normal.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
            row.push(scalar_fn(n));
        }
        field.push(row);
    }
    field
}

/// Find the (min, max) of a 2D scalar field.
fn field_range(field: &[Vec<f64>]) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for row in field {
        for &v in row {
            if v < lo {
                lo = v;
            }
            if v > hi {
                hi = v;
            }
        }
    }
    if lo > hi { (0.0, 0.0) } else { (lo, hi) }
}

/// Extract iso-contours from a sampled scalar field using marching squares.
///
/// For each requested level, walks every grid cell. If the level crosses
/// any cell edge (i.e. the values at the two endpoints straddle the level),
/// the crossing point is found by linear interpolation and recorded.
fn extract_contours(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
    field: &[Vec<f64>],
    grid_density: usize,
    levels: &[f64],
) -> Vec<Isophote> {
    let du = (u_range.1 - u_range.0) / (grid_density as f64 - 1.0);
    let dv = (v_range.1 - v_range.0) / (grid_density as f64 - 1.0);

    let mut result = Vec::with_capacity(levels.len());

    for &level in levels {
        let mut samples: Vec<(f64, f64, Point3)> = Vec::new();

        for i in 0..grid_density - 1 {
            for j in 0..grid_density - 1 {
                // Cell corners: (i,j), (i+1,j), (i+1,j+1), (i,j+1)
                let v00 = field[i][j];
                let v10 = field[i + 1][j];
                let v11 = field[i + 1][j + 1];
                let v01 = field[i][j + 1];

                let u0 = u_range.0 + du * i as f64;
                let u1 = u_range.0 + du * (i + 1) as f64;
                let v0 = v_range.0 + dv * j as f64;
                let v1 = v_range.0 + dv * (j + 1) as f64;

                // Check each of the 4 cell edges for a crossing.
                // Bottom edge: (i,j) -> (i+1,j)
                if let Some((u, v)) = interpolate_edge(v00, v10, u0, u1, v0, v0, level) {
                    samples.push((u, v, surface.point_at(u, v)));
                }
                // Right edge: (i+1,j) -> (i+1,j+1)
                if let Some((u, v)) = interpolate_edge(v10, v11, u1, u1, v0, v1, level) {
                    samples.push((u, v, surface.point_at(u, v)));
                }
                // Top edge: (i,j+1) -> (i+1,j+1)
                if let Some((u, v)) = interpolate_edge(v01, v11, u0, u1, v1, v1, level) {
                    samples.push((u, v, surface.point_at(u, v)));
                }
                // Left edge: (i,j) -> (i,j+1)
                if let Some((u, v)) = interpolate_edge(v00, v01, u0, u0, v0, v1, level) {
                    samples.push((u, v, surface.point_at(u, v)));
                }
            }
        }

        result.push(Isophote { level, samples });
    }
    result
}

/// Linear interpolation along a cell edge to find the crossing point.
///
/// Returns `Some((u, v))` if the level lies strictly between `val_a`
/// and `val_b`, otherwise `None`.
fn interpolate_edge(
    val_a: f64,
    val_b: f64,
    u_a: f64,
    u_b: f64,
    v_a: f64,
    v_b: f64,
    level: f64,
) -> Option<(f64, f64)> {
    let min_val = val_a.min(val_b);
    let max_val = val_a.max(val_b);
    if level < min_val || level > max_val {
        return None;
    }
    let denom = val_b - val_a;
    if denom.abs() < 1e-30 {
        return None;
    }
    let t = (level - val_a) / denom;
    let u = u_a + t * (u_b - u_a);
    let v = v_a + t * (v_b - v_a);
    Some((u, v))
}

#[cfg(test)]
#[path = "visual_analysis_tests.rs"]
mod tests;
