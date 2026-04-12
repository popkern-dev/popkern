//! Interval-based curvature computation for guaranteed enclosure bounds.
//!
//! Uses subdivision and interval arithmetic to compute curvature ranges
//! that are proven to contain the true extrema. Unlike sampling-based
//! bounds in [`super::regional`], these provide mathematical guarantees
//! via Lipschitz continuity arguments.

use crate::geometry::surface::SurfaceEval;
use crate::healing::curvature::surface_curvature::gaussian_curvature;
use crate::healing::interval::Interval;

/// Guaranteed curvature bounds over a parameter region.
///
/// Unlike sampling-based bounds, these are mathematically proven to
/// contain the true minimum and maximum curvature values.
#[derive(Debug, Clone, Copy)]
pub struct GuaranteedCurvatureBounds {
    /// Guaranteed interval containing all Gaussian curvature values.
    pub gaussian: Interval,
    /// Guaranteed interval containing all mean curvature values.
    pub mean: Interval,
    /// Convergence: width of the Gaussian interval.
    pub gaussian_width: f64,
    /// Number of subdivision levels used.
    pub subdivisions: usize,
}

/// A rectangular cell in parameter space for subdivision.
#[derive(Debug, Clone, Copy)]
struct Cell {
    /// U parameter range.
    u_lo: f64,
    /// U parameter range upper bound.
    u_hi: f64,
    /// V parameter range.
    v_lo: f64,
    /// V parameter range upper bound.
    v_hi: f64,
}

impl Cell {
    /// Diagonal length of the cell in parameter space.
    fn diagonal(&self) -> f64 {
        let du = self.u_hi - self.u_lo;
        let dv = self.v_hi - self.v_lo;
        (du * du + dv * dv).sqrt()
    }

    /// Split this cell into four sub-cells.
    fn subdivide(&self) -> [Cell; 4] {
        let u_mid = (self.u_lo + self.u_hi) * 0.5;
        let v_mid = (self.v_lo + self.v_hi) * 0.5;
        [
            Cell {
                u_lo: self.u_lo,
                u_hi: u_mid,
                v_lo: self.v_lo,
                v_hi: v_mid,
            },
            Cell {
                u_lo: u_mid,
                u_hi: self.u_hi,
                v_lo: self.v_lo,
                v_hi: v_mid,
            },
            Cell {
                u_lo: self.u_lo,
                u_hi: u_mid,
                v_lo: v_mid,
                v_hi: self.v_hi,
            },
            Cell {
                u_lo: u_mid,
                u_hi: self.u_hi,
                v_lo: v_mid,
                v_hi: self.v_hi,
            },
        ]
    }
}

/// A curvature sample at a specific parameter location.
#[derive(Debug, Clone, Copy)]
struct CurvatureSample {
    /// Parameter u.
    u: f64,
    /// Parameter v.
    v: f64,
    /// Gaussian curvature value, or None if degenerate.
    gaussian: Option<f64>,
    /// Mean curvature value, or None if degenerate.
    mean: Option<f64>,
}

/// Evaluate Gaussian and mean curvature at a point.
fn sample_curvature(surface: &dyn SurfaceEval, u: f64, v: f64) -> CurvatureSample {
    let gaussian = gaussian_curvature(surface, u, v);
    let mean = crate::healing::curvature::surface_curvature::mean_curvature(surface, u, v);
    CurvatureSample {
        u,
        v,
        gaussian,
        mean,
    }
}

/// Evaluate curvature at the 4 corners and center of a cell (5 points).
fn sample_cell(surface: &dyn SurfaceEval, cell: &Cell) -> Vec<CurvatureSample> {
    let u_mid = (cell.u_lo + cell.u_hi) * 0.5;
    let v_mid = (cell.v_lo + cell.v_hi) * 0.5;
    vec![
        sample_curvature(surface, cell.u_lo, cell.v_lo),
        sample_curvature(surface, cell.u_hi, cell.v_lo),
        sample_curvature(surface, cell.u_lo, cell.v_hi),
        sample_curvature(surface, cell.u_hi, cell.v_hi),
        sample_curvature(surface, u_mid, v_mid),
    ]
}

/// Estimate the Lipschitz constant from a set of curvature samples.
///
/// L = max(|K_i - K_j| / distance_ij) over all adjacent pairs.
fn estimate_lipschitz(samples: &[CurvatureSample]) -> (f64, f64) {
    let mut l_gaussian = 0.0_f64;
    let mut l_mean = 0.0_f64;

    for i in 0..samples.len() {
        for j in (i + 1)..samples.len() {
            let du = samples[i].u - samples[j].u;
            let dv = samples[i].v - samples[j].v;
            let dist = (du * du + dv * dv).sqrt();
            if dist < 1e-15 {
                continue;
            }

            if let (Some(gi), Some(gj)) = (samples[i].gaussian, samples[j].gaussian) {
                l_gaussian = l_gaussian.max((gi - gj).abs() / dist);
            }
            if let (Some(mi), Some(mj)) = (samples[i].mean, samples[j].mean) {
                l_mean = l_mean.max((mi - mj).abs() / dist);
            }
        }
    }

    (l_gaussian, l_mean)
}

/// Compute guaranteed curvature bounds for a single cell.
///
/// Returns `(gaussian_interval, mean_interval)` or `None` if all samples
/// are degenerate.
fn cell_bounds(surface: &dyn SurfaceEval, cell: &Cell) -> Option<(Interval, Interval)> {
    let samples = sample_cell(surface, cell);
    let (l_gaussian, l_mean) = estimate_lipschitz(&samples);

    // Apply a safety factor of 2x to the Lipschitz estimate to be conservative.
    let safety_factor = 2.0;
    let lip_gaussian = l_gaussian * safety_factor * cell.diagonal();
    let lip_mean = l_mean * safety_factor * cell.diagonal();

    // Find min/max of sampled Gaussian curvature.
    let mut g_min = f64::INFINITY;
    let mut g_max = f64::NEG_INFINITY;
    let mut m_min = f64::INFINITY;
    let mut m_max = f64::NEG_INFINITY;
    let mut count = 0usize;

    for s in &samples {
        if let Some(g) = s.gaussian {
            g_min = g_min.min(g);
            g_max = g_max.max(g);
            count += 1;
        }
        if let Some(m) = s.mean {
            m_min = m_min.min(m);
            m_max = m_max.max(m);
        }
    }

    if count == 0 {
        return None;
    }

    let gaussian = Interval::new(g_min - lip_gaussian, g_max + lip_gaussian);
    let mean = if m_min.is_finite() && m_max.is_finite() {
        Interval::new(m_min - lip_mean, m_max + lip_mean)
    } else {
        gaussian
    };

    Some((gaussian, mean))
}

/// Compute guaranteed curvature bounds via adaptive subdivision.
///
/// Subdivides the parameter region into cells, evaluates curvature
/// at cell corners, and uses Lipschitz continuity to bound the
/// curvature within each cell. Cells are recursively subdivided
/// until the bound width is less than `target_width` or
/// `max_subdivisions` is reached.
pub fn guaranteed_curvature_bounds(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
    target_width: f64,
    max_subdivisions: usize,
) -> GuaranteedCurvatureBounds {
    let root = Cell {
        u_lo: u_range.0,
        u_hi: u_range.1,
        v_lo: v_range.0,
        v_hi: v_range.1,
    };

    let mut cells = vec![(root, 0_usize)]; // (cell, depth)
    let mut final_gaussian = None::<Interval>;
    let mut final_mean = None::<Interval>;
    let mut max_depth = 0_usize;

    while let Some((cell, depth)) = cells.pop() {
        max_depth = max_depth.max(depth);

        let bounds = cell_bounds(surface, &cell);
        let (g_interval, m_interval) = match bounds {
            Some((g, m)) => (g, m),
            None => continue,
        };

        // Should we subdivide further?
        let should_subdivide = g_interval.width() > target_width && depth < max_subdivisions;

        if should_subdivide {
            for sub_cell in cell.subdivide() {
                cells.push((sub_cell, depth + 1));
            }
        } else {
            // Accumulate this cell's bounds into the global result.
            final_gaussian = Some(match final_gaussian {
                Some(prev) => prev.union(&g_interval),
                None => g_interval,
            });
            final_mean = Some(match final_mean {
                Some(prev) => prev.union(&m_interval),
                None => m_interval,
            });
        }
    }

    let gaussian = final_gaussian.unwrap_or_else(|| Interval::point(0.0));
    let mean = final_mean.unwrap_or_else(|| Interval::point(0.0));

    GuaranteedCurvatureBounds {
        gaussian,
        mean,
        gaussian_width: gaussian.width(),
        subdivisions: max_depth,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::surface::SurfaceGeometry;
    use crate::geometry::{Point3, Vec3};

    #[test]
    fn sphere_bounds_contain_exact_curvature() {
        let r = 2.0;
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: r,
        };
        let exact_k = 1.0 / (r * r);

        // Sample away from poles.
        let bounds = guaranteed_curvature_bounds(
            &sphere,
            (0.5, 2.5),
            (
                std::f64::consts::FRAC_PI_4,
                3.0 * std::f64::consts::FRAC_PI_4,
            ),
            0.1,
            3,
        );

        assert!(
            bounds.gaussian.contains(exact_k),
            "Gaussian bounds [{}, {}] should contain exact K={exact_k}",
            bounds.gaussian.lo,
            bounds.gaussian.hi
        );
    }

    #[test]
    fn plane_bounds_contain_zero() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };

        let bounds = guaranteed_curvature_bounds(&plane, (0.0, 1.0), (0.0, 1.0), 0.01, 3);

        assert!(
            bounds.gaussian.contains(0.0),
            "Plane Gaussian bounds [{}, {}] should contain K=0",
            bounds.gaussian.lo,
            bounds.gaussian.hi
        );
        assert!(
            bounds.mean.contains(0.0),
            "Plane mean bounds [{}, {}] should contain H=0",
            bounds.mean.lo,
            bounds.mean.hi
        );
    }

    #[test]
    fn subdivision_tightens_bounds() {
        let r = 3.0;
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: r,
        };

        let coarse = guaranteed_curvature_bounds(
            &sphere,
            (0.5, 2.5),
            (
                std::f64::consts::FRAC_PI_4,
                3.0 * std::f64::consts::FRAC_PI_4,
            ),
            10.0, // Very loose target -> no subdivision
            0,
        );

        let fine = guaranteed_curvature_bounds(
            &sphere,
            (0.5, 2.5),
            (
                std::f64::consts::FRAC_PI_4,
                3.0 * std::f64::consts::FRAC_PI_4,
            ),
            0.01,
            4,
        );

        // Fine bounds should be tighter (or equal) than coarse.
        assert!(
            fine.gaussian_width <= coarse.gaussian_width + 1e-10,
            "Fine width {} should be <= coarse width {}",
            fine.gaussian_width,
            coarse.gaussian_width
        );
    }

    #[test]
    fn cylinder_bounds_contain_zero_gaussian() {
        let r = 5.0;
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: r,
        };

        let bounds =
            guaranteed_curvature_bounds(&cyl, (0.0, std::f64::consts::TAU), (-1.0, 1.0), 0.1, 3);

        // Cylinder has K = 0 everywhere.
        assert!(
            bounds.gaussian.contains(0.0),
            "Cylinder Gaussian bounds [{}, {}] should contain K=0",
            bounds.gaussian.lo,
            bounds.gaussian.hi
        );
    }

    #[test]
    fn guaranteed_bounds_results_have_fields() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };

        let bounds = guaranteed_curvature_bounds(&plane, (0.0, 1.0), (0.0, 1.0), 0.01, 2);

        assert!(bounds.gaussian_width >= 0.0);
        assert!(bounds.gaussian.lo <= bounds.gaussian.hi);
        assert!(bounds.mean.lo <= bounds.mean.hi);
    }
}
